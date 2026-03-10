import os
import sys
import pandas
import subprocess

from allegro.utils.shell_colors import bcolors

def reverse_complement(string):
    """
    Compute the reverse complement of a DNA sequence.

    This function:
      * Uppercases the input.
      * Uses the DNA complement rules A<->T and C<->G.
      * Reverses the complemented string.

    Args:
        string (str): DNA sequence consisting of characters A/C/G/T (case-insensitive).

    Returns:
        str: Reverse-complemented sequence.
    """
    
    s = ''
    for c in string.upper():
        if c == 'A': s += 'T'
        elif c == 'C': s += 'G'
        elif c == 'G': s += 'C'
        elif c == 'T': s += 'A'
        else:
            print(f'{bcolors.RED}>{bcolors.RESET} Sequence {string} contains non-standard ({c}) base.')
            print(f'{bcolors.RED}>{bcolors.RESET} Exiting.')
            sys.exit(1)
    return s[::-1]

def check_if_file_with_cached_index_was_modified(
    cache_index_dir: str,
    index_base_name: str,
    file_path: str
) -> bool | str:
    """
    Determine whether a FASTA file has changed since a cached Bowtie index was created.
    """

    EPSILON = 1e-6  # seconds tolerance for float comparison

    if os.path.exists(file_path):
        creation_record_file_path = os.path.join(
            cache_index_dir, f"{index_base_name}_last_modified_date.txt"
        )

        if os.path.exists(creation_record_file_path):
            with open(creation_record_file_path, "r") as record_file:
                recorded_creation_date = float(record_file.readline().strip())

            last_modified_time = os.path.getmtime(file_path)

            if abs(last_modified_time - recorded_creation_date) > EPSILON:
                return True
        else:
            return "N/A"

    return False

def record_creation_date(cache_index: str, index_base_name: str, file_path: str) -> None:
    """
    Record the current modification time of a FASTA file for Bowtie index caching.

    Writes the file's modification timestamp (mtime) to:

        {cache_index}/{index_base_name}_last_modified_date.txt

    This timestamp is later used by `check_if_file_with_cached_index_was_modified()`
    to decide whether a cached Bowtie index should be invalidated and rebuilt.

    Args:
        cache_index (str): Directory where the record file should be written
            (typically the Bowtie index cache directory).
        index_base_name (str): Base name used to create the record filename.
        file_path (str): Path to the FASTA file whose mtime should be recorded.

    Returns:
        None

    Raises:
        FileNotFoundError: If `file_path` does not exist.
        OSError: If the record file cannot be written.
    """
    creation_record_file_path = os.path.join(cache_index, f'{index_base_name}_last_modified_date.txt')

    # Get the modified time of the file
    modified_time = os.path.getmtime(file_path)

    with open(creation_record_file_path, 'w') as record_file:
        record_file.write(str(modified_time))

class OfftargetFinder:
    _self = None
    _base_path = None
    _cache_path = None

    # Singleton
    def __new__(self):
        """
        Create or return the singleton instance.

        OfftargetFinder is implemented as a singleton so that:
          * The working directory base path is captured once.
          * Bowtie cache directories are created once and reused across calls.

        Returns:
            OfftargetFinder: The singleton instance.
        """
        if self._self is None:
            self._self = super().__new__(self)
        return self._self
    
    def __init__(self) -> None:
        """
        Initialize cache paths and ensure Bowtie cache directories exist.

        Sets:
          * `_base_path` to the current working directory (os.getcwd()).
          * `_cache_path` to `allegro_cache/bowtie`.

        Creates (if missing) the following directories under `_cache_path`:
          * index/       cached Bowtie index files (.ebwt)
          * reads/       generated FASTQ "reads" from guide sequences
          * alignments/  temporary Bowtie alignment outputs (.sam)

        Side Effects:
            Creates directories on disk.
        """
        if self._base_path is None:
            self._base_path = os.getcwd()

        if self._cache_path is None:
            self._cache_path = os.path.join("allegro_cache", "bowtie")

        index_dir = os.path.join(self._base_path, self._cache_path, "index")
        reads_dir = os.path.join(self._base_path, self._cache_path, "reads")
        alignments_dir = os.path.join(self._base_path, self._cache_path, "alignments")

        self._index_dir = index_dir
        self._reads_dir = reads_dir
        self._alignments_dir = alignments_dir

        os.makedirs(index_dir, exist_ok=True)
        os.makedirs(reads_dir, exist_ok=True)
        os.makedirs(alignments_dir, exist_ok=True)
    
    def purge_cached_index(self, index_base_name: str) -> None:
        """
        Delete cached Bowtie index files for a given index basename.

        Removes the standard Bowtie v1 index files (if present):

            {index_base_name}.1.ebwt
            {index_base_name}.2.ebwt
            {index_base_name}.3.ebwt
            {index_base_name}.4.ebwt
            {index_base_name}.rev.1.ebwt
            {index_base_name}.rev.2.ebwt

        Args:
            index_base_name (str): Index basename used when building the Bowtie index.

        Returns:
            None
        """
        index_dir = os.path.join(self._base_path, self._cache_path, "index")

        suffixes = [
            ".1.ebwt",
            ".2.ebwt",
            ".3.ebwt",
            ".4.ebwt",
            ".rev.1.ebwt",
            ".rev.2.ebwt",
        ]

        for suffix in suffixes:
            try:
                os.remove(os.path.join(index_dir, f"{index_base_name}{suffix}"))
            except FileNotFoundError:
                pass

    def write_guides_as_reads(self, species_name: str, guide_seqs_list: list[str]) -> None:
        """
        Write guide sequences as synthetic FASTQ reads for Bowtie alignment.

        For each protospacer in `guide_seqs_list`, this function generates
        4 sequences by appending each possible 'NGG' PAM instantiation:

            AGG, CGG, TGG, GGG

        These sequences are then written as FASTQ reads to:

            {cache_path}/reads/{species_name}_reads.fq

        Each read is assigned a dummy quality string of 'I' characters.

        Args:
            species_name (str): Species name used to name the FASTQ file.
            guide_seqs_list (list[str]): List of protospacer sequences (typically length 20).

        Returns:
            None

        Side Effects:
            Writes/overwrites a FASTQ file on disk.
        """
        guides_w_pam = dict()

        for seq in guide_seqs_list:
            with_pam = [seq + pam for pam in ['AGG', 'CGG', 'TGG', 'GGG']]

            for wp in with_pam:
                guides_w_pam[wp] = seq

        reads_path = f'{self._cache_path}/reads/{species_name}_reads.fq'
        with open(reads_path, 'w') as f:
            for idx, guide in enumerate(guides_w_pam.keys()):
                f.write(f'@READ_{idx+1}\n{guide}\n+\nIIIIIIIIIIIIIIIIIIIIIII\n')

    def run_bowtie_build(self, species_name: str, path_to_background_fasta: str, gene_or_genome_str: str) -> str:
        """
        Build (or reuse) a cached Bowtie index for a background FASTA.

        Returns:
            str: The Bowtie index basename (full path without .ebwt suffix).
        """
        index_base_name = f"{species_name}_{gene_or_genome_str}_idx"

        index_dir = os.path.join(self._base_path, self._cache_path, "index")
        idx_base = os.path.join(index_dir, index_base_name)
        fasta_path = os.path.join(self._base_path, path_to_background_fasta)

        one_exists = os.path.exists(f"{idx_base}.1.ebwt")
        two_exists = os.path.exists(f"{idx_base}.2.ebwt")
        three_exists = os.path.exists(f"{idx_base}.3.ebwt")
        four_exists = os.path.exists(f"{idx_base}.4.ebwt")
        rev_one_exists = os.path.exists(f"{idx_base}.rev.1.ebwt")
        rev_two_exists = os.path.exists(f"{idx_base}.rev.2.ebwt")

        if one_exists and two_exists and three_exists and four_exists and rev_one_exists and rev_two_exists:
            modified = check_if_file_with_cached_index_was_modified(index_dir, index_base_name, fasta_path)

            if modified:
                if modified != "N/A":
                    print(
                        f"{bcolors.BLUE}>{bcolors.RESET} {fasta_path} was modified since it was last cached. Rebuilding index."
                    )
                self.purge_cached_index(index_base_name)
            else:
                # Use cached indices.
                return idx_base

        # Rebuild indices.
        bowtie_build_command = ["bowtie-build", "--quiet", "-f", fasta_path, idx_base]

        process = subprocess.Popen(bowtie_build_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_build, stderr_build = process.communicate()

        if process.returncode != 0:
            print(
                f"{bcolors.RED}> Error{bcolors.RESET}: bowtie-build failed (exit {process.returncode})."
            )
            if stderr_build:
                print(stderr_build.decode(errors="replace"))
            if stdout_build:
                print(stdout_build.decode(errors="replace"))
            sys.exit(1)

        record_creation_date(index_dir, index_base_name, fasta_path)
        return idx_base

    def run_bowtie_against_other(
        self,
        this_species_name: str,
        that_species_name: str,
        gene_or_genome_str: str,
        num_of_mismatches: int
    ) -> pandas.DataFrame:
        """
        Align one species' guide "reads" against another species' Bowtie index.

        Returns:
            pandas.DataFrame: alignments (possibly empty).
        """
        mm = str(int(num_of_mismatches))

        idx_base = os.path.join(self._index_dir, f"{that_species_name}_{gene_or_genome_str}_idx")
        reads_fq = os.path.join(self._reads_dir, f"{this_species_name}_reads.fq")
        alignments_sam = os.path.join(self._alignments_dir, f"{this_species_name}_against_{that_species_name}_{mm}mm_alignment.sam")

        bowtie_command = [
            "bowtie",
            "-a",
            "-v", mm,
            "--quiet",
            "-x", idx_base,
            reads_fq,
            alignments_sam,
        ]

        process = subprocess.Popen(bowtie_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_bt, stderr_bt = process.communicate()

        if process.returncode != 0:
            print(f"{bcolors.RED}> Error{bcolors.RESET}: bowtie failed (exit {process.returncode}).")
            if stderr_bt:
                print(stderr_bt.decode(errors="replace"))
            if stdout_bt:
                print(stdout_bt.decode(errors="replace"))
            sys.exit(1)

        df_mm_genomic = pandas.read_csv(
            alignments_sam,
            sep="\t",
            names=[
                "query_name",
                "strand",
                "reference_name",
                "start_position",
                "aligned_seq",
                "mapping_quality",
                "idk",
                "mismatch",
            ],
        )

        if len(df_mm_genomic) == 0:
            # Clean up temp file if it exists.
            try:
                os.remove(aln_sam)
            except FileNotFoundError:
                pass
            return df_mm_genomic

        guide_w_pam: list[str] = []
        for _, row in df_mm_genomic.iterrows():
            aligned = row["aligned_seq"]
            if row["strand"] == "-":
                guide_w_pam.append(reverse_complement(aligned))
            else:
                guide_w_pam.append(aligned)

        df_mm_genomic["pam"] = [s[-3:] for s in guide_w_pam]
        df_mm_genomic["sequence"] = [s[:-3] for s in guide_w_pam]

        # Remove mismatches occurring at the N base of the NGG PAM
        df_mm_genomic["mismatch"] = df_mm_genomic["mismatch"].fillna("N/A")
        df_mm_genomic = df_mm_genomic[~df_mm_genomic["mismatch"].str.contains("20:", na=False)]
        df_mm_genomic.drop(columns=["idk", "mapping_quality"], inplace=True)

        try:
            os.remove(alignments_sam)
        except FileNotFoundError:
            pass

        return df_mm_genomic.reset_index(drop=True)
        