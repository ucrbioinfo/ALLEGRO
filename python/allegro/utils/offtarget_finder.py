import os
import re
import sys
import pandas
import subprocess

from allegro.utils.iupac import expand_pam
from allegro.utils.shell_colors import bcolors

def reverse_complement(string: str) -> str:
    s = ''
    for c in string.upper():
        if c == 'A': s += 'T'
        elif c == 'C': s += 'G'
        elif c == 'G': s += 'C'
        elif c == 'T': s += 'A'
        else:
            print(f'{bcolors.RED}>{bcolors.RESET} Sequence {string} contains non-standard ({c}) base. Exiting.')
            sys.exit(1)
    return s[::-1]

def check_if_file_with_cached_index_was_modified(
    cache_index_dir: str,
    index_base_name: str,
    file_path: str) -> bool | str:
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
    _initialized: bool = False
    _base_path_cache: str = ''
    _cache_path_cache: str = ''
    _index_dir_cache: str = ''
    _reads_dir_cache: str = ''
    _alignments_dir_cache: str = ''
    _protospacer_length_cache: int = 0
    _expanded_pams_cache: list[str] = list()
    _pam_degenerate_offsets_cache: list[int] = list()

    def __init__(self) -> None:
        if not self.__class__._initialized:
            raise RuntimeError(
                'OfftargetFinder has not been initialized. '
                'Call OfftargetFinder.initialize(base_path, protospacer_length, pam) first.'
            )

    @property
    def base_path(self) -> str:
        return self.__class__._base_path_cache

    @property
    def cache_path(self) -> str:
        return self.__class__._cache_path_cache

    @property
    def index_dir(self) -> str:
        return self.__class__._index_dir_cache

    @property
    def reads_dir(self) -> str:
        return self.__class__._reads_dir_cache

    @property
    def alignments_dir(self) -> str:
        return self.__class__._alignments_dir_cache

    @property
    def protospacer_length(self) -> int:
        return self.__class__._protospacer_length_cache

    @property
    def expanded_pams(self) -> list[str]:
        return self.__class__._expanded_pams_cache

    @property
    def pam_degenerate_offsets(self) -> list[int]:
        return self.__class__._pam_degenerate_offsets_cache

    @classmethod
    def initialize(cls, base_path: str | None, protospacer_length: int, pam: str) -> None:
        if base_path is None:
            base_path = os.getcwd()

        base_path = os.path.abspath(base_path)

        if cls._initialized and cls._base_path_cache == base_path:
            return

        cache_path = os.path.join('allegro_cache', 'bowtie')

        index_dir = os.path.join(base_path, cache_path, 'index')
        reads_dir = os.path.join(base_path, cache_path, 'reads')
        alignments_dir = os.path.join(base_path, cache_path, 'alignments')

        os.makedirs(index_dir, exist_ok=True)
        os.makedirs(reads_dir, exist_ok=True)
        os.makedirs(alignments_dir, exist_ok=True)

        cls._pam_degenerate_offsets_cache = [
            idx for idx, c in enumerate(pam.upper())
            if c not in {'A', 'C', 'G', 'T'}
        ]

        cls._base_path_cache = base_path
        cls._cache_path_cache = cache_path
        cls._index_dir_cache = index_dir
        cls._reads_dir_cache = reads_dir
        cls._alignments_dir_cache = alignments_dir
        cls._expanded_pams_cache = expand_pam(pam)
        cls._protospacer_length_cache = protospacer_length
        cls._initialized = True
    
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
        """

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
                os.remove(os.path.join(self.index_dir, f"{index_base_name}{suffix}"))
            except FileNotFoundError:
                pass

    def write_guides_as_reads(self, species_name: str, guide_seqs_list: list[str]) -> None:
        """
        Write guide sequences as synthetic FASTQ reads for Bowtie alignment.
        """
        guides_w_pam = dict()

        for seq in guide_seqs_list:
            with_pam = [seq + pam for pam in self.expanded_pams]

            for wp in with_pam:
                guides_w_pam[wp] = seq

        reads_path = os.path.join(self.reads_dir, f'{species_name}_reads.fq')
        with open(reads_path, 'w') as f:
            for idx, guide in enumerate(guides_w_pam.keys()):
                I_str = 'I' * len(guide)
                f.write(f'@READ_{idx+1}\n{guide}\n+\n{I_str}\n')

    def run_bowtie_build(self, species_name: str, path_to_background_fasta: str, gene_or_genome_str: str) -> str:
        """
        Build (or reuse) a cached Bowtie index for a background FASTA.

        Returns:
            str: The Bowtie index basename (full path without .ebwt suffix).
        """
        index_base_name = f"{species_name}_{gene_or_genome_str}_idx"

        index_dir = os.path.join(self.base_path, self.cache_path, "index")
        idx_base = os.path.join(index_dir, index_base_name)
        fasta_path = os.path.join(self.base_path, path_to_background_fasta)

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

        idx_base = os.path.join(self.index_dir, f"{that_species_name}_{gene_or_genome_str}_idx")
        reads_fq = os.path.join(self.reads_dir, f"{this_species_name}_reads.fq")
        alignments_bowtie = os.path.join(self.alignments_dir, f"{this_species_name}_against_{that_species_name}_{mm}mm_alignment.sam")

        bowtie_command = [
            "bowtie",
            "-a",
            "-v", mm,
            "--quiet",
            "-x", idx_base,
            reads_fq,
            alignments_bowtie,
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

        try:
            df_mm_genomic = pandas.read_csv(
                alignments_bowtie,
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
        except pandas.errors.EmptyDataError:
            try:
                os.remove(alignments_sam)
            except FileNotFoundError:
                pass
            return pandas.DataFrame(columns=[
                "query_name",
                "strand",
                "reference_name",
                "start_position",
                "aligned_seq",
                "mapping_quality",
                "idk",
                "mismatch",
            ])

        if df_mm_genomic.empty:
            # Clean up temp file if it exists.
            try:
                os.remove(alignments_bowtie)
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

        df_mm_genomic["pam"] = [s[self.protospacer_length:] for s in guide_w_pam]
        df_mm_genomic["sequence"] = [s[:self.protospacer_length] for s in guide_w_pam]

        # Remove mismatches occurring at the first PAM base.
        df_mm_genomic["mismatch"] = df_mm_genomic["mismatch"].fillna("N/A")
        pam_degenerate_positions = {
            self.protospacer_length + offset
            for offset in self.pam_degenerate_offsets
        }

        def has_only_allowed_pam_mismatches(mismatch_str: str) -> bool:
            if mismatch_str == "N/A":
                return True

            mismatch_positions = [int(num) for num in re.findall(r'\d+', mismatch_str)]

            for pos in mismatch_positions:
                if pos >= self.protospacer_length and pos not in pam_degenerate_positions:
                    return False

            return True

        df_mm_genomic = df_mm_genomic[
            df_mm_genomic["mismatch"].apply(has_only_allowed_pam_mismatches)
        ]

        try:
            os.remove(alignments_bowtie)
        except FileNotFoundError:
            pass

        return df_mm_genomic.reset_index(drop=True)
        