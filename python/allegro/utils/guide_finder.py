import re
import os
import sys
import time
import pandas
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
from itertools import product

from allegro.scorers.scorer_base import Scorer
from allegro.utils.shell_colors import bcolors
from allegro.utils.iupac_dict import iupac_dict

def build_pam_regex(pam: str) -> str:
    pattern = ''

    for c in pam.upper():
        if c not in iupac_dict:
            raise ValueError(f'{bcolors.RED}>{bcolors.RESET} Dev error. Unknown IUPAC code: {c}')

        bases = ''.join(iupac_dict[c])
        pattern += f'[{bases}]'

    return rf'(?=({pattern}))'

def calculate_gc_content(sequence):
    """
    Compute GC content of a nucleotide sequence.

    GC content is defined as:

        (count('G') + count('C')) / len(sequence)

    The input is uppercased for counting.

    Args:
        sequence (str): Nucleotide sequence (DNA/RNA). Must be non-empty.

    Returns:
        float: Fraction in [0, 1] representing GC content.

    Raises:
        ZeroDivisionError: If `sequence` is empty.
    """
    
    return (sequence.upper().count('G') + sequence.upper().count('C')) / len(sequence)

def count_kmers(sequence, k):
    """
    Count all k-mers in a sequence using a sliding window.

    For each position i in [0, len(sequence) - k], the substring sequence[i:i+k]
    is counted.

    Args:
        sequence (str): Input sequence.
        k (int): k-mer length (must be >= 1).

    Returns:
        dict[str, int]: Mapping from k-mer string to its count.

    Notes:
        * If k > len(sequence), returns an empty dict.
        * This function does not uppercase the sequence; k-mers preserve the
          original casing.
    """
    kmers = dict()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers[kmer] = kmers.setdefault(kmer, 0) + 1

    return kmers

def align_guides_to_seq_bowtie(
    name: str,
    output_align_temp_path: str,
    file_path: str,
    output_directory: str,
    ) -> tuple[list[str], list[str], list[str], list[str], list[str], float]:
    """
    Align guide sequences to a FASTA reference using Bowtie (v1) and return hits.

    This function:
      1) Builds a Bowtie index from the provided FASTA file via `bowtie-build`.
      2) Aligns all query sequences (FASTA) in `output_align_temp_path` to the
         index via `bowtie` with exact matching (`-v 0`) and reports all hits (`-a`).
      3) Parses Bowtie tabular output into per-hit fields.
      4) Reads the reference FASTA and extracts `orthologous_to_gene` annotations
         from record descriptions when present.

    Args:
        name (str): Base name used to create the Bowtie index basename
            `{output_directory}/{name}_idx`.
        output_align_temp_path (str): Path to a FASTA file containing query
            guide sequences to align. This path is interpreted relative to the
            current working directory (because `os.getcwd()` is prepended).
        file_path (str): Path to the reference FASTA used to build the Bowtie
            index. Also used for extracting ortholog annotations.
        output_directory (str): Directory where the Bowtie index files are written.

    Returns:
        tuple[
            list[str],  # sequences (query names as emitted by Bowtie)
            list[str],  # strands ('+' or '-')
            list[str],  # reference_names (FASTA record IDs)
            list[str],  # orthos (parsed [orthologous_to_gene=...] or 'N/A')
            list[str],  # start_positions (0-based start coordinate from Bowtie)
            float       # elapsed_seconds (CPU process time)
        ]: Per-hit fields and runtime.

    Exits: 
        sys.exit(1) on non-zero Bowtie exit code

    Side Effects:
        Prints errors and exits the process with `sys.exit(1)` if Bowtie
        exits with a return code not 0.
        Writes Bowtie index files into `output_directory`.

    Requirements:
        * `bowtie-build` and `bowtie` must be available on PATH.
        * `pandas` and Biopython (`Bio.SeqIO`) must be installed.
    """

    base_path = os.getcwd()
    ref_fa = os.path.join(base_path, file_path)
    index_basename = os.path.join(output_directory, f'{name}_idx')
    idx_base = os.path.join(base_path, index_basename)
    query_fa = os.path.join(base_path, output_align_temp_path)

    bowtie_build_command = [
        'bowtie-build',
        '--quiet',
        '-f', 
        ref_fa, idx_base
    ]

    bowtie_command = [
        "bowtie",
        "-a",
        "-v", "0",
        "--quiet",
        "--suppress", "5,6,7,8",
        "-f",
        "-x", idx_base, query_fa
    ]

    start_time = time.process_time()

    # BUILD IDX
    process = subprocess.Popen(bowtie_build_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_build, stderr_build = process.communicate()

    if process.returncode != 0:
        print(f"{bcolors.RED}> Error{bcolors.RESET}: Dev error. bowtie-build failed (exit {process.returncode}).")
        if stderr_build:
            print(stderr_build.decode(errors="replace"))
        sys.exit(1)

    # ALIGN
    process = subprocess.Popen(bowtie_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"{bcolors.RED}> Error{bcolors.RESET}: Dev error. bowtie failed (exit {process.returncode}).")
        if stderr:
            print(stderr.decode(errors="replace"))
        sys.exit(1)

    end_time = time.process_time()
    # Calculate the elapsed time in seconds
    elapsed_seconds = end_time - start_time

    text = stdout.decode(errors="replace").strip()

    if text == "":
        return [], [], [], [], [], elapsed_seconds

    df = pandas.read_csv(
        StringIO(text),
        sep='\t',
        names=['query_name', 'strand', 'reference_name', 'start_position'],
        dtype={'reference_name': str, 'query_name': str}
    )

    records = list(SeqIO.parse(ref_fa, "fasta"))

    map_ref_id_to_ortho = dict()
    for record in records:
        match = re.search(r'\[orthologous_to_gene=(.*?)\]', record.description)
        map_ref_id_to_ortho[record.id] = match.group(1) if match else 'N/A'

    sequences = df['query_name'].tolist()
    strands = df['strand'].tolist()
    reference_names = df['reference_name'].tolist()
    orthos = [map_ref_id_to_ortho[ref_id] for ref_id in reference_names]
    start_positions = df['start_position'].tolist()

    return sequences, strands, reference_names, orthos, start_positions, elapsed_seconds

class GuideFinder:
    _initialized: bool = False
    _filter_by_gc_cache: bool = False
    _gc_min_cache: float = 0.0
    _gc_max_cache: float = 0.0
    _guide_score_threshold_cache = 0.0
    _protospacer_length_cache: int = 0
    _pam_dict_cache: dict[str, re.Pattern] = dict()
    _exclusion_list_cache: list[str] = list()
    _exclusion_patterns_combinations_cache: list[str] = list()
    _guide_scorer_cache: Scorer

    def __init__(self) -> None:
        if not self.__class__._initialized:
            raise RuntimeError(
                'GuideFinder has not been initialized. ' + \
                'Call GuideFinder.initialize(blocklist_path' + \
                ', pam_list' + \
                ', filter_by_gc' + \
                ', gc_min' + \
                ', gc_max' + \
                ', protospacer_length' + \
                ', patterns_to_exclude' + \
                ', guide_scorer) first.')

    @property
    def pam_dict(self) -> dict[str, re.Pattern]:
        return self.__class__._pam_dict_cache

    @property
    def protospacer_length(self) -> int:
        return self.__class__._protospacer_length_cache

    @property
    def filter_by_gc(self) -> bool:
        return self.__class__._filter_by_gc_cache

    @property
    def gc_min(self) -> float:
        return self.__class__._gc_min_cache

    @property
    def gc_max(self) -> float:
        return self.__class__._gc_max_cache

    @property
    def exclusion_list(self) -> list[str]:
        return self.__class__._exclusion_list_cache

    @property
    def exclusion_patterns_combinations(self) -> list[str]:
        return self.__class__._exclusion_patterns_combinations_cache

    @property
    def guide_scorer(self) -> Scorer:
        return self.__class__._guide_scorer_cache

    @property
    def guide_score_threshold(self) -> float:
        return self.__class__._guide_score_threshold_cache

    @classmethod
    def initialize(
        cls,
        blocklist_path: str,
        pam_list: list[str],
        filter_by_gc: bool,
        gc_min: float,
        gc_max: float,
        protospacer_length: int,
        patterns_to_exclude: list[str],
        guide_scorer: Scorer,
        guide_score_threshold: float) -> None:
        normalized_pams = [pam.strip().upper() for pam in pam_list if pam.strip()]
        normalized_patterns = [pattern.strip().upper() for pattern in patterns_to_exclude if pattern.strip()]

        cls._pam_dict_cache = {
            pam: re.compile(build_pam_regex(pam))
            for pam in normalized_pams
        }

        if os.path.exists(blocklist_path):
            with open(blocklist_path, 'r') as file:
                cls._exclusion_list_cache = file.read().splitlines()
            if cls._exclusion_list_cache:
                print(f'{bcolors.BLUE}>{bcolors.RESET} Excluding the guides in {bcolors.BLACK}{blocklist_path}{bcolors.RESET}')
        else:
            cls._exclusion_list_cache = list()

        combinations: list[str] = list()
        for pattern in normalized_patterns:
            for code in pattern:
                if code not in iupac_dict:
                    raise ValueError(f'Unknown IUPAC code in exclusion pattern: {code}')

            lists = [iupac_dict[code] for code in pattern]
            combinations.extend(''.join(pair) for pair in product(*lists))

        cls._protospacer_length_cache = protospacer_length
        cls._filter_by_gc_cache = filter_by_gc
        cls._gc_min_cache = gc_min
        cls._gc_max_cache = gc_max
        cls._exclusion_patterns_combinations_cache = list(dict.fromkeys(combinations))
        cls._guide_scorer_cache = guide_scorer
        cls._guide_score_threshold_cache = guide_score_threshold
        cls._initialized = True

    def contains_iupac_pattern(self, sequence: str, patterns: list[str]) -> bool:
        sequence = sequence.strip().upper()

        for pattern in patterns:
            pattern = pattern.strip().upper()

            if len(pattern) == len(sequence):
                if pattern == sequence:
                    return True

            elif len(pattern) < len(sequence):
                # For example, if pattern = 'MK', if any of the ['GA', 'GC', 'TA', 'TC'] 
                # appear in sequence, return True
                if pattern in sequence:
                    return True

            else:
                print(f'{bcolors.ORANGE}>{bcolors.RESET} Pattern {pattern} is longer than {sequence}. Ignoring this pattern for this sequence.')
                continue

        return False

    def score_guides_from_list(self, guides_context_list: list[str]) -> list[float]:
        return self.guide_scorer.score_sequence(guides_context_list)

    def identify_guides_and_indicate_strand(
        self, 
        sequence: str) -> tuple[list[str], list[str], list[str], list[int], list[float]]:
        guides_list: list[str] = list()
        guides_context_list: list[str] = list()
        strands_list: list[str] = list()
        locations_list: list[int] = list()
        scores_list: list[float] = list()

        valid_bases = {'A', 'C', 'G', 'T'}
        
        for strand in ['+', '-']:
            seq = sequence if strand == '+' else str(Seq(sequence).reverse_complement())

            for pam_name, pam_regex in self.pam_dict.items():
                for match in pam_regex.finditer(seq):
                    position = match.start()
                    if position - self.protospacer_length < 0:
                        continue

                    guide = seq[position-self.protospacer_length:position]

                    # Skip guides with non-standard nucleotides, or the gap delimiter: | 
                    if any(c not in valid_bases for c in guide):
                        continue

                    if guide in self.exclusion_list:
                        continue

                    if self.exclusion_patterns_combinations:
                        if self.contains_iupac_pattern(guide, self.exclusion_patterns_combinations):
                            continue

                    if self.filter_by_gc:
                        gc = calculate_gc_content(guide)
                        if (gc > self.gc_max) or (gc < self.gc_min):
                            continue

                    # There is enough context on both sides of the PAM
                    guide_with_context = ''
                    if (
                        position - self.protospacer_length - self.guide_scorer.context_toward_five_prime >= 0 
                        and position + len(pam_name) + self.guide_scorer.context_toward_three_prime < len(seq) + 1
                    ):
                        guide_with_context = seq[
                            position - self.protospacer_length - self.guide_scorer.context_toward_five_prime:
                            position + len(pam_name) + self.guide_scorer.context_toward_three_prime
                        ]

                        # If we're not using uCRISPR, or a scorer that doesn't need a context, context_toward_five_prime and 
                        # context_toward_three_prime should be set to 0 by the configurator in settings. That way, 
                        # guides with | around them (but not in them) won't be discarded.
                        if '|' in guide_with_context:
                            continue
                    else:
                        continue

                    guide_score = self.score_guides_from_list([guide_with_context])

                    if guide_score[0] < self.guide_score_threshold:
                        continue

                    guides_list.append(guide)
                    guides_context_list.append(guide_with_context)
                    strands_list.append(strand)
                    locations_list.append(position)
                    scores_list.append(guide_score[0])

        return guides_list, guides_context_list, strands_list, locations_list, scores_list
    
    def locate_guides_in_sequence(
        self,
        sequence: str,
        file_path: str,
        to_upper: bool=True) -> list[tuple[str, str, int, int, str, str]]:

        chrom_strand_start_end_misc: list[tuple[str, str, int, int, str, str]] = list()

        seq = sequence.upper() if to_upper else sequence

        with open(file_path, 'r') as handle:
            records = list(SeqIO.parse(handle, 'fasta'))
        
        for strand in ['+', '-']:
            for record in records:
                dna = str(record.seq).upper() if to_upper else str(record.seq)
                if strand == "-":
                    dna = str(record.seq.reverse_complement()).upper() if to_upper else str(record.seq.reverse_complement())

                ortho_to = "N/A"
                match = re.search(r"\[orthologous_to_gene=(.*?)\]", record.description)
                if match:
                    ortho_to = match.group(1)

                misc = 'N/A'
                parts = record.description.split()
                if len(parts) > 1:
                    misc = parts[1]

                for pam_name in self.pam_dict.keys():
                    seq_and_pam = re.compile(re.escape(seq) + build_pam_regex(pam_name))

                    for match in seq_and_pam.finditer(dna):
                        match_start = match.start()
                        match_end = match.end()

                        if strand == "-":
                            match_start = len(record.seq) - match.start() - len(seq)
                            match_end = len(record.seq) - match.end() + len(seq)

                        chrom_strand_start_end_misc.append(
                            (record.id, strand, match_start, match_end, ortho_to, misc)
                        )

        return chrom_strand_start_end_misc
