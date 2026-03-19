import os
import re
import pandas
from concurrent.futures import ProcessPoolExecutor, as_completed

from allegro.utils.shell_colors import bcolors
from allegro.utils.offtarget_finder import OfftargetFinder

_WORKER_EXPERIMENT_NAME = None
_WORKER_SPECIES_TO_FILE = None
_WORKER_OUTPUT_LIBRARY = None
_WORKER_GENE_OR_GENOME_STR = None
_WORKER_INPUT_SPECIES_OFFTARGET_DIR = None
_WORKER_NUM_MISMATCHES = None
_WORKER_SEED_REGION_IS_N_UPSTREAM_OF_PAM = None

def _init_offtarget_worker(
    output_directory: str,
    species_to_file: dict[str, str],
    output_library: pandas.DataFrame,
    experiment_name: str,
    gene_or_genome_str: str,
    input_species_offtarget_dir: str,
    num_mismatches: int,
    seed_region_is_n_upstream_of_pam: int,
    protospacer_length: int,
    pam: str,
) -> None:
    global _WORKER_SPECIES_TO_FILE
    global _WORKER_OUTPUT_LIBRARY
    global _WORKER_GENE_OR_GENOME_STR
    global _WORKER_INPUT_SPECIES_OFFTARGET_DIR
    global _WORKER_NUM_MISMATCHES
    global _WORKER_SEED_REGION_IS_N_UPSTREAM_OF_PAM
    global _WORKER_EXPERIMENT_NAME

    OfftargetFinder.initialize(None, protospacer_length, pam)

    _WORKER_EXPERIMENT_NAME = experiment_name
    _WORKER_SPECIES_TO_FILE = species_to_file
    _WORKER_OUTPUT_LIBRARY = output_library
    _WORKER_GENE_OR_GENOME_STR = gene_or_genome_str
    _WORKER_INPUT_SPECIES_OFFTARGET_DIR = input_species_offtarget_dir
    _WORKER_NUM_MISMATCHES = num_mismatches
    _WORKER_SEED_REGION_IS_N_UPSTREAM_OF_PAM = seed_region_is_n_upstream_of_pam
    
def _find_targets_worker(target_species: str):
    return find_targets(
        target_species=target_species,
        experiment_name=_WORKER_EXPERIMENT_NAME,
        species_to_file=_WORKER_SPECIES_TO_FILE,
        output_library=_WORKER_OUTPUT_LIBRARY,
        gene_or_genome_str=_WORKER_GENE_OR_GENOME_STR,
        input_species_offtarget_dir=_WORKER_INPUT_SPECIES_OFFTARGET_DIR,
        num_mismatches=_WORKER_NUM_MISMATCHES,
        seed_region_is_n_upstream_of_pam=_WORKER_SEED_REGION_IS_N_UPSTREAM_OF_PAM,
    )

def hamming_distance(str1: str, str2: str, length: int) -> int:
    """
    Calculates the Hamming distance between two strings up to a specified length.
    """
    if len(str1) < length or len(str2) < length:
        raise ValueError('Strings must have a length of at least {n}'.format(n=length))

    return sum(ch1 != ch2 for ch1, ch2 in zip(str1[:length], str2[:length]))

def cluster_strings(strings: list[str], req_match_len: int, mm_allowed: int) -> list[list[str]]:
    """
    Clusters a list of strings based on a mismatch in the first `mm_allowed` letters.
    """
    clusters = list()
    for string in strings:
        # Check if the string belongs to any existing cluster.
        found_cluster = False

        for cluster in clusters:
            agreement = [False for _ in range(len(cluster))]

            for idx, existing_string in enumerate(cluster):
                if string[req_match_len:] == existing_string[req_match_len:] and hamming_distance(string, existing_string, req_match_len) <= mm_allowed:
                    agreement[idx] = True
                else:
                    break

            if all(agreement):
                cluster.append(string)
                found_cluster = True
                break

        if not found_cluster:
            # Create a new cluster for the string.
            clusters.append([string])

    return clusters

def cluster_solution(experiment_name: str,
                     output_directory: str,
                     req_match_len: int,
                     mm_allowed: int) -> int:
        
    print(f'{bcolors.BLUE}>{bcolors.RESET} Clustering guides: Guides in the same cluster have an identical sequence for the first {req_match_len} nucleotides after the PAM (3\' to 5\').')
    print(f'{bcolors.BLUE}>{bcolors.RESET} Guides in the same cluster may mismatch up to {mm_allowed} nucleotides after the seed region.')
    
    solution_path = os.path.join(output_directory, experiment_name + '.csv')

    df = pandas.read_csv(solution_path)
    df['cluster'] = 0
    seqs = df.sequence.unique().tolist()

    clusters = cluster_strings(seqs, req_match_len, mm_allowed)

    for idx, cluster in enumerate(clusters):
        df.loc[df['sequence'].isin(cluster), 'cluster'] = idx

    df.to_csv(solution_path, index=False)

    print(f'{bcolors.BLUE}>{bcolors.RESET} Done clustering. Added a new column \'cluster\' to {solution_path}')
    print(f'{bcolors.BLUE}>{bcolors.RESET} The output guide RNA set contains {len(clusters)} clusters.')

    return len(clusters)

def extract_complete_digits(s: str) -> list[int]:
        return [int(num) for num in re.findall(r'\d+', s)]

def determine_ot_in_nonseed(
    digits: list[int],
    guide_len: int,
    seed_region_is_n_upstream_of_pam: int) -> bool:
        return any(digit + 1 + seed_region_is_n_upstream_of_pam > guide_len for digit in digits)

# Function to use with multiprocessing
def find_targets(
    target_species: str,
    experiment_name: str,
    species_to_file: dict[str, str],
    output_library: pandas.DataFrame,
    gene_or_genome_str: str,
    input_species_offtarget_dir: str,
    num_mismatches: int,
    seed_region_is_n_upstream_of_pam: int) -> pandas.DataFrame:

    target_species_offtarget_file = os.path.join(
        input_species_offtarget_dir,
        species_to_file[target_species]
    )

    OTF = OfftargetFinder()
    OTF.run_bowtie_build(target_species, target_species_offtarget_file, gene_or_genome_str)

    all_targets = OTF.run_bowtie_against_other(
        f'{experiment_name}_output_lib',
        target_species,
        gene_or_genome_str,
        num_mismatches
    )

    if all_targets.empty:
        return all_targets

    all_targets['reference_name'] = all_targets['reference_name'].astype(str)
    all_targets['target_species'] = target_species
    all_targets['is_off_target'] = '1'
    all_targets['self_off_targets'] = '0'

    on_targets = pandas.merge(
        output_library[output_library['target'] == target_species],
        all_targets,
        on=['sequence', 'reference_name', 'strand', 'start_position']
    )

    cols = [
        'is_off_target', 'self_off_targets', 'sequence', 'pam', 'mismatch', 'aligned_seq',
        'target_species', 'strand', 'reference_name', 'orthologous_to', 'start_position'
    ]

    on_targets = on_targets[cols].copy()
    on_targets['is_off_target'] = '0'

    merged_df = all_targets.merge(
        on_targets,
        on=['sequence', 'reference_name', 'strand', 'start_position'],
        how='left',
        suffixes=('', '_new')
    )

    all_targets['orthologous_to'] = 'N/A'
    all_targets['is_off_target'] = merged_df['is_off_target_new'].combine_first(all_targets['is_off_target'])
    all_targets['orthologous_to'] = merged_df['orthologous_to'].combine_first(all_targets['orthologous_to'])
    all_targets = all_targets[cols]

    all_targets['digits'] = all_targets['mismatch'].apply(extract_complete_digits)
    all_targets['ot_in_seed'] = all_targets.apply(
        lambda row: determine_ot_in_nonseed(
            row['digits'],
            len(row['sequence']),
            seed_region_is_n_upstream_of_pam
        ),
        axis=1
    )

    all_targets = all_targets[~all_targets['ot_in_seed']]
    all_targets = all_targets.drop(columns=['digits', 'ot_in_seed'])

    OT = all_targets[all_targets['is_off_target'] == '1']
    for seq in OT['sequence'].unique():
        seq_off_target_species = set(OT[OT['sequence'] == seq]['target_species'])
        all_target_species = set(all_targets[all_targets['sequence'] == seq]['target_species'])

        if seq_off_target_species.intersection(all_target_species):
            all_targets.loc[all_targets['sequence'] == seq, 'self_off_targets'] = '1'

    return all_targets

def report_offtargets(
    input_species_path: str,
    output_directory: str,
    experiment_name: str,
    input_species_offtarget_dir: str,
    input_species_offtarget_column: str,
    num_mismatches: int,
    seed_region_is_n_upstream_of_pam: int,
    protospacer_length: int,
    pam: str) -> None:

    pam_length = len(pam)
    solution_path = os.path.join(output_directory, experiment_name + '.csv')
    gene_or_genome_str = 'genome' if 'genome' in input_species_offtarget_column else 'genes'
    
    OfftargetFinder.initialize(None, protospacer_length, pam)
    OTF = OfftargetFinder()

    created_dfs: list[pandas.DataFrame] = list()
    species_df = pandas.read_csv(input_species_path)
    output_library = pandas.read_csv(solution_path)

    output_library.loc[output_library['strand'] == '-', 'start_position'] -= pam_length

    output_library = output_library.rename(columns={'sequence': 'seq_with_pam'}).copy()
    output_library['reference_name'] = output_library['reference_name'].astype(str)
    output_library['sequence'] = output_library['seq_with_pam'].str[:-pam_length]

    species_to_file = dict(
        zip(species_df['species_name'], species_df[input_species_offtarget_column])
    )

    seqs: list[str] = list()
    library_path = os.path.join(output_directory, experiment_name + '_library.txt')
    with open(library_path, 'r') as f:
        for l in f.readlines():
            seqs.append(l.strip())

    OTF.write_guides_as_reads(f'{experiment_name}_output_lib', seqs)

    species_names = list(species_df['species_name'].unique())
    count = 0

    num_cpus = os.cpu_count()
    max_workers = max(1, (num_cpus - 1) if num_cpus is not None else 1)

    print(f'{bcolors.BLUE}>{bcolors.RESET} Looking for off-targets...')
    print(f'{bcolors.BLUE}>{bcolors.RESET} Done with {count}/{len(species_names)} species...', end='\r')

    with ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=_init_offtarget_worker,
        initargs=(
            output_directory,
            species_to_file,
            output_library,
            experiment_name,
            gene_or_genome_str,
            input_species_offtarget_dir,
            num_mismatches,
            seed_region_is_n_upstream_of_pam,
            protospacer_length,
            pam)
    ) as executor:

        future_to_species = {
            executor.submit(_find_targets_worker, target_species): target_species
            for target_species in species_names
        }

        for future in as_completed(future_to_species):
            target_species = future_to_species[future]

            try:
                result = future.result()
            except Exception as e:
                print()
                print(f'{bcolors.RED}> Error{bcolors.RESET}: failed for species "{target_species}": {e}')
                raise

            if not result.empty:
                created_dfs.append(result)

            count += 1
            print(f'{bcolors.BLUE}>{bcolors.RESET} Done with {count}/{len(species_names)} species...', end='\r')

    print()

    final_df = pandas.concat(created_dfs, ignore_index=True) if created_dfs else pandas.DataFrame()
    final_df.to_csv(f'{output_directory}/targets.csv', index=False)
    print(f'{bcolors.BLUE}>{bcolors.RESET} Done. Check {output_directory}/targets.csv for your off-targets report.')
