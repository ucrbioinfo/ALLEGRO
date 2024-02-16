import os
import re
import pandas
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed

from utils.shell_colors import bcolors
from utils.offtarget_finder import OfftargetFinder


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


def extract_complete_digits(s):
        return [int(num) for num in re.findall(r'\d+', s)]


def determine_ot_in_nonseed(digits, guide_len, seed_region_is_n_upstream_of_pam):
        return any(digit + 1 + seed_region_is_n_upstream_of_pam > guide_len for digit in digits)


# Function to use with multiprocessing
def find_targets(target_species: str,
                 experiment_name: str,
                 species_df: pandas.DataFrame,
                 output_library: pandas.DataFrame,
                 background_source: str,
                 input_species_offtarget_dir: str,
                 input_species_offtarget_column: str,
                 num_mismatches: int,
                 seed_region_is_n_upstream_of_pam: int) -> None:

    
    target_species_offtarget_dir = os.path.join(input_species_offtarget_dir, species_df[species_df['species_name'] == target_species][input_species_offtarget_column].values[0])
    
    # Create a Bowtie index for the target species if it doesn't exist (checks in the function)
    OTF = OfftargetFinder()
    OTF.run_bowtie_build(target_species, target_species_offtarget_dir, background_source)
    
    # Get all hits -- Needs to wait for index creation in the previous line
    all_targets = OTF.run_bowtie_against_other(f'{experiment_name}_output_lib', target_species, background_source, num_mismatches)
    
    all_targets['target_species'] = target_species
    all_targets['is_off_target'] = '1'  # initially, all hits are off-targets until proven otherwise
    all_targets['self_off_targets'] = '0'

    # These are on-targets -- everything else is off-target
    on_targets = pandas.merge(output_library[output_library['target'] == target_species],
                            all_targets,
                            left_on=['sequence', 'reference_name', 'strand', 'start_position'],
                            right_on=['sequence', 'reference_name', 'strand', 'start_position'])
    
    cols = ['is_off_target', 'self_off_targets', 'sequence', 'pam', 'mismatch', 'aligned_seq',
            'target_species', 'strand', 'reference_name', 'orthologous_to', 'start_position']
    
    on_targets = on_targets[cols]  # Only retain certain columns
    on_targets['is_off_target'] = '0'

    # Merging all_targets with on_targets on columns 'sequence',
    # 'reference_name', 'strand', 'start_position' to identify matching rows
    merged_df = all_targets.merge(on_targets,
                                on=['sequence', 'reference_name', 'strand', 'start_position'],
                                how='left', suffixes=('', '_new'))

    # Updating columns 'is_off_target' and 'orthologous_to' in all_targets where there's a match
    all_targets['orthologous_to'] = 'N/A'
    all_targets['is_off_target'] = merged_df['is_off_target_new'].combine_first(all_targets['is_off_target'])
    all_targets['orthologous_to'] = merged_df['orthologous_to'].combine_first(all_targets['orthologous_to'])
    all_targets = all_targets[cols]

    # make a list out of the mismatch locations
    all_targets['digits'] = all_targets['mismatch'].apply(extract_complete_digits)

    all_targets['ot_in_seed'] = all_targets.apply(lambda row: determine_ot_in_nonseed(row['digits'], len(row['sequence']), seed_region_is_n_upstream_of_pam), axis=1)
    
    # Remove rows where off-target occurs in the seed region. Remove useless columns.
    all_targets = all_targets[~all_targets['ot_in_seed'] == True]
    all_targets = all_targets.drop(columns=['digits', 'ot_in_seed'])

    OT = all_targets[all_targets['is_off_target'] == '1']
    for seq in OT['sequence'].unique():
        seq_off_target_species = set(OT[OT['sequence'] == seq]['target_species'])

        all_target_species = set(all_targets[all_targets['sequence'] == seq]['target_species'])

        if len(seq_off_target_species.intersection(all_target_species)) != 0:
            all_targets.loc[all_targets['sequence'] == seq, 'self_off_targets'] = '1'

    return all_targets


def report_offtargets(input_species_path: str,
                      output_directory: str,
                      experiment_name: str,
                      input_species_offtarget_dir: str,
                      input_species_offtarget_column: str,
                      num_mismatches: int,
                      seed_region_is_n_upstream_of_pam: int,
                      pam_length: int = 3):
    
    solution_path = os.path.join(output_directory, experiment_name + '.csv')

    background_source = 'genome' if 'genome' in input_species_offtarget_column else 'genes'

    OTF = OfftargetFinder()
    final_df = pandas.DataFrame()
    created_dfs: list[pandas.DataFrame] = list()  # For threads to deposit their results
    species_df = pandas.read_csv(input_species_path)
    output_library = pandas.read_csv(solution_path)

    # Pull back the start position of guides on negative strand
    # by PAM length -- to comply with bowtie hits and consider the PAM
    output_library.loc[output_library['strand'] == '-', 'start_position'] -= pam_length
    
    # Get all the library guides
    seqs = list()
    library_path = os.path.join(output_directory, experiment_name + '_library.txt')
    with open(library_path, 'r') as f:
        for l in f.readlines():
            seqs.append(l.strip())

    # Create reads out of guide sequences for Bowtie alignment
    OTF.write_guides_as_reads(f'{experiment_name}_output_lib', seqs)

    count = 0
    num_cpus = os.cpu_count()
    max_workers = num_cpus - 1 if num_cpus is not None else 1  # Ensure at least 1 worker is used
    
    print(f'{bcolors.BLUE}>{bcolors.RESET} Looking for off-targets...')
    print(f'Done with {count}/{len(output_library["target"].unique())} species...', end='\r')
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Prepare the partial function to include all necessary parameters except 'row'
        partial_find_targets = partial(find_targets,
                       experiment_name=experiment_name,
                       species_df=species_df,
                       output_library=output_library,
                       background_source=background_source,
                       input_species_offtarget_dir=input_species_offtarget_dir,
                       input_species_offtarget_column=input_species_offtarget_column,
                       num_mismatches=num_mismatches,
                       seed_region_is_n_upstream_of_pam=seed_region_is_n_upstream_of_pam
                       )
        
        # Submit tasks
        futures = [executor.submit(partial_find_targets, target_species) for target_species in output_library['target'].unique()]

        # Collect results as tasks complete
        for future in as_completed(futures):
            result = future.result()

            created_dfs.append(result)

            # time_elapsed += result["time_elapsed"]
            count += 1
            print(f'{bcolors.BLUE}>{bcolors.RESET} Done with {count}/{len(output_library["target"].unique())} species...', end='\r')
    print()

    final_df = pandas.concat(created_dfs, ignore_index=True)
    final_df.to_csv(f'{output_directory}/targets.csv', index=False)