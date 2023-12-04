# Functions imported by ALLEGRO. No need to run it manually.
import os
import pandas

from utils.offtarget_finder import OfftargetFinder
from utils.shell_colors import bcolors


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


def cluster_solution(solution_path: str, req_match_len: int, mm_allowed: int) -> int:
    print(f'{bcolors.BLUE}>{bcolors.RESET} Clustering guides: Guides in the same cluster have an identical sequence for the first {req_match_len} nucleotides after the PAM (3\' to 5\').')
    print(f'{bcolors.BLUE}>{bcolors.RESET} Guides in the same cluster may mismatch up to {mm_allowed} nucleotides after the seed region.')
    
    df = pandas.read_csv(solution_path)
    df['cluster'] = 0
    seqs = df.sequence.unique().tolist()

    clusters = cluster_strings(seqs, req_match_len, mm_allowed)

    # REPRESENTATIVE STRING CONSTRUCTION IS CURRENTLY DISABLED.
    #
    #
    # rep_strings = list()
    # scores = list()  # rep string has to be scored
    # new_clusters = list()
    # synthetics = list()
    # targets = list()

    for idx, cluster in enumerate(clusters):
        df.loc[df['sequence'].isin(cluster), 'cluster'] = idx

    #     synthetic = False
    #     representative_string = cluster[0]

    #     if len(cluster) > 1:
    #         # Step 1 & 2: Find the string with minimum average Hamming distance to all other strings
    #         avg_distances = [average_distance_to_others(cluster, i) for i in range(len(cluster))]
    #         min_avg_distance_string = cluster[avg_distances.index(min(avg_distances))]

    #         # Step 3: Construct a potential representative string
    #         constructed_representative = construct_representative(cluster)

    #         # Step 4: Compare and choose the representative string
    #         constructed_avg_distance = average_distance_to_others([constructed_representative] + cluster, 0)

    #         if constructed_avg_distance < min(avg_distances):
    #             representative_string = constructed_representative
    #         else:
    #             representative_string = min_avg_distance_string
            
    #         if representative_string not in df['sequence'].tolist():
    #             synthetic = True

    #     rep_strings.append(representative_string)
    #     targets.append(' | '.join([' - '.join(i) for i in zip(df[df['cluster'] == idx]['target'].tolist(), df[df['cluster'] == idx]['chromosome_or_gene'].tolist())]))
        
    #     if not synthetic:
    #         scores.append(df[df['sequence'] == representative_string]['score'].values[0])
    #     else:
    #         # TODO
    #         scores.append(1)

    #     synthetics.append(synthetic)
    #     new_clusters.append(idx)

    # rep_path = solution_path[:solution_path.find('.csv')] + '_rep_guides.csv'
    # pandas.DataFrame.from_dict({
    #     'sequence': rep_strings,
    #     'score': scores,
    #     'cluster': new_clusters,
    #     'synthetic': synthetics,
    #     'targets': targets
    # }).to_csv(rep_path, index=False)

    df.to_csv(solution_path, index=False)

    print(f'{bcolors.BLUE}>{bcolors.RESET} Done clustering. Added a new column \'cluster\' to {solution_path}')
    print(f'{bcolors.BLUE}>{bcolors.RESET} The output guide RNA set contains {len(clusters)} clusters.')

    return len(clusters)


def report_offtargets(input_species_path: str,
                      solution_path: str,
                      output_dir: str,
                      input_species_offtarget_dir: str,
                      input_species_offtarget_column: str,
                      experiment_name: str,
                      seed_region_is_n_from_pam: int,
                      num_mismatches: int,
                      pam_length: int = 3):
    
    background_source = 'genome' if 'genome' in input_species_offtarget_column else 'genes'
    
    OTF = OfftargetFinder()
    species_df = pandas.read_csv(input_species_path)
    output_library = pandas.read_csv(solution_path)

    # Pull back the start position of guides on negative strand
    # by PAM length -- to comply with bowtie hits and consider the PAM
    output_library.loc[output_library['strand'] == '-', 'start_position'] -= pam_length 
    seqs = output_library.sequence.unique().tolist()
    final_df = pandas.DataFrame()
    
    # Create reads out of guide sequences for Bowtie alignment
    OTF.write_guides_as_reads(f'{experiment_name}_output_lib', seqs)

    for target_species in output_library['target'].unique():
        target_species_offtarget_dir = os.path.join(input_species_offtarget_dir, species_df[species_df['species_name'] == target_species][input_species_offtarget_column].values[0])

        OTF.run_bowtie_build(target_species, target_species_offtarget_dir, background_source)

        # Get all hits
        all_targets = OTF.run_bowtie_against_other(f'{experiment_name}_output_lib', target_species, background_source, seqs, seed_region_is_n_from_pam, num_mismatches)
        
        all_targets['target_species'] = target_species
        all_targets['on_off_target'] = 'Off-Target'  # initially, all hits are offtargets until proven otherwise

        # These are on-targets -- everything else is off-target
        on_targets = pandas.merge(output_library[output_library['target'] == target_species],
                                 all_targets,
                                 left_on=['sequence', 'misc', 'strand', 'start_position'],
                                 right_on=['sequence', 'reference_name', 'strand', 'start_position'])
        
        cols = ['on_off_target', 'sequence', 'guide_w_pam', 'num_mutations', 'mismatch', 'aligned_seq',
                'target_species', 'strand', 'reference_name', 'orthologous_to', 'start_position']
        
        on_targets = on_targets[cols]  # Only retain certain columns
        on_targets['on_off_target'] = 'On-Target'

        # Merging all_targets with on_targets on columns 'sequence', 'reference_name', 'strand', 'start_position'
        #  to identify matching rows
        merged_df = all_targets.merge(on_targets,
                                    on=['sequence', 'reference_name', 'strand', 'start_position'],
                                    how='left', suffixes=('', '_new'))

        # Updating column 'on_off_target' in all_targets where there's a match
        all_targets['orthologous_to'] = 'N/A'
        all_targets['on_off_target'] = merged_df['on_off_target_new'].combine_first(all_targets['on_off_target'])
        all_targets['orthologous_to'] = merged_df['orthologous_to'].combine_first(all_targets['orthologous_to'])
        all_targets = all_targets[cols]
        
        final_df = pandas.concat([final_df, all_targets], ignore_index=True)

    final_df.to_csv(f'{output_dir}/targets.csv', index=False)
