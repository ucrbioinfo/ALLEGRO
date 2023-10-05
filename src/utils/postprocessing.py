# Functions imported by ALLEGRO. No need to run it manually.
import pandas
from collections import Counter

from utils.shell_colors import bcolors


def hamming_distance(str1: str, str2: str, length: int) -> int:
    """
    Calculates the Hamming distance between two strings up to a specified length.
    """
    if len(str1) < length or len(str2) < length:
        raise ValueError('Strings must have a length of at least {n}'.format(n=length))

    return sum(ch1 != ch2 for ch1, ch2 in zip(str1[:length], str2[:length]))


# def hamming_distance_full_length(s1: str, s2: str) -> int:
#     """Compute the Hamming distance between two strings."""
#     return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


# def average_distance_to_others(strings: list[str], idx: int) -> float:
#     """Compute the average Hamming distance of a string to all other strings."""
#     total_distance = sum(hamming_distance_full_length(strings[idx], s) for i, s in enumerate(strings) if i != idx)
#     return total_distance / (len(strings) - 1)


# def construct_representative(strings: list[str]) -> str:
#     """Construct a representative string by choosing the most common character at each position."""
#     representative = []
#     for i in range(len(strings[0])):
#         chars = [s[i] for s in strings]
#         most_common_char = Counter(chars).most_common(1)[0][0]
#         representative.append(most_common_char)
#     return ''.join(representative)


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