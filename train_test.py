import os
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from sklearn.model_selection import KFold
from src.utils.guide_finder import GuideFinder


def find_first_match(list1, list2):
    for index, item in enumerate(list2):
        if item in list1:
            return index
    return -1  # If no match found.


def hamming_distance(str1: str, str2: str, length: int) -> int:
    """
    Calculates the Hamming distance between two strings up to a specified length.
    """
    if len(str1) < length or len(str2) < length:
        raise ValueError('Strings must have a length of at least {n}'.format(n=length))

    return sum(ch1 != ch2 for ch1, ch2 in zip(str1[:length], str2[:length]))


def get_record_metadata(record):
    # This gene's own name -- Usually N/A for unannotated CDS files or chromosomes/scaffolds
    gene_match = re.search(gene_regex, record.description)
    gene_name = gene_match.group(1) if gene_match is not None else 'N/A'

    # This gene's own protein id e.g., XP_022674739.1
    protein_id_match = re.search(protein_id_regex, record.description)
    protein_id = protein_id_match.group(1) if protein_id_match is not None else 'N/A'

    # For example, [orthologous_to_ref_protein=XP_022674739.1], extracts XP_022674739.1
    ortho_prot_to_match = re.search(orthologous_protein_regex, record.description)
    ortho_prot_id = ortho_prot_to_match.group(1) if ortho_prot_to_match is not None else 'N/A'

    # For example, [orthologous_to_gene=HIS7], extracts HIS7
    ortho_gene_to_match = re.search(orthologous_name_regex, record.description)
    ortho_gene_name = ortho_gene_to_match.group(1) if ortho_gene_to_match is not None else 'N/A'

    return gene_name, protein_id, ortho_prot_id, ortho_gene_name


gene_regex = r'\[gene=(.*?)\]'
tag_regex = r'\[locus_tag=(.*?)\]'
protein_id_regex = r'\[protein_id=(.*?)\]'
reference_species_regex = r'\[ref_species=(.*?)\]'
orthologous_name_regex = r'\[orthologous_to_gene=(.*?)\]'
orthologous_protein_regex = r'\[orthologous_to_ref_protein=(.*?)\]'

config = '''
# ---
experiment_name: {experiment_name}
# ---

# ---
input_directory: {input_directory}
input_species_path: {input_species_path}
input_species_path_column: {input_species_path_column}
# ---

# ---
track: {track}
# ---

# ---
multiplicity: {mult}
# ---

# ---
beta: {beta}
# ---

# ---
scorer: {scorer}
# ---

# ---
filter_repetitive: {filter_repetitive}
# ---

# ---
mp_threshold: {mp_threshold}
# ---

# ---
num_trials: 1000
# ---

# ---
cluster_guides: True
seed_region_is_n_from_pam: 10
mismatches_allowed_after_seed_region: 4
# ---
'''

df = pd.read_csv('data/input/final_ncbi_concat_1k_input_species.csv')

# Set the random seed for reproducibility.
seed = 42
np.random.seed(seed)

# Shuffle the data.
df = df.sample(frac=1, random_state=seed).reset_index(drop=True)

# Initialize the KFold cross-validator.
kfold = KFold(n_splits=10)

split = 1
for train_index, test_index in kfold.split(df):
    train, test = df.iloc[train_index], df.iloc[test_index]

    new_exp_name = f'all_cds_split_{split}'
    train_new_path = f'data/input/train_{new_exp_name}.csv'
    train.to_csv(train_new_path, index=False)

    test_new_path = f'data/input/test_{new_exp_name}.csv'
    test.to_csv(test_new_path, index=False)

    context = {
        'experiment_name': new_exp_name,
        'input_species_path': train_new_path,
        'input_species_path_column': 'cds_file_name',
        # 'input_species_path_column': 'genome_file_name',
        'input_directory': 'data/input/cds/orthogroups/',
        # 'input_directory': 'data/input/genomes/',
        'track': 'track_a',
        'scorer': 'dummy',
        'beta': 0,
        'mult': '1',
        'filter_repetitive': True,
        'mp_threshold': 0
    }

    config_name = f'temp_config_{split}.yaml'
    with open(config_name, 'w') as f:
        f.write(config.format(**context))

    # Call ALLEGRO
    os.system('python src/main.py --config ' + config_name)

    split_out = pd.read_csv(f'data/output/{new_exp_name}/{new_exp_name}.csv')

    library = split_out.sequence.unique()

    test_df = pd.read_csv(f'data/input/test_{new_exp_name}.csv')
    test_df = test_df.drop(columns=['cds_file_name'])

    df_covered_list = list()
    df_strands = list()
    df_mismatch = list()
    df_locations = list()
    df_gene_names = list()
    df_protein_ids = list()
    df_ortho_prot_ids = list()
    df_ortho_gene_names = list()
    df_lib_partial_match = list()

    gf = GuideFinder()
    base_path = context['input_directory']

    for label, row in test_df.iterrows():
        species_is_covered = False
        records_path = base_path + row.file_name
        records = list(SeqIO.parse(open(records_path), 'fasta'))

        for record in records:
            (guides_list,
            _guides_context_list,
            strands_list,
            locations_list) = gf.identify_guides_and_indicate_strand(
                pam='NGG',
                sequence=str(record.seq).upper(),
                protospacer_length=20,
                context_toward_five_prime=0,
                context_toward_three_prime=0,
                filter_repetitive=False
            )

            # Attempt to find an exact match.
            match = find_first_match(library, guides_list)
            if match != -1:
                gene_name, protein_id, ortho_prot_id, ortho_gene_name = get_record_metadata(record)

                df_covered_list.append(guides_list[match])
                df_locations.append(locations_list[match])
                df_strands.append(strands_list[match])
                df_gene_names.append(gene_name)
                df_protein_ids.append(protein_id)
                df_ortho_prot_ids.append(ortho_prot_id)
                df_ortho_gene_names.append(ortho_gene_name)
                df_lib_partial_match.append('N/A')
                df_mismatch.append('N/A')

                species_is_covered = True
                break  # This species has been covered.
        
        if species_is_covered == False:
            partial_matches = list()

            # Attempt to find partial matches.
            mm_allowed = 2
            req_match_len = 10
            
            best_record = None
            smallest_mm = mm_allowed + 1
            best_lib_guide = ''
            best_test_guide = ''
            best_test_guide_loc = 0
            best_test_guide_strand = ''

            for record in records:
                (guides_list,
                _guides_context_list,
                strands_list,
                locations_list) = gf.identify_guides_and_indicate_strand(
                    pam='NGG',
                    sequence=str(record.seq).upper(),
                    protospacer_length=20,
                    context_toward_five_prime=0,
                    context_toward_three_prime=0,
                    filter_repetitive=False
                )

                for idx, test_guide in enumerate(guides_list):
                    for lib_guide in library:
                        
                        distance_after_seed = hamming_distance(lib_guide, test_guide, req_match_len)
                        
                        if lib_guide[req_match_len:] == test_guide[req_match_len:] and distance_after_seed < mm_allowed and distance_after_seed < smallest_mm:
                            
                            best_record = record
                            best_lib_guide = lib_guide
                            best_test_guide = test_guide
                            smallest_mm = distance_after_seed
                            best_test_guide_loc = locations_list[idx]
                            best_test_guide_strand = strands_list[idx]

                            species_is_covered = True

            if species_is_covered:
                gene_name, protein_id, ortho_prot_id, ortho_gene_name = get_record_metadata(best_record)

                df_covered_list.append(best_test_guide)
                df_locations.append(best_test_guide_loc)
                df_strands.append(best_test_guide_strand)
                df_gene_names.append(gene_name)
                df_protein_ids.append(protein_id)
                df_ortho_prot_ids.append(ortho_prot_id)
                df_ortho_gene_names.append(ortho_gene_name)
                df_mismatch.append(smallest_mm)
                df_lib_partial_match.append(best_lib_guide)
                

            # No match.
            elif species_is_covered == False:
                df_covered_list.append('Not covered')
                df_strands.append('Not covered')
                df_locations.append('Not covered')
                df_gene_names.append('Not covered')
                df_protein_ids.append('Not covered')
                df_ortho_prot_ids.append('Not covered')
                df_ortho_gene_names.append('Not covered')
                df_mismatch.append('Not covered')
                df_lib_partial_match.append('Not covered')
            
    test_df['covered_by'] = df_covered_list
    test_df['partial_match_in_library'] = df_lib_partial_match
    test_df['mismatch'] = df_mismatch
    test_df['strand'] = df_strands
    test_df['location'] = df_locations
    test_df['gene_name'] = df_gene_names
    test_df['protein_id'] = df_protein_ids
    test_df['ortho_prot_id'] = df_ortho_prot_ids
    test_df['ortho_gene_name'] = df_ortho_gene_names

    test_df.to_csv(f'test_{new_exp_name}_results.csv', index=False)

    split += 1