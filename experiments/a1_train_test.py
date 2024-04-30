import os
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
import multiprocessing
from sklearn.model_selection import KFold
from src.utils.test_guide_finder import GuideFinder


def find_all_matches(list1, list2):
    matches = list()

    for index, item in enumerate(list2):
        if item in list1:
            matches.append(index)
        
    return matches


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
    gene_name = gene_match.group(1) if gene_match is not None else record.id

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
#================================================================
# General Settings
#================================================================
experiment_name: {experiment_name}
# ---

# ---------------------------------------------------------------
# Path Settings
# ---------------------------------------------------------------

input_directory: {input_directory}
input_species_path: {input_species_path}
input_species_path_column: {input_species_path_column}
# ---------------------------------------------------------------

track: {track}
# ---

multiplicity: {mult}
# ---

beta: {beta}
# ---

scorer: 'dummy'
# ---

early_stopping_patience: 60
# ---

filter_by_gc: True
gc_max: 0.7
gc_min: 0.3
# ---

#================================================================
# Advanced Settings
# ===============================================================

patterns_to_exclude: {patterns_to_exclude}

output_offtargets: False
report_up_to_n_mismatches: 3  # This may be 0, 1, 2, or 3
seed_region_is_n_upstream_of_pam: 12

input_species_offtarget_dir: 'data/input/cds/cds_from_gff'
input_species_offtarget_column: 'cds_file_name'
# ---

cluster_guides: False
mismatches_allowed_after_seed_region: 2  # Integer value, default: 2

enable_solver_diagnostics: True
# ---
'''

def threaded_split(args):
    train_index, test_index, split, run_range = args
    train, test = df.iloc[train_index], df.iloc[test_index]

    train_new_path = f'data/input/train_a1_cv_split_{split}.csv'
    train.to_csv(train_new_path, index=False)

    test_new_path = f'data/input/test_a1_cv_split_{split}.csv'
    test.to_csv(test_new_path, index=False)

    for run in run_range:
        new_exp_name = f'a1_cv_split_{split}_run_{run}'

        if os.path.exists(f'data/output/test_{new_exp_name}_results.csv'):
            return
        
        context = {
            'experiment_name': new_exp_name,
            'input_species_path': train_new_path,
            'input_directory': 'data/input/cds/ortho_from_gff/',
            'input_species_path_column': 'ortho_file_name',
            'track': 'track_a',
            'beta': 0,
            'mult': '1',
            'patterns_to_exclude': ['TTTT'],
        }

        config_name = f'temp_config_a1_cv_{split}_{run}.yaml'
        with open(config_name, 'w') as f:
            f.write(config.format(**context))

        # Call ALLEGRO
        os.system('python src/main.py --config ' + config_name)
        
        library = list()

        with open(f'data/output/{new_exp_name}/{new_exp_name}_library.txt', 'r') as f:
            for line in f.readlines():
                library.append(line.strip())

        test_df = pd.read_csv(f'data/input/test_a1_cv_split_{split}.csv')[['species_name', context['input_species_path_column']]]

        df_species_name_list = list()
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
            species_name = row['species_name']

            species_is_covered_n_times = 0
            
            records_path = base_path + row[context['input_species_path_column']]
            records = list(SeqIO.parse(open(records_path), 'fasta'))

            for record in records:
                gene_name, protein_id, ortho_prot_id, ortho_gene_name = get_record_metadata(record)

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
                matches = find_all_matches(library, guides_list)
                for match in matches:
                    df_species_name_list.append(species_name)
                    df_covered_list.append(guides_list[match])
                    df_locations.append(locations_list[match])
                    df_strands.append(strands_list[match])
                    df_gene_names.append(gene_name)
                    df_protein_ids.append(protein_id)
                    df_ortho_prot_ids.append(ortho_prot_id)
                    df_ortho_gene_names.append(ortho_gene_name)
                    df_lib_partial_match.append('N/A')
                    df_mismatch.append('N/A')
                    species_is_covered_n_times += 1
            
            # Attempt to find partial matches.
            mm_allowed = 2
            req_match_len = 10

            for record in records:
                gene_name, protein_id, ortho_prot_id, ortho_gene_name = get_record_metadata(record)

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
                        
                        if (lib_guide[req_match_len:] == test_guide[req_match_len:]) and \
                            (distance_after_seed <= mm_allowed) and \
                                (distance_after_seed > 0):

                            df_species_name_list.append(species_name)
                            df_covered_list.append(test_guide)
                            df_locations.append(locations_list[idx])
                            df_strands.append(strands_list[idx])
                            df_gene_names.append(gene_name)
                            df_protein_ids.append(protein_id)
                            df_ortho_prot_ids.append(ortho_prot_id)
                            df_ortho_gene_names.append(ortho_gene_name)
                            df_mismatch.append(distance_after_seed)
                            df_lib_partial_match.append(lib_guide)
                            species_is_covered_n_times += 1

            # No match.
            if species_is_covered_n_times == 0:
                df_species_name_list.append(species_name)
                df_covered_list.append('Not covered')
                df_strands.append('Not covered')
                df_locations.append('Not covered')
                df_gene_names.append('Not covered')
                df_protein_ids.append('Not covered')
                df_ortho_prot_ids.append('Not covered')
                df_ortho_gene_names.append('Not covered')
                df_mismatch.append('Not covered')
                df_lib_partial_match.append('Not covered')

        test_df = pd.DataFrame.from_dict({
            'species': df_species_name_list,
            'covered_by': df_covered_list,
            'partial_match_in_library': df_lib_partial_match,
            'mismatch': df_mismatch,
            'strand': df_strands,
            'location': df_locations,
            'gene_name': df_gene_names,
            'protein_id': df_protein_ids,
            'ortho_prot_id': df_ortho_prot_ids,
            'ortho_gene_name': df_ortho_gene_names,
        })

        test_df.to_csv(f'data/output/test_{new_exp_name}_results.csv', index=False)

        os.remove(config_name)


df = pd.read_csv('data/input/fourdbs_input_species.csv')

# Set the random seed for reproducibility.
seed = 42
np.random.seed(seed)

# Shuffle the data.
df = df.sample(frac=1, random_state=seed).reset_index(drop=True)

# Initialize the KFold cross-validator.
kfold = KFold(n_splits=10)

split = 1
args_list = list()

# runs = [[i for i in range(1 * n - 1, 1 * n)] for n in range(1, 101)]

runs = [[0]]

with multiprocessing.Pool(processes=1) as pool:
    for train_index, test_index in kfold.split(df):

        for r in runs:
            payload = (train_index, test_index, split, r)
            pool.apply_async(threaded_split, args=(payload,))

        split += 1
    
    # Close the pool and wait for all the processes to finish
    pool.close()
    pool.join()
