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
# ---

track: {track}
# ---

multiplicity: {mult}
# ---

beta: {beta}
# ---

#================================================================
# Advanced Settings
# ===============================================================

output_offtargets: False
report_up_to_n_mismatches: 3  # This may be 0, 1, 2, or 3
seed_region_is_n_upstream_of_pam: 12


input_species_offtarget_dir: 'data/input/cds/cds'
input_species_offtarget_column: 'cds_file_name'

max_threads: 128
# ---

scorer: {scorer}
# ---

filter_repetitive: {filter_repetitive}
# ---

mp_threshold: {mp_threshold}
# ---

num_trials: 10000
# ---

cluster_guides: False
mismatches_allowed_after_seed_region: 2
# ---
'''

import os
import multiprocessing

def threaded_split(run_range):
    for run in run_range:
        new_exp_name = f'a1_run_{run}'
        
        context = {
            'experiment_name': new_exp_name,
            'input_species_path': 'data/input/fourdbs_hi_gff_input_species.csv',
            'input_species_path_column': 'cds_file_name',
            'input_directory': 'data/input/cds/cds_from_gff/',
            'track': 'track_a',
            'scorer': 'dummy',
            'beta': 0,
            'mult': '1',
            'filter_repetitive': True,
            'mp_threshold': 0
        }

        config_name = f'temp_config_a1_{run}.yaml'
        with open(config_name, 'w') as f:
            f.write(config.format(**context))

        # Call ALLEGRO
        os.system('python src/main.py --config ' + config_name)

        os.remove(config_name)


args_list = list()
runs = [[i for i in range(1 * n - 1, 1 * n)] for n in range(1, 101)]

with multiprocessing.Pool(processes=70) as pool:
    for r in runs:
        pool.apply_async(threaded_split, args=(r,))

    # Close the pool and wait for all the processes to finish
    pool.close()
    pool.join()
