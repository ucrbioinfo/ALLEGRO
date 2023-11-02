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
num_trials: 10000
# ---

# ---
cluster_guides: True
seed_region_is_n_from_pam: 10
mismatches_allowed_after_seed_region: 2
# ---
'''

import os
import multiprocessing

def threaded_split(run_range):
    for run in run_range:
        new_exp_name = f'e1_run_{run}'
        
        context = {
            'experiment_name': new_exp_name,
            'input_species_path': 'data/input/fourdbs_hi_input_species.csv',
            'input_species_path_column': 'cds_file_name',
            'input_directory': 'data/input/cds/orthogroups/',
            'track': 'track_e',
            'scorer': 'dummy',
            'beta': 0,
            'mult': '1',
            'filter_repetitive': True,
            'mp_threshold': 0
        }

        config_name = f'temp_config_e1_{run}.yaml'
        with open(config_name, 'w') as f:
            f.write(config.format(**context))

        # Call ALLEGRO
        os.system('python src/main.py --config ' + config_name)


args_list = list()
runs = [[i for i in range(10 * n - 10, 10 * n)] for n in range(1, 11)]

with multiprocessing.Pool(processes=10) as pool:
    for r in runs:
        pool.apply_async(threaded_split, args=(r))

    # Close the pool and wait for all the processes to finish
    pool.close()
    pool.join()
