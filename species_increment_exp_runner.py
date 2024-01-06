config_text = '''
# ---
experiment_name: '{experiment_name}'
# ---

# ---
input_directory: 'data/input/cds/orthogroups/'
input_species_path: '{input_species_path}'
input_species_path_column: 'cds_file_name'
# ---

# ---
track: '{track}'
# ---

# ---
multiplicity: {multiplicity}
# ---

# ---
beta: {beta}
# ---

# ---
scorer: 'dummy'
# ---

# ---
filter_repetitive: True
# ---

# ---
mp_threshold: 0
# ---

# ---
num_trials: 10000
# ---

# ---
cluster_guides: False
seed_region_is_n_from_pam: 10
mismatches_allowed_after_seed_region: 2
# ---
'''

import os
import pandas as pd
import multiprocessing

# Manually writing this out below b/c this generates out of order numbers with set cast.
# increments = set(np.geomspace(1, 2434, 30).astype(int))
increments = [1, 2, 3, 5, 6, 8, 11, 14, 19, 25, 32, 43, 56, 73, 96, 126, 165, 216, 283, 370, 484, 634, 830, 1086, 1421, 1860, 2434]

# Tasks for A7 and E1 to be assigned to processes.
tasks = list()
for i in increments: tasks.append(('track_a', 7, i))
for i in increments: tasks.append(('track_e', 1, i))

# Read in all species.
all_species = pd.read_csv('data/input/fourdbs_hi_input_species.csv')
all_species = all_species.sample(frac=1, random_state=2647)  # Initial totally unnecessary shuffle.

def run_exp(track, mult, increment):
    for run in range(100):  # 100 runs.
        species_to_try = all_species.sample(increment, replace=False)  # randomly sample n=increment species.

        # Generate the input CSV file with this increment's species.
        new_csv_path = f'data/input/increment_exp_{track}{mult}_{increment}_input_species.csv'
        species_to_try.to_csv(new_csv_path, index=False)
        
        # Configure.
        new_exp_name = f'{track}{mult}_ns_{increment}_species_run_{run}'
        context = {
            'experiment_name': new_exp_name,
            'track': track,
            'multiplicity': mult,
            'beta': 0,
            'input_species_path': new_csv_path,
        }

        # Write a new config file for ALLEGRO.
        config_name = f'{new_exp_name}_config.yaml'
        with open(config_name, 'w') as f:
            f.write(config_text.format(**context))

        # Run ALLEGRO.
        os.system('python src/main.py --config ' + config_name)

        # Clean up.
        os.remove(config_name)
        os.remove(new_csv_path)

# Reserve 40 processes and assign the next task as each process is done with its task.
with multiprocessing.Pool(processes=40) as pool:    
    for task in tasks:
        pool.apply_async(run_exp, args=(task[0], task[1], task[2]))

    pool.close()
    pool.join()