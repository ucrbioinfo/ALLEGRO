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

tasks = list()

for i in range(100, 2100, 100):
    tasks.append(('track_a', 7, i))

for i in range(100, 4600, 100):
    tasks.append(('track_e', 1, i))


def run_exp(track, multi, beta):
    ae = track.split('_')[1]

    for run in range(100):
        new_exp_name = f'{ae}{multi}_b{beta}_run_{run}'

        context = {
            'experiment_name': new_exp_name,
            'track': track,
            'multiplicity': multi,
            'beta': beta
        }

        config_name = f'{new_exp_name}_config.yaml'
        with open(config_name, 'w') as f:
            f.write(config.format(**context))

        os.system('python src/main.py --config ' + config_name)
        os.remove(config_name)

        if len(os.listdir(f'data/output/{new_exp_name}')) == 1:
            break


with multiprocessing.Pool(processes=40) as pool:
    for task in tasks:
        pool.apply_async(run_exp, args=(task[0], task[1], task[2]))

    # Close the pool and wait for all the processes to finish
    pool.close()
    pool.join()
