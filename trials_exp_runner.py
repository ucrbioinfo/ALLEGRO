config_text = '''
# ---
experiment_name: {experiment_name}
# ---

# ---
input_directory: 'data/input/cds/orthogroups/'
input_species_path: 'data/input/fourdbs_hi_input_species.csv'
input_species_path_column: 'cds_file_name'
# ---

# ---
track: {track}
# ---

# ---
multiplicity: {multiplicity}
# ---

# ---
beta: 0
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
num_trials: {num_trials}
# ---

# ---
cluster_guides: False
seed_region_is_n_from_pam: 10
mismatches_allowed_after_seed_region: 0
# ---
'''

import os
import multiprocessing

tracks = ['track_a', 'track_a', 'track_e']
multiplicities = [1, 7, 1]
trials = [1, 500, 1000, 1500, 2000, 5000, 10000]

def worker(args):
    track, multi, trial, run = args
    ae = track.split('_')[1]

    new_exp_name = f'trials_exp_t{ae}_m{multi}_ns_run{run}_{trial}'

    context = {
        'experiment_name': new_exp_name,
        'track': track,
        'multiplicity': multi,
        'num_trials': trial
    }

    config_name = f'{new_exp_name}_mm_config.yaml'
    with open(config_name, 'w') as f:
        f.write(config_text.format(**context))

    os.system('python src/main.py --config ' + config_name)

    os.remove(config_name)


for idx, track in enumerate(tracks):
    for trial in trials:
        for batch in range(5, 10):  # 5 batches of 20 runs each
            args_list = [(track, multiplicities[idx], trial, run) for run in range(batch * 10, (batch + 1) * 10)]
            
            with multiprocessing.Pool(processes=10) as pool:
                pool.map(worker, args_list)