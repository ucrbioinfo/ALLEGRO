import os
import subprocess

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

scorer: 'ucrispr'
# ---

early_stopping_patience: {early_stopping_patience}
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

enable_solver_diagnostics: False
# ---
'''

# E1 ---------------------------------
e1_times = [10]
c = 60
while(True):
    if e1_times[-1] >= 10800:
        break

    e1_times.append(e1_times[-1] + c)
# ------------------------------------

def find_smallest_time(min_time, max_time):
    left = min_time
    right = max_time
    while left < right:
        mid = left + (right - left) // 2
        if response(mid):
            right = mid  # Move left if the response is True to find the smaller time
        else:
            left = mid + 1  # Move right if the response is False
    return left

# Define the maximum and minimum time
min_time = 10
max_time = 10800

tracks = [('track_e', 1)]
beta_lists = [e1_times]

for idx, t in enumerate(tracks):
    track, mult = t
    new_exp_name = f'{track}{mult}_time{time}'

    context = {
        'experiment_name': new_exp_name,
        'input_species_path': 'data/input/fourdbs_input_species.csv',
        'input_directory': 'data/input/cds/ortho_from_gff/',
        'input_species_path_column': 'ortho_file_name',
        'track': track,
        'beta': 0,
        'mult': mult,
        'patterns_to_exclude': ['TTTT'],
        'early_stopping_patience': time
    }

    config_name = f'temp_config_{new_exp_name}.yaml'
    with open(config_name, 'w') as f:
        f.write(config.format(**context))

    # Call ALLEGRO
    os.system('python src/main.py --config ' + config_name)

    command = ['python', 'src/main.py', '--config', config_name]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()



    os.remove(config_name)