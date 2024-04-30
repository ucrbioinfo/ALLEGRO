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

multiplicity: {multiplicity}
# ---

beta: 0
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

patterns_to_exclude: ['TTTT']

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

import os
import pandas as pd

increments = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2059]

tasks = [('track_a', 1), ('track_a', 7), ('track_e', 1)]

# Read in all species.
all_species = pd.read_csv('data/input/fourdbs_input_species.csv')
all_species = all_species.sample(frac=1, random_state=2647)  # Initial totally unnecessary shuffle.

def run_exp(track, mult, species_to_try):
    # Generate the input CSV file with this increment's species.
    new_csv_path = f'data/input/increment_exp_{track}{mult}_{len(species_to_try)}_species.csv'
    species_to_try.to_csv(new_csv_path, index=False)
    
    # Configure.
    new_exp_name = f'{track}{mult}_{len(species_to_try)}_species'
    context = {
        'experiment_name': new_exp_name,
        'input_species_path_column': 'ortho_file_name',
        'input_directory': 'data/input/cds/ortho_from_gff',
        'input_species_path': new_csv_path,
        'track': track,
        'multiplicity': mult,
    }

    # Write a new config file for ALLEGRO.
    config_name = f'{new_exp_name}_config.yaml'
    with open(config_name, 'w') as f:
        f.write(config.format(**context))

    # Run ALLEGRO.
    os.system('python src/main.py --config ' + config_name + ' -od data/output/species_increments --align_solution_to_input False')

    # Clean up.
    os.remove(config_name)


for inc in increments:
    for task in tasks:
        for _ in range(100):
            species_to_try = all_species.sample(inc, replace=False)  # randomly sample n=increment species.

            run_exp(task[0], task[1], species_to_try)