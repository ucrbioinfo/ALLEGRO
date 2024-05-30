import os

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

enable_solver_diagnostics: True
# ---
'''

# a6 ---------------------------------
a6_betas = [1300]
c = 100
while(True):
    if a6_betas[-1] >= 12330:
        break

    a6_betas.append(a6_betas[-1] + c)
    c += 100
a6_betas = [1269] + a6_betas[:-1] + [12330]
# ------------------------------------ 

tracks = [('track_a', 6)]
beta_lists = [a6_betas]

for idx, t in enumerate(tracks):
    track, mult = t

    for beta in beta_lists[idx]:
        new_exp_name = f'{track}{mult}_b{beta}'

        if os.path.exists(f'data/output/{new_exp_name}'):
            continue

        context = {
            'experiment_name': new_exp_name,
            'input_species_path': 'data/input/fourdbs_input_species.csv',
            'input_directory': 'data/input/cds/ortho_from_gff/',
            'input_species_path_column': 'ortho_file_name',
            'track': track,
            'beta': beta,
            'mult': mult,
            'patterns_to_exclude': ['TTTT'],
            'early_stopping_patience': 300
        }

        config_name = f'temp_config_{new_exp_name}.yaml'
        with open(config_name, 'w') as f:
            f.write(config.format(**context))

        # Call ALLEGRO
        os.system('python src/main.py --config ' + config_name)

        os.remove(config_name)