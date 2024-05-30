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

# A1 ---------------------------------
a1_betas = [180]
c = 20
while(True):
    if a1_betas[-1] >= 2055:
        break

    a1_betas.append(a1_betas[-1] + c)
    c += 20
a1_betas = [177] + a1_betas[:-1] + [2055]
# ------------------------------------
   
tracks = [('track_a', 1)]
beta_lists = [a1_betas]

for idx, t in enumerate(tracks):
    track, mult = t

    for beta in beta_lists[idx]:
        new_exp_name = f'{track}{mult}_b{beta}'

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