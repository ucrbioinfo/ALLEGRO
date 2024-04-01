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

# E1 ---------------------------------
e1_betas = [2700]
c = 100
while(True):
    if e1_betas[-1] >= 13001:
        break

    e1_betas.append(e1_betas[-1] + c)
    c += 100
e1_betas = [2633] + e1_betas
# ------------------------------------

# A1 ---------------------------------
a1_betas = [170]
c = 20
while(True):
    if a1_betas[-1] >= 2059:
        break

    a1_betas.append(a1_betas[-1] + c)
    c += 20
a1_betas = [166] + a1_betas
# ------------------------------------

# A7 ---------------------------------
a7_betas = [1400]
c = 100
while(True):
    if a7_betas[-1] >= 14413:
        break

    a7_betas.append(a7_betas[-1] + c)
    c += 100
a7_betas = [1355] + a7_betas
# ------------------------------------

            # 2633          # 166              # 1355
tracks = [('track_e', 1), ('track_a', 1), ('track_a', 7)]
beta_lists = [e1_betas, a1_betas, a7_betas]

for idx, t in enumerate(tracks):
    track, mult = t

    for beta_list in beta_lists[idx]:
        for beta in beta_list:
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