config_text = '''
# --------------------------------------------------------------------
# General settings
# --------------------------------------------------------------------
experiment_name: {experiment_name}
mode: {mode}  # 'from_genome', 'from_cds'
# objective: 'min'  # default: 'max' uses beta -- 'min' does not use beta
scorer: 'dummy'  # options: default: 'chopchop', 'dummy'
                 # dummy assigns a score of 1.0 to all guides
# solver_engine: 'GLOP'
# beta: 0  # Integer value: Used with objective: 'max'
         # The final size of the guides set must be <= than this. Think of it as your budget.

# (Affects performance)
# A higher number increases running time and reduces memory consumption.
# Integer value: Search all combinatorial possibilities unless there are more
# feasible solutions than this (then use randomized rounding)
# exhaustive_threshold: 5

# (Affects performance)
# False reduces running time and reduces memory consumption.
# If False, discards guides with 5 or more repeated 2-mers.
# For example, this cas9 guide will not be in the output: ACCACCACCACCACCACCAC
#  since it contains 7 AC 2-mers.
# Also discards guides containing repeating 4- or 5-mers such as AAAAA or TTTT.
include_repetitive: {include_repetitive}

# (Affects performance)
# A higher number increases running time while reducing memory consumption.
# Pre-select guides that hit only up to this number of species to act as representatives
# for these species. Set to 0 to disable (saves all guides to memory).
mp_threshold: {mp_threshold}

# (Affects performance)
# A higher number increases running time, but attempts to find a smaller cover set.
# Integer value: Only used if the number of feasible guides is above the exhaustive_threshold.
# How many times to run the randomized rounding algorithm? Any value under 1 disables trials.
num_trials: 10000

# graph: False
# output_csv: True
# --------------------------------------------------------------------


# --------------------------------------------------------------------
# Endonuclease settings
# --------------------------------------------------------------------
cas: 'cas9'  # default: 'cas9'
pam: 'NGG'  # default for cas9: 'NGG'
            # Supports 'NGG' which looks for 'A/C/T/G' + 'GG' 
            # Any other PAM will be treated as a literal string.
protospacer_length: 20  # Integer, default: 20
context_toward_five_prime: 0  # Integer, default: 5
context_toward_three_prime: 0  # Integer, default: 3
# --------------------------------------------------------------------


# --------------------------------------------------------------------
# Path Settings
# --------------------------------------------------------------------
output_directory: 'data/output/'
input_species_path: {input_species_path}
input_genomes_directory: 'data/input/genomes/'

# Required if mode is set to 'from_cds' -- ignored otherwise -- default: input_cds_directory: 'data/input/cds/orthogroups/''
input_cds_directory: 'data/input/cds/'  # Files inside must end with _cds.fna or _cds.faa or _cds.fasta

# Required if scorer is set to 'chopchop' -- ignored otherwise
# chopchop_scoring_method: 'DOENCH_2016'
# absolute_path_to_chopchop: '/home/amohs002/projects/research/chopchop/'  # directory only
# absolute_path_to_genomes_directory: '/home/amohs002/projects/research/ALLEGRO/data/input/genomes/'
'''

import os
import pandas as pd
import numpy as np
import subprocess

input_csv = 'data/input/final_standard_ncbi_input_species.csv'

df = pd.read_csv(input_csv)
dfs = np.array_split(df, len(df) // 50)  # split input into smaller dataframes of at most 50 rows each

processes: list = list()

for i, d in enumerate(dfs):
    new_exp_name = 'no_filter_ncbi_split_' + str(i)
    new_path = 'data/input/' + new_exp_name + '.csv'
    d.to_csv(new_path, index=False)  # type: ignore

    context = {
        'experiment_name': new_exp_name,
        'input_species_path': new_path,
        'mode': 'from_genome',
        'include_repetitive': True,
        'mp_threshold': 0
    }

    config_name = 'no_filter_temp_config_' + str(i) + '.yaml'
    with open(config_name, 'w') as f:
        f.write(config_text.format(**context))

    # cmd = 'python src/main.py --config ' + config_name

    # process = subprocess.Popen(cmd, shell=True)
    # processes.append(process)

    os.system('python src/main.py --config ' + config_name)

# for process in processes:
#     process.wait()


# for i, d in enumerate(dfs):
#     new_exp_name = 'filter_ncbi_split_' + str(i)
#     new_path = 'data/input/' + new_exp_name + '.csv'
#     d.to_csv(new_path, index=False)  # type: ignore

#     context = {
#         'experiment_name': new_exp_name,
#         'input_species_path': new_path,
#         'mode': 'from_genome',
#         'include_repetitive': False,
#         'mp_threshold': 4
#     }

#     config_name = 'filter_temp_config_' + str(i) + '.yaml'
#     with open(config_name, 'w') as f:
#         f.write(config_text.format(**context))

#     os.system('python src/main.py --config ' + config_name)

#     cmd = 'python src/main.py --config ' + config_name

#     process = subprocess.Popen(cmd, shell=True)
#     processes.append(process)

#     # os.system('python src/main.py --config ' + config_name)

# for process in processes:
#     process.wait()


# # ALL CDS
# context = {
#         'experiment_name': 'all_ncbi_cds',
#         'mode': 'from_cds',
#         'input_species_path': 'data/input/final_standard_ncbi_input_species.csv',
#         'include_repetitive': False,
#         'mp_threshold': 4
#         }

# config_name = 'all_cds_ncbi_temp_config.yaml'
# with open(config_name, 'w') as f:
#     f.write(config_text.format(**context))

# os.system('python src/main.py --config ' + config_name)