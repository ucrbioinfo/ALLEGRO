config_text = '''
# Some attributes are required conditionally. These are annotated via comments.

# Required settings
experiment_name: '339_ncbi'
mode: 'from_orthogroups'  # 'from_genome', 'from_orthogroups'
objective: 'min'  # default: 'max' uses beta -- 'min' does not use beta
scorer: 'chopchop'  # options: default: 'chopchop', 'random', 'dummy' -- dummy assigns a score of 1 to all guides
solver_engine: 'GLOP'
cas: 'cas9'  # default: 'cas9'
pam_regex: r'(?=(GG))'  # default for cas9: r'(?=(GG))'

beta: 10  # Used with objective: 'max' -- The final size of the guides set must be <= than this. Think of it as your budget.

input_species_path: 'data/input/final_standard_ncbi_input_species.csv'
input_genomes_directory: 'data/input/genomes/'
output_directory: 'data/339_output_ortho/'

# Search all combinatorial possibilities unless there are more
# feasible solutions than this (then use randomized rounding)
exhaustive_threshold: 0

# Only used if the number of feasible guides is above the exhaustive_threshold.
# How many times to run the randomized rounding algorithm?
num_trials: 10000

graph: False
output_csv: True
traceback: False

# Required if mode is set to 'from_orthogroups' -- ignored otherwise
input_cds_directory: 'data/input/cds/orthogroups/'  # Files inside must end with _cds.fna or _cds.faa or _cds.fasta
orthogroups_path: 'data/input/orthogroups/Orthogroups.tsv'

# Required if scorer is set to 'chopchop' -- ignored otherwise
chopchop_scoring_method: 'DOENCH_2016'
absolute_path_to_chopchop: '/home/amohs002/projects/research/chopchop/chopchop.py'
absolute_path_to_bowtie_build: '/home/amohs002/projects/research/chopchop/bowtie/bowtie-build'
absolute_path_to_genomes_directory: '/home/amohs002/projects/research/ALLEGRO/data/input/genomes/'
absolute_path_to_chopchop_config_local_json: '/home/amohs002/projects/research/chopchop/config_local.json'
'''

name = 'test'
n_species = 336

import os
import pandas as pd

for beta in range(0, n_species + 1):
    context = {
        'name': name,
        'beta': beta,
        'exhaustive_threshold': -1
    }
    with open('temp_config.yaml', 'w') as f:
        f.write(config_text.format(**context))

    # RUN main.py
    os.system('python src/main.py --config temp_config.yaml')
    # df = pd.read_csv('data/output/test_metrics.csv')
    # df.loc[df.index[-1], 'solved_with_exhaustive'] = 'No'
    # df.to_csv('data/output/test_metrics.csv', index=False)

    # context['exhaustive_threshold'] = 50
    # with open('temp_config.yaml', 'w') as f:
    #     f.write(config_text.format(**context))

    # os.system('python src/main.py --config temp_config.yaml')
    # df = pd.read_csv('data/output/test_metrics.csv')
    # df.loc[df.index[-1], 'solved_with_exhaustive'] = 'Yes'
    # df.to_csv('data/output/test_metrics.csv', index=False)