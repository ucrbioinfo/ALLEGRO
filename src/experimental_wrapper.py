config_text = '''experiment_name: {name}
mode: 'from_orthogroups'
scorer: 'chopchop'
solver_engine: 'GLOP'
cas: 'cas9'
pam_regex: r'(?=(GG))'

beta: {beta}
exhaustive_threshold: {exhaustive_threshold}

num_trials: 10000

traceback: False
graph: False
output_csv: True

input_species_path: 'data/input/input_species.csv'
input_genomes_directory: 'data/input/genomes/'
output_directory: 'data/output/'

orthogroups_path: 'data/input/orthogroups/Orthogroups.tsv'
input_cds_directory: 'data/input/cds/orthogroups/'  # Files inside must end with _cds.fna or _cds.faa or _cds.fasta

absolute_path_to_chopchop: '/home/amohs002/projects/research/chopchop/chopchop.py'
absolute_path_to_bowtie_build: '/home/amohs002/projects/research/chopchop/bowtie/bowtie-build'
absolute_path_to_genomes_directory: '/home/amohs002/projects/research/application_v3/data/input/genomes/'
chopchop_scoring_method: 'DOENCH_2016'
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
    df = pd.read_csv('data/output/test_metrics.csv')
    df.loc[df.index[-1], 'solved_with_exhaustive'] = 'No'
    df.to_csv('data/output/test_metrics.csv', index=False)

    # context['exhaustive_threshold'] = 50
    # with open('temp_config.yaml', 'w') as f:
    #     f.write(config_text.format(**context))

    # os.system('python src/main.py --config temp_config.yaml')
    # df = pd.read_csv('data/output/test_metrics.csv')
    # df.loc[df.index[-1], 'solved_with_exhaustive'] = 'Yes'
    # df.to_csv('data/output/test_metrics.csv', index=False)