# Functions imported by ALLEGRO. No need to run it manually.
import os
import sys
import yaml
import argparse

from utils.shell_colors import bcolors


def greet() -> None:
    print(f'{bcolors.BLUE}>{bcolors.RESET} Welcome to {bcolors.ORANGE}ALLEGRO{bcolors.RESET}.')
    print(f'{bcolors.BLUE}>{bcolors.RESET} All unspecified command-line arguments default to the values in config.yaml.')


# TO use CHOPCHOP, ALLEGRO needs a conda environment called 'chopchop' with all the appropriate
# CHOPCHOP dependencies and python 2.7. If not found, an error is printed out.
def conda_env_exists(env_name: str) -> bool:
    if os.system(f'conda env list | grep {env_name} > /dev/null') == 0:
        return True
    else:
        return False
    

def parse_configurations() -> argparse.Namespace:
    config_parser = argparse.ArgumentParser(add_help=False)

    help = '''
    - The config file to use. Must be placed in the root folder.
    '''
    config_parser.add_argument(
        '-c',
        '--config',
        type=argparse.FileType(mode='r'),
        default='config.yaml',
        help=help,
    )
    
    config_args, remaining_args = config_parser.parse_known_args()

    config_arg_dict = vars(config_args)
    if config_args.config:
        config_arg_dict.update(yaml.load(config_args.config, Loader=yaml.FullLoader))

    parser = argparse.ArgumentParser(
        prog='ALLEGRO: An algorithm for a linear program to enhance guide RNA optimization',
        description='Find the smallest set of guide RNAs that cut through all species.',
        epilog="For more info, visit https://github.com/AmirUCR/allegro",
        parents=[config_parser],
    )

    help = '''
    - Name of the experiment. Output directory and file will be labeled with this.
    '''
    parser.add_argument(
        '-n',
        '--experiment_name',
        type=str,
        help=help
    )

    help = '''
    - Which scoring method to use? Default: 'dummy'
    - Options are 'chopchop', 'ucrispr', 'dummy' where 'dummy' assigns a score of 1.0 to all guides.
    '''
    parser.add_argument(
        '--scorer',
        type=str,
        default='dummy',
        help=help,
    )

    help = '''
    - Which track to use?
    - Options are 'track_e', 'track_a'
    - track_a: any of the genes in container_names can be targeted. There will be multiplicity targets per guide container (species or genes).
    - track_e: each species/gene has to be targeted at least multiplicity times,
    '''
    parser.add_argument(
        '--track',
        type=str,
        help=help,
    )

    help = '''
    - (Affects performance) Default: True
    - True reduces running time and reduces memory consumption.
    - If True, discards guides with 5 or more repeated 2-mers.
    - For example, this cas9 guide will not be in the output: ACCACCACCACCACCACCAC since it contains 7 'AC' 2-mers.
    - Also discards guides containing repeating 4- or 5-mers such as AAAAA or TTTT.
    '''
    parser.add_argument(
        '--filter_repetitive',
        type=bool,
        default=True,
        help=help
    )

    # help = '''
    # - Which cas endonuclease to use? Default: 'cas9'. Options are: 'cas9'.
    # '''
    # parser.add_argument(
    #     '--cas',
    #     type=str,
    #     default='cas9',
    #     help=help,
    # )

    help = '''
    - Ensures to cut each species/gene at least this many times. Default: 1
    '''
    parser.add_argument(
        '--multiplicity',
        type=int,
        default=1,
        help=help,
    )

    help = '''
    - (Affects performance) Post-processing. Default: True
    - Compresses the output guide set by clustering similar guides; adds a new column to output.csv
    '''
    parser.add_argument(
        '--cluster_guides',
        type=bool,
        default=True,
        help=help,
    )
    
    help = '''
    - Only used when cluster_guides is True. Default: 10
    - Choose how many nucleotides toward the 5' is considered seeds region.
    - For example, when set to 10 for 5'- AAACTGGTACTGACTGACCGNGG -3', TGACTGACCG is considered the seed region.
    '''
    parser.add_argument(
        '--seed_region_is_n_from_pam',
        type=int,
        default=10,
        help=help,
    )

    help = '''
    - Only used when cluster_guides is True. Default: 5
    - Choose how many mismatches after the seed region is allowed for placing guides with an identical seed region in the same cluster.
    - For example, guides AAACTGGTACTGACTGACCG and AAACCCCTACTGACTGACCG are placed in the same cluster when this number is 3 or higher.
    '''
    parser.add_argument(
        '--mismatches_allowed_after_seed_region',
        type=int,
        default=4,
        help=help,
    )

    help = '''
    - Beta represents a loose budge or threshold for the maximum number of guides you would like in the final covering set.
    - This is *not* a guaranteed maximum due to the hardness of the problem. See the topic on Integrality Gap.
    '''
    parser.add_argument(
        '-b',
        '--beta',
        type=int,
        help=help,
    )

    help = '''
    - (Affects performance) Default: True
    - True decreases subsequent running time but increases memory consumption.
    - If True, saves guides and their scores to data/cache/saved_guides.pickle for future lookups.
    - If False, each guide has to be scored again every time ALLEGRO is run.
    '''
    parser.add_argument(
        '--use_secondary_memory',
        type=bool,
        default=True,
        help=help
    )

    help = '''
    - The .csv file that includes three columns: species_name, genome_file_name, cds_file_name. genome_file_name rows start with species_name + _genomic.fna. cds_file_name rows start with species_name + _cds.fna.
    - Genome files must be placed in data/input/genomes while cds files must be placed in data/input/cds.
    '''
    parser.add_argument(
        '--input_species_path',
        type=str,
        help=help,
    )

    help = '''
    - The desired column name of the .csv file that includes three columns: species_name, genome_file_name, cds_file_name. genome_file_name rows start with species_name + _genomic.fna. cds_file_name rows start with species_name + _cds.fna.
    - Genome files must be placed in data/input/genomes while cds files must be placed in data/input/cds.
    '''
    parser.add_argument(
        '--input_species_path_column',
        type=str,
        help=help,
    )

    parser.add_argument(
        '--output_directory',
        type=str,
        default="data/output/"
    )

    help = '''
    - Files in this directory must end with .fna.
    '''
    parser.add_argument(
        '--input_directory',
        type=str,
        help=help,
    )

    help = '''
    - Only used in chopchop is selected as the scorer. Which chopchop scoring model to use?
    '''
    parser.add_argument(
        '--chopchop_scoring_method',
        type=str,
        help=help,
    )

    help = '''
    - Point to the directory where chopchop is located.
    - Download chopchop from here https://bitbucket.org/valenlab/chopchop/src/master/"
    '''
    parser.add_argument(
        '--absolute_path_to_chopchop',
        type=str,
        help=help,
    )

    help = '''
    - Points to the directory where all species' genomes are placed.
    '''
    parser.add_argument(
        '--absolute_path_to_genomes_directory',
        type=str,
        help=help,
    )

    # help = '''
    # - Search all combinatorial possibilities unless there are more feasible solutions than this (then uses randomized rounding).
    # - The number of calculations grows exponentially.
    # - ONLY increase this parameter when the number of guides in under about twenty.
    # '''
    # parser.add_argument(
    #     '--exhaustive_threshold',
    #     type=int,
    #     help=help,
    # )

    help = '''
    - Only used if the number of feasible guides is above the exhaustive_threshold.
    - How many times to run the randomized rounding algorithm?
    '''
    parser.add_argument(
        '--num_trials',
        type=int,
        help=help,
    )

    help = '''
    - A higher number increases running time while decreasing memory consumption.
    - Pre-select guides that hit only up to this number of species/genes to act as representatives for them.
    '''
    parser.add_argument(
        '--mp_threshold',
        type=int,
        help=help,
    )
    
    parser.set_defaults(**config_arg_dict)
    args = parser.parse_args(remaining_args)

    return args


def check_and_fix_configurations(args: argparse.Namespace) -> tuple[argparse.Namespace, dict]:
    # If CHOPCHOP is the selected scorer, set the chopchop scoring method.
    if 'chopchop' in args.scorer or 'CHOPCHOP' in args.scorer:
        # Set paths for CHOPCHOP
        args.absolute_path_to_chopchop = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scorers/chopchop/'))
        args.absolute_path_to_genomes_directory = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../..', 'data/input/genomes/'))
        
        if conda_env_exists('chopchop') == False:
            print(f'{bcolors.RED}> Warning{bcolors.RESET}: You have selected CHOPCHOP as the guide RNA scorer. {bcolors.ORANGE}ALLEGRO{bcolors.RESET} will attempt to run CHOPCHOP with a conda environment called "chopchop" that must include python 2.7 and all the other python libraries for running CHOPCHOP.\n')
            print(f'{bcolors.BLUE}>{bcolors.RESET} For more info, see here https://bitbucket.org/valenlab/chopchop/src/master/\n')
            print(f'{bcolors.BLUE}>{bcolors.RESET} {bcolors.ORANGE}ALLEGRO{bcolors.RESET} ships with CHOPCHOP and Bowtie so you do not need to download the repository or set any paths manually. You only need to create a conda environment called "chopchop" with python 2.7, and install any required scorer libraries in it such as scikit-learn, keras, theano, and etc.\n')
            print(f'{bcolors.BLUE}>{bcolors.RESET} When CHOPCHOP is selected as the scorer, you need to place the genome fasta files of every input species in data/input/genomes/ to be used with Bowtie.')
            sys.exit(1)

        split = args.scorer.split('_')
        join = '_'.join(split[1:])

        if join == '':
            print(f'{bcolors.RED}> Warning{bcolors.RESET}: Please select a scorer for CHOPCHOP in config.yaml. Exiting.')
            sys.exit(1)

        args.chopchop_scoring_method = join.upper()
        args.scorer = 'chopchop'

    if args.track not in ['track_a', 'track_e']:
        print(f'Unknown track {args.track} selected in config.yaml. Exiting.')
        sys.exit(1)

    if args.multiplicity < 1:
        print(f'{bcolors.RED}> Warning{bcolors.RESET}: Multiplicity is set to {args.multiplicity}, a value smaller than 1. Auto adjusting multiplicity to 1.')
        args.multiplicity = 1

    if args.mp_threshold <= 0:
        print(f'{bcolors.BLUE}>{bcolors.RESET} mp_threshold is set to {args.mp_threshold} and thus disabled. Saving all guides to memory.')
        args.mp_threshold = 0

    if args.mp_threshold > 0 and args.mp_threshold < args.multiplicity:
        print(f'{bcolors.RED}> Warning{bcolors.RESET}: mp_threshold is set to {args.mp_threshold}, a positive value smaller than the multiplicity {args.multiplicity}.')
        print(f'{bcolors.ORANGE}ALLEGRO{bcolors.RESET} cannot remove all but {args.mp_threshold} guides from each container and still ensure each container is targeted at least {args.multiplicity} times.')
        print(f'{bcolors.BLUE}>{bcolors.RESET} Auto adjusting mp_threshold to be equal to multiplicity. You may also set mp_threshold to 0 to disable this feature.')
        args.mp_threshold = args.multiplicity

    if args.beta <= 0:
        args.beta = 0
        print(f'{bcolors.BLUE}>{bcolors.RESET} Beta is set to {args.beta} and thus disabled.')
        print(f'{bcolors.BLUE}>{bcolors.RESET} {bcolors.ORANGE}ALLEGRO{bcolors.RESET} will minimize the set size.')

        if args.scorer != 'dummy':
            print(f'{bcolors.BLUE}>{bcolors.RESET} Scorer is set to {args.scorer}. {bcolors.ORANGE}ALLEGRO{bcolors.RESET} will score the guides for information only and will not use them in calculations.')
        
        
    if args.scorer == 'dummy' and args.beta > 0:
        # No feasible solutions if there are fewer guides than beta
        # Say there are 5 species, 5 guides total, and beta is set to 1. Say that none of the species share any guides.
        # This will ask ALLEGRO to find 1 guide out of 5 to cover all 5 species. There is no solution.
        print(f'{bcolors.BLUE}>{bcolors.RESET} The scorer is set to dummy and beta to the positive value {args.beta}. {bcolors.ORANGE}ALLEGRO{bcolors.RESET} will try to find (approximately) {args.beta} guides to cover all guide containers.')
        print(f'{bcolors.RED}> Warning{bcolors.RESET}: {bcolors.ORANGE}ALLEGRO{bcolors.RESET} may find that there is no feasible solution if the number of shared guides is fewer than beta.')

    if args.beta > 0 and args.beta < args.multiplicity:
        print(f'{bcolors.RED}> Warning{bcolors.RESET}: Beta is set to {args.beta}, a positive value smaller than the multiplicity {args.multiplicity}')
        print(f'{bcolors.BLUE}>{bcolors.RESET} {bcolors.ORANGE}ALLEGRO{bcolors.RESET} cannot find a total of {args.beta} guides while each guide container is required to be targeted at least {args.multiplicity} times.')
        print(f'{bcolors.BLUE}>{bcolors.RESET} Auto adjusting beta to be equal to the multiplicity. You may also set beta to 0. See the manual for more details.')
        args.beta = args.multiplicity

    if args.filter_repetitive == True:
        print(f'{bcolors.BLUE}>{bcolors.RESET} filter_repetitive is set to True. Filtering guides with repetitive sequences.')

    if args.num_trials <= 0:
        print(f'{bcolors.BLUE}>{bcolors.RESET} num_trials is <= 0. Running randomized rounding only once. Note that the solution may not be the one with the smallest size.')
        args.num_trials = 1

    scorer_settings = configure_scorer_settings(args)

    # Create the output folder using the output directory and experiment name
    args.output_directory = create_output_directory(args.output_directory, args.experiment_name)

    return args, scorer_settings


def log_args(args: argparse.Namespace) -> None:
    output_txt_path = os.path.join(args.output_directory, args.experiment_name + '_config_used.txt')

    with open(output_txt_path, 'w') as f:
        f.write(f'Config used for experiment {args.experiment_name}\n')
        for key, value in vars(args).items():
            f.writelines(key + ': ' + str(value) + '\n')


def create_output_directory(output_directory: str, experiment_name: str) -> str:
    dir_name = os.path.join(output_directory, experiment_name)

    # Check if the directory already exists
    if os.path.exists(dir_name):
        # If it exists, append a number to the directory name
        i = 1

        while os.path.exists(dir_name + '_' + str(i)):
            i += 1

        dir_name = dir_name + '_' + str(i)

    print(f'{bcolors.BLUE}>{bcolors.RESET} Creating directory', dir_name)
    os.makedirs(dir_name)

    return dir_name


def configure_scorer_settings(args: argparse.Namespace) -> dict:
    scorer_settings = dict()

    match args.scorer:
            case 'chopchop':
                scorer_settings = {
                    'experiment_name': args.experiment_name,
                    'output_directory': args.output_directory,
                    'chopchop_scoring_method': args.chopchop_scoring_method,
                    'absolute_path_to_chopchop': args.absolute_path_to_chopchop,
                    'absolute_path_to_genomes_directory': args.absolute_path_to_genomes_directory,
                    'input_species_csv_file_path': args.input_species_path
                }

            case 'dummy':
                scorer_settings = {
                    'pam': 'NGG',
                    'protospacer_length': 20,
                    'filter_repetitive': args.filter_repetitive,
                    'context_toward_five_prime': 0,
                    'context_toward_three_prime': 0,
                }

            case 'uCRISPR' | 'ucrispr':
                scorer_settings = {
                    'pam': 'NGG',
                    'use_secondary_memory': args.use_secondary_memory,
                    'protospacer_length': 20,
                    'filter_repetitive': args.filter_repetitive,
                    'context_toward_five_prime': 4,  # example: ACAATTTAAAGCTTGCCTCTAACTTGGCCA
                    'context_toward_three_prime': 3,  # 4 refers to ACAA, followed by a 20mer,
                                                        # then TGG PAM, then 3 bases CCA 
                }

            case _:
                print(f'{bcolors.RED}> Warning{bcolors.RESET}: Unknown scorer {args.scorer} selected in config.yaml. Exiting.')
                sys.exit(1)
    
    return scorer_settings