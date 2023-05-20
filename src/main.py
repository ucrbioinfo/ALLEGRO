import os
import sys
import yaml
import argparse

import utils.write_solution_to_file as write_solution
from cython_libs.coverset2 import CoversetsCython as coverset  # type: ignore

def parse_arguments() -> argparse.Namespace:
    config_parser = argparse.ArgumentParser(add_help=False)

    config_parser.add_argument(
        '-c',
        '--config',
        type=argparse.FileType(mode='r'),
        default='config.yaml',
        help='The config file to use. Must be placed in the root folder.',
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

    parser.add_argument(
        '-n',
        '--experiment_name',
        type=str,
        help='Name of the experiment. Output file(s) will be labeled with this.'
    )

    # parser.add_argument(
    #     '-m',
    #     '--mode',
    #     type=str,
    #     help="'from_genome' or 'from_orthogroups'",
    # )

    parser.add_argument(
        '--objective',
        type=str,
        help=("'min' or 'max'. max and uses the --beta argument. " +
        "Select 'min' to disregard beta and find the smallest set of guides to " +
        "cut all species with."),
    )

    parser.add_argument(
        '--scorer',
        type=str,
        help=("Which scoring method to use? Options are 'chopchop', 'dummy' " +
        "where 'dummy' assigns a score of 1.0 to all guides."),
    )

    parser.add_argument(
        '--solver_engine',
        type=str,
        help=("Which Google OR-Tools pywraplp solver to use? For a list see " +
        "here https://developers.google.com/optimization/lp/lp_advanced"),
    )

    parser.add_argument(
        '--cas',
        type=str,
        default='cas9',
        help="Which cas endonuclease to use? Defaults to 'cas9'. Options are: 'cas9'.",
    )

    parser.add_argument(
        '-b',
        '--beta',
        type=int,
        help=("Beta represents a loose budge or threshold for the maximum " + 
        "number of guides you would like in the final covering set. This is " +
        "NOT a guaranteed maximum due to the hardness of the problem. See the " +
        "topic on Integrality Gap."),
    )

    parser.add_argument(
        '--input_species_path',
        type=str,
        help=(".csv file that includes three columns: species_name, " +
        "genome_file_name, cds_file_name. genome_file_name rows start " +
        "with species_name + _genomic.fna. cds_file_name rows start with " +
        "species_name + _cds.fna. Genome files must be placed in data/input/genomes " +
        "while cds files must be placed in data/input/cds. These directories may be changed " +
        "with other config options."),
    )

    parser.add_argument(
        '--input_genomes_directory',
        type=str,
    )

    parser.add_argument(
        '--output_directory',
        type=str,
    )

    parser.add_argument(
        '--input_cds_directory',
        type=str,
        help="Files in this directory must end with _cds.fna or _cds.faa or _cds.fasta",
    )

    # parser.add_argument(
    #     '--orthogroups_path',
    #     type=str,
    # )

    parser.add_argument(
        '--chopchop_scoring_method',
        type=str,
        help='Only used in chopchop is selected as the scorer. Which chopchop scoring model to use?',
    )

    parser.add_argument(
        '--absolute_path_to_chopchop',
        type=str,
        help=("Point to the directory where chopchop is located. Download " +
        "chopchop from here https://bitbucket.org/valenlab/chopchop/src/master/"),
    )

    parser.add_argument(
        '--absolute_path_to_genomes_directory',
        type=str,
        help="Points to the directory where all species' genomes are placed.",
    )

    parser.add_argument(
        '--exhaustive_threshold',
        type=int,
        help=("Search all combinatorial possibilities unless there are more " +
        "feasible solutions than this (then uses randomized rounding). The number " +
        "of calculations grows exponentially. ONLY increase this parameter when the " +
        "number of guides in under about twenty or you will have to wait years for it to complete :)."),
    )

    parser.add_argument(
        '--num_trials',
        type=int,
        help=("Only used if the number of feasible guides is above the exhaustive_threshold. " +
        "How many times to run the randomized rounding algorithm?"),
    )

    parser.add_argument(
        '--mp_threshold',
        type=int,
        help=("A higher number increases running time while decreasing memory consumption." +
        " Pre-select guides that hit only up to this number of species to act as representatives for these species."),
    )
    
    parser.set_defaults(**config_arg_dict)
    args = parser.parse_args(remaining_args)

    return args


def check_and_fix_configurations(args: argparse.Namespace) -> argparse.Namespace:
    if args.multiplicity < 1:
        print('WARNING: Multiplicity is set to {m}, a value smaller than 1. Auto adjusting multiplicity to 1.')
        args.multiplicity = 1

    if args.mp_threshold <= 0:
        print('mp_threshold is set to {mp} and thus disabled. Saving all guides to memory.'.format(
            mp=args.mp_threshold
        ))
        args.mp_threshold = 0

    if args.mp_threshold > 0 and args.mp_threshold < args.multiplicity:
        print('WARNING: mp_threshold is set to {mp}, a non-zero value smaller than the multiplicity {mult}.'.format(
            mp=args.mp_threshold,
            mult=args.multiplicity))
        print('ALLEGRO cannot remove all but {mp} guides from each container and still ensure each container is targeted at least {m} times.'.format(
            mp=args.mp_threshold,
            mult=args.multiplicity))
        print('Auto adjusting mp_threshold to be equal to multiplicity. You may also set mp_threshold to 0 to disable this memory-saving feature.')
        args.mp_threshold = args.multiplicity

    if args.beta <= 0:
        args.beta = 0
        print('Beta is set to {b} and thus disabled.'.format(b=args.beta))

        if args.scorer != 'dummy':
            print('ALLEGRO will find the guides with the best efficiency. It will not be minimizing the set size.')
        else:
            print('ALLEGRO will minimize the set size.')
        
    if args.scorer == 'dummy' and args.beta > 0:
        # No feasible solutions if there are fewer guides than beta
        # Say there are 5 species, 5 guides total, and beta is set to 1. Say that none of the species share any guides.
        # This will ask ALLEGRO to find 1 guide out of 5 to cover all 5 species. There is no solution.
        print('The scorer is set to dummy and beta to non-zero {b}. ALLEGRO will try to find (approximately) {b} guides to cover all guide containers.'.format(b=args.beta))
        print('WARNING: ALLEGRO may find that there is no feasible solution if the number of shared guides is fewer than beta.')

    if args.beta > 0 and args.beta < args.multiplicity:
        print('WARNING: Beta is set to {b}, a non-zero value smaller than the multiplicity {mp}'.format(
            b=args.beta,
            mp=args.multiplicity
        ))
        print('ALLEGRO cannot find a total of {b} guides while each guide container is required to be targeted at least {m} times.'.format(
            b=args.beta,
            mp=args.multiplicity
        ))
        print('Auto adjusting beta to be equal to the multiplicity. You may also set beta to 0. See the documentation for more details.')
        args.beta = args.multiplicity

    if args.include_repetitive == True:
        print('include_repetitive is set to True. Filtering guides with repetitive sequences.')

    if args.num_trials < 0:
        print('num_trials is set to 0. Running randomized rounding only once. Note that the solution may not be the one with the smallest size.')
        args.num_trials = 0

    return args


def log_args(args: argparse.Namespace) -> None:
    output_txt_path = os.path.join(args.output_directory, args.experiment_name + '_config_used.txt')

    with open(output_txt_path, 'w') as f:
        f.write('Config used for experiment ' + args.experiment_name + '\n')
        for key, value in vars(args).items():
            f.writelines(key + ': ' + str(value) + '\n')


def create_output_directory(output_directory: str, experiment_name: str) -> str:
    dir_name = os.path.join(output_directory, experiment_name)

    # Check if the directory already exists
    if os.path.exists(dir_name):
        # If it exists, append a number to the directory name
        i = 1

        while os.path.exists(dir_name + "_" + str(i)):
            i += 1

        dir_name = dir_name + "_" + str(i)

    print('Creating directory', dir_name)
    os.makedirs(dir_name)

    return dir_name


def make_scorer_settings(args: argparse.Namespace):
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
                    'pam': args.pam,
                    'protospacer_length': args.protospacer_length,
                    'include_repetitive': args.include_repetitive,
                    'context_toward_five_prime': args.context_toward_five_prime,
                    'context_toward_three_prime': args.context_toward_three_prime,
                }

            case _:
                print('Unknown scorer selected. Aborting.')
                raise ValueError
    
    return scorer_settings


def main() -> int:
    print('Welcome to ALLEGRO. All unspecified command-line arguments default to the values in config.yaml')
    
    args = parse_arguments()
    args = check_and_fix_configurations(args)
    args.output_directory = create_output_directory(args.output_directory, args.experiment_name)
    log_args(args)

    scorer_settings = make_scorer_settings(args)

    # coversets_obj = coverset(
    #     beta=args.beta,
    #     num_trials=args.num_trials,
    #     cas_variant=args.cas,
    #     guide_source=args.mode,
    #     guide_length=args.protospacer_length,
    #     scorer_name=args.scorer,
    #     scorer_settings=scorer_settings,
    #     monophonic_threshold=args.mp_threshold,
    #     output_directory=args.output_directory,
    #     input_cds_directory=args.input_cds_directory,
    #     input_species_csv_file_path=args.input_species_path,
    #     input_genome_directory=args.input_genomes_directory,
    # )

    coversets_obj = coverset(
        num_trials=args.num_trials,
        cas_variant=args.cas,
        guide_source=args.mode,
        guide_length=args.protospacer_length,
        scorer_name=args.scorer,
        scorer_settings=scorer_settings,
        mp_threshold=args.mp_threshold,
        beta=args.beta,
        gene_cut_multiplicity=args.multiplicity,
        output_directory=args.output_directory,
        input_cds_directory=args.input_cds_directory,
        input_species_csv_file_path=args.input_species_path,
        input_genome_directory=args.input_genomes_directory,
    )

    # if args.graph and len(solution) > 0 and not solver.solved_with_exhaustive:
    #     graph_size_dist(
    #         beta=solver.beta,
    #         exp_name=args.experiment_name, 
    #         output_dir=args.output_directory,
    #         size_of_solutions_for_n_trials=solver.set_size_for_each_trial,
    #     )

    #     graph_score_dist(
    #         beta=solver.beta,
    #         exp_name=args.experiment_name, 
    #         output_dir=args.output_directory,
    #         average_scores_for_n_trials=solver.average_score_for_each_trial,
    #     )

    d = {
        'input_sequence_directory': '',
        'paths_csv_column_name': ''
    }
    match args.mode:
        case 'from_genome':
            d['input_sequence_directory'] = args.input_genomes_directory
            d['paths_csv_column_name'] = 'genome_file_name'

        case 'from_cds':
            d['input_sequence_directory'] = args.input_cds_directory
            d['paths_csv_column_name'] = 'cds_file_name'
        
        case _:
            print('Unknown mode selected. Aborting.')
            raise ValueError
    

    if args.mode == 'from_genome':
        write_solution.write_solution_to_file(
            beta=args.beta,
            species_names=coversets_obj.species_names,
            solution=coversets_obj.solution,
            experiment_name=args.experiment_name,
            input_csv_path=args.input_species_path,
            input_sequence_directory=d['input_sequence_directory'],
            paths_csv_column_name=d['paths_csv_column_name'],
            species_names_csv_column_name='species_name',
            output_directory=args.output_directory
        )
    elif args.mode == 'from_cds':
            write_solution.write_cds_solution_to_file(
            species_names=coversets_obj.species_names,
            gene_names=['LYS2', 'LYS5', 'MET17', 'TRP1', 'URA3', 'FCY1', 'GAP1', 'CAN1'],
            multiplicity=args.multiplicity,
            solution=coversets_obj.solution,
            experiment_name=args.experiment_name,
            input_csv_path=args.input_species_path,
            input_sequence_directory=d['input_sequence_directory'],
            paths_csv_column_name=d['paths_csv_column_name'],
            species_names_csv_column_name='species_name',
            output_directory=args.output_directory
        )

    # if args.output_csv:
    #     output_csv(
    #         solver=solver,
    #         beta=args.beta,
    #         experiment_name=args.experiment_name,
    #         output_directory=args.output_directory,
    #     )

    return 0


if __name__ == '__main__':
    sys.exit(main())