import os
import sys
import yaml
import numpy
import pandas
import argparse
import matplotlib.pyplot

from utils.guide_finder import GuideFinder
from coverset_parsers.coverset import CoversetsCython as coverset # type: ignore

matplotlib.pyplot.rcParams['figure.dpi'] = 300


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
    
    parser.set_defaults(**config_arg_dict)
    args = parser.parse_args(remaining_args)

    return args


# def validate_arguments(args) -> None:
# TODO: Implement 


def graph_size_dist(
    beta: int,
    exp_name: str, 
    output_dir: str,
    size_of_solutions_for_n_trials: list[int],
    ) -> None:

    print('Drawing size distribution graph...')
    output_path = os.path.join(output_dir, exp_name + '_b' + str(beta) + '_size_hist.png')

    y = numpy.asarray(size_of_solutions_for_n_trials)
    average_size_over_all_trials = numpy.mean(y)

    n = len(size_of_solutions_for_n_trials)
    x = numpy.arange(n)

    title_1 = 'Set size distribution for {n} trials'.format(n=len(size_of_solutions_for_n_trials))
    title_2 = 'Average size over all trials: {n:.2f}'.format(n=average_size_over_all_trials)
    title_3 = 'Beta: {b}'.format(b=beta)

    mininum_size = min(size_of_solutions_for_n_trials)
    red_label = 'Smallest size: {n}'.format(n=mininum_size)
    blue_mask = y != mininum_size
    red_mask = y == mininum_size

    matplotlib.pyplot.figure(figsize=(20, 3))

    matplotlib.pyplot.bar(x=x[blue_mask], height=y[blue_mask])
    matplotlib.pyplot.bar(x=x[red_mask], height=y[red_mask], color='red', label=red_label)

    matplotlib.pyplot.legend()
    matplotlib.pyplot.title(title_1 + '\n' + title_2 + '\n' + title_3)
    matplotlib.pyplot.xlabel('Trial')
    matplotlib.pyplot.ylabel('Set size')
    matplotlib.pyplot.margins()

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(output_path)


def graph_score_dist(
    beta: int,
    exp_name: str, 
    output_dir: str,
    average_scores_for_n_trials: list[float],
    ) -> None:

    print('Drawing score distribution graph...')
    output_path = os.path.join(output_dir, exp_name + '_b' + str(beta) + '_avg_score_hist.png')

    title = 'Solution guides\' average score distribution over {n} trials'.format(n=len(average_scores_for_n_trials))

    matplotlib.pyplot.figure(figsize=(10, 3))

    matplotlib.pyplot.hist(average_scores_for_n_trials, edgecolor='black', linewidth=1.2)
    
    matplotlib.pyplot.title(title)
    matplotlib.pyplot.xlabel('Average score')
    matplotlib.pyplot.ylabel('Count')

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(output_path)


def write_solution_to_file(
    beta: int,
    species_names: list[str],
    solution: list[tuple[str, str]],
    experiment_name: str,
    input_csv_path: str,
    input_sequence_directory: str,
    paths_csv_column_name: str,
    species_names_csv_column_name: str,
    output_directory: str,
    ) -> None:

    # TODO maybe move to own module -- main.py is too crowded...
    output_txt_path = os.path.join(output_directory, experiment_name + '_b{b}.txt'.format(b=beta))
    output_csv_path = os.path.join(output_directory, experiment_name + '_b{b}.csv'.format(b=beta))

    i = 1  # If an output file with the same name already exists, make a new numbered one.
    while os.path.isfile(output_txt_path):
        new_experiment_name = experiment_name + '_b{b}_'.format(b=beta) + str(i)
        output_txt_path = os.path.join(output_directory, new_experiment_name + '.txt')
        output_csv_path = os.path.join(output_directory, new_experiment_name + '.csv')
        i += 1


    # Writing the new output file.
    print('Writing to file:', output_txt_path)
    with open(output_txt_path, 'w') as f:
        f.write('We can cut the following {n} species: {species}.\n'.format(
            n=len(species_names),
            species=str(species_names),
            )
        )

        f.write('Using the following {n} guides: {guides}.\n'.format(
            n=str(len(solution)),
            guides=str(solution),
            )
        )

    paths: list[str] = list()
    species_list: list[str] = list()
    scores: list[int] = list()
    strands: list[str] = list()
    sequences: list[str] = list()
    end_positions: list[int] = list()
    start_positions: list[int] = list()
    chromosomes_or_genes: list[str] = list()

    guide_finder = GuideFinder()
    input_df = pandas.read_csv(input_csv_path)[[species_names_csv_column_name, paths_csv_column_name]]

    for pair in solution:
        seq = pair[0]
        hit_species = pair[1]

        for species in hit_species:
            df_file_path = input_df[input_df[species_names_csv_column_name] == species][paths_csv_column_name].values[0]
            full_path = os.path.join(input_sequence_directory, df_file_path)

            list_of_tuples = guide_finder.locate_guides_in_sequence(sequence=seq, file_path=full_path, to_upper=True)

            for tuple in list_of_tuples:
                chromosome = tuple[0]
                strand = tuple[1]
                start_pos = tuple[2]
                end_pos = tuple[3]

                sequences.append(seq)
                paths.append(df_file_path)
                scores.append(1)  # TODO fix for other than 1
                strands.append(strand)
                start_positions.append(start_pos)
                end_positions.append(end_pos)
                chromosomes_or_genes.append(chromosome)
                species_list.append(species)
            
    
    pandas.DataFrame(list(zip(
        sequences,
        species_list,
        scores,
        chromosomes_or_genes,
        strands,
        start_positions,
        end_positions,
        paths,
        )),
    columns=['sequence','targets', 'score', 'chromosome_or_gene',
    'strand', 'start_position', 'end_position', 'path']).to_csv(output_csv_path, index=False)
    
    print('Done. Check {path} for the output.'.format(path=output_csv_path))


def output_csv(
    beta: int,
    experiment_name: str,
    output_directory: str, 
    solver,
    ) -> None:

    output_path = os.path.join(output_directory, experiment_name + '_metrics.csv'.format(b=beta))

    print('Writing to file:', output_path)

    avg_num_while_iters_for_n_trials = 0
    if len(solver.num_while_iters_for_each_trial) > 0:
        avg_num_while_iters_for_n_trials = numpy.mean(solver.num_while_iters_for_each_trial).round(2)

    df = pandas.DataFrame({
        'budget': [beta],
        'avg_score_over_all_trials': [solver.average_score_for_all_trials],
        'avg_score_over_each_trial': [solver.average_score_for_each_trial],
        'avg_num_while_iters_for_n_trials': [avg_num_while_iters_for_n_trials],
        'num_while_iters_for_each_trial': [solver.num_while_iters_for_each_trial],
        'fractional_glop_vals': [solver.fractional_glop_vals],
        'solver_time': [solver.solver_time],
        'exhaustive_time': [solver.exhaustive_time],
        'random_rounding_time': [solver.randomized_rounding_time],
        'num_non_zero_feasible': [solver.num_non_zero_feasible],
        'num_non_zero_feasible_under_one': [solver.num_feasible_guides_with_prob_lt_one],
        'num_exhausted_combos': [solver.num_exhausted_combos],
        'solved_with_exhaustive': ['Yes' if solver.solved_with_exhaustive else 'No'],
        'size_of_solutions_for_n_trials': [solver.set_size_for_each_trial],
        'num_species_constraints': [len(solver.species)],
        'num_guide_variables': [solver.k],
        'num_rounding_trials': [solver.num_trials],
    })

    write_header_if_file_doesnt_exist = not os.path.exists(output_path)  # otherwise append w/ mode='a'

    df.to_csv(output_path, mode='a', header=write_header_if_file_doesnt_exist, index=False)
    print('Done. Check {path} for the output.'.format(path=output_path))


def main() -> int:
    print('Welcome to ALLEGRO. All unspecified command-line arguments default to the values in config.yaml')
    
    args = parse_arguments()
    scorer_settings = dict()

    match args.scorer:
        case 'chopchop':
            scorer_settings = {
                'output_directory': args.output_directory,
                'chopchop_scoring_method': args.chopchop_scoring_method,
                'absolute_path_to_chopchop': args.absolute_path_to_chopchop,
                'absolute_path_to_genomes_directory': args.absolute_path_to_genomes_directory,
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

    coversets_obj = coverset(
        num_trials=args.num_trials,
        cas_variant=args.cas,
        guide_source=args.mode,
        guide_length=args.protospacer_length,
        scorer_name=args.scorer,
        scorer_settings=scorer_settings,
        input_cds_directory=args.input_cds_directory,
        input_species_csv_file_path=args.input_species_path,
        input_genome_directory=args.input_genomes_directory,
    )

    # solver = Solver(
    #     beta=args.beta,
    #     objective=args.objective,
    #     species=coversets_obj.species_set,
    #     coversets=coversets_obj.coversets,
    #     num_trials=args.num_trials,
    #     solver_engine=args.solver_engine,
    #     exhaustive_threshold=args.exhaustive_threshold,
    # )

    # solution = solver.solve()

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
    
    write_solution_to_file(
        beta=args.beta,
        species_names=coversets_obj.species_names,
        solution=coversets_obj.solution,
        experiment_name=args.experiment_name,
        input_csv_path=args.input_species_path,
        input_sequence_directory=args.input_genomes_directory,
        paths_csv_column_name='genome_file_name',
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