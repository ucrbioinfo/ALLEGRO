import os
import sys
import yaml
import numpy
import pandas
import argparse
import matplotlib.pyplot

from solvers.solver import Solver
from cover_set_parsers.coversets import Coversets


matplotlib.pyplot.rcParams['figure.dpi'] = 300


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--config',
        type=argparse.FileType(mode='r'),
        default='config.yaml', 
        help='The config file to use. Must be placed in the root folder.',
    )

    args = parser.parse_args()
    if args.config:
        data = yaml.load(args.config, Loader=yaml.FullLoader)
        arg_dict = args.__dict__

        for key, value in data.items():
            arg_dict[key] = value

    return args


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


def print_solution(solution: set[str], parser: Coversets) -> None:
    for guide_seq in solution:
        guide_objects = parser.get_guides_from_seq(guide_seq)
        
        for guide_object in guide_objects:
            guide_object.print_info()
            print()


def write_solution_to_file(
    beta: int,
    parser: Coversets,
    solution: set[str], 
    experiment_name: str,
    output_directory: str,
    ) -> None:

    output_path = os.path.join(output_directory, experiment_name + '_b{b}.txt'.format(b=beta))
    output_csv_path = os.path.join(output_directory, experiment_name + '_b{b}.txt'.format(b=beta))


    i = 1  # If an output file with the same name already exists, make a new numbered one.
    while os.path.isfile(output_path):
        new_experiment_name = experiment_name + '_b{b}_'.format(b=beta) + str(i)
        output_path = os.path.join(output_directory, new_experiment_name + '.txt')
        output_csv_path = os.path.join(output_directory, new_experiment_name + '.csv')
        i += 1


    # Writing the new output file.
    print('Writing to file:', output_path)
    with open(output_path, 'w') as f:
        f.write('We can cut the following {n} species: {species}.\n'.format(
            n=len(parser.species_names),
            species=str(parser.species_names),
            )
        )

        f.write('Using the following {n} guides: {guides}.\n'.format(
            n=str(len(solution)),
            guides=str(solution),
            )
        )

    list_of_attributes_dicts: list[dict] = list()

    for guide_seq in solution:
        guide_objects = parser.get_guides_from_seq(guide_seq)
        
        for guide_object in guide_objects:
            list_of_attributes_dicts.append(guide_object.get_attributes_dict())

    aggregate_dict = dict()
    for key in list_of_attributes_dicts[0].keys():
        aggregate_dict[key] = [d[key] for d in list_of_attributes_dicts]    
    
    pandas.DataFrame.from_dict(aggregate_dict).to_csv(output_csv_path)
    
    print('Done. Check {path} for the output.'.format(path=output_csv_path))


def output_csv(
    beta: int,
    experiment_name: str,
    output_directory: str, 
    solver,
    ) -> None:
    
    output_path = os.path.join(output_directory, experiment_name + '_metrics.csv'.format(b=beta))

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


def main() -> int:
    args = parse_arguments()

    scorer_settings = dict()
    match args.scorer:
        case 'chopchop':
            scorer_settings = {
                'output_directory': args.output_directory,
                'chopchop_scoring_method': args.chopchop_scoring_method,
                'absolute_path_to_chopchop': args.absolute_path_to_chopchop,
                'absolute_path_to_bowtie_build': args.absolute_path_to_bowtie_build,
                'absolute_path_to_genomes_directory': args.absolute_path_to_genomes_directory,
            }

    coversets = Coversets(
        guide_source=args.mode,
        cas_variant = args.cas,
        scorer_name=args.scorer,
        scorer_settings=scorer_settings,
        input_cds_directory=args.input_cds_directory,
        input_species_csv_file_path=args.input_species_path,
        input_genome_directory=args.input_genomes_directory,
    )

    solver = Solver(
        beta=args.beta,
        objective=args.objective,
        coverset_parser=coversets,
        num_trials=args.num_trials,
        solver_engine=args.solver_engine,
        exhaustive_threshold=args.exhaustive_threshold,
    )

    solution = solver.solve()

    if args.traceback:
        print_solution(
            parser=coversets,
            solution=solution,
        )

    if args.graph and len(solution) > 0 and not solver.solved_with_exhaustive:
        graph_size_dist(
            beta=solver.beta,
            exp_name=args.experiment_name, 
            output_dir=args.output_directory,
            size_of_solutions_for_n_trials=solver.set_size_for_each_trial,
        )

        graph_score_dist(
            beta=solver.beta,
            exp_name=args.experiment_name, 
            output_dir=args.output_directory,
            average_scores_for_n_trials=solver.average_score_for_each_trial,
        )
    
    write_solution_to_file(
        beta=args.beta,
        parser=coversets,
        solution=solution,
        experiment_name=args.experiment_name,
        output_directory=args.output_directory,
    )

    if args.output_csv:
        output_csv(
            solver=solver,
            beta=args.beta,
            experiment_name=args.experiment_name,
            output_directory=args.output_directory,
        )

    return 0


if __name__ == '__main__':
    sys.exit(main())