# Utility functions imported by main.py.
import os
import shutil
import pandas
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed

from utils.shell_colors import bcolors
from utils.guide_finder import align_guides_to_seq_bowtie


def process_row(row, species_names_csv_column_name, paths_csv_column_name, input_directory, output_align_temp_path, output_bowtie_temp_dir, solution_dict, pam):
    name = row[species_names_csv_column_name]
    df_file_path = row[paths_csv_column_name]
    full_path = os.path.join(input_directory, df_file_path)

    (bowtie_sequences,
     bowtie_strands,
     reference_names,
     orthos,
     bowtie_start_positions,
     time_elapsed) = align_guides_to_seq_bowtie(
        name=name,
        output_align_temp_path=output_align_temp_path,
        file_path=full_path,
        output_directory=output_bowtie_temp_dir)

    bowtie_scores = [solution_dict[s[:-len(pam)]][0] for s in bowtie_sequences]

    return {
        "sequences": bowtie_sequences,
        "paths": [df_file_path] * len(bowtie_sequences),
        "scores": bowtie_scores,
        "strands": bowtie_strands,
        "start_positions": bowtie_start_positions,
        "chromosomes_or_genes": reference_names,
        "species_list": [name] * len(bowtie_sequences),
        "ortho_list": orthos,
        "time_elapsed": time_elapsed
    }

def process_dataframe(input_df, species_names_csv_column_name, paths_csv_column_name, input_directory, output_align_temp_path, output_bowtie_temp_dir, solution_dict, pam):
    count = 0
    num_cpus = os.cpu_count()
    max_workers = num_cpus - 1 if num_cpus is not None else 1  # Ensure at least 1 worker is used
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Prepare the partial function to include all necessary parameters except 'row'
        func = partial(process_row,
                       species_names_csv_column_name=species_names_csv_column_name,
                       paths_csv_column_name=paths_csv_column_name,
                       input_directory=input_directory,
                       output_align_temp_path=output_align_temp_path,
                       output_bowtie_temp_dir=output_bowtie_temp_dir,
                       solution_dict=solution_dict,
                       pam=pam)
        
        # Submit tasks
        futures = [executor.submit(func, row) for _, row in input_df.iterrows()]

        # Initialize lists for results outside the loop
        sequences = []
        paths = []
        scores = []
        strands = []
        start_positions = []
        chromosomes_or_genes = []
        species_list = []
        ortho_list = []
        time_elapsed = 0.0

        # Collect results as tasks complete
        for future in as_completed(futures):
            result = future.result()
            sequences.extend(result["sequences"])
            paths.extend(result["paths"])
            scores.extend(result["scores"])
            strands.extend(result["strands"])
            start_positions.extend(result["start_positions"])
            chromosomes_or_genes.extend(result["chromosomes_or_genes"])
            species_list.extend(result["species_list"])
            ortho_list.extend(result["ortho_list"])
            time_elapsed += result["time_elapsed"]
            count += 1
            print(f'{bcolors.BLUE}>{bcolors.RESET} Done with {count + 1}/{len(input_df)} species...', end='\r')

        return sequences, paths, scores, strands, start_positions, chromosomes_or_genes, species_list, ortho_list, time_elapsed


def write_solution_to_file_bowtie(
    pam: str,
    species_names: list[str],
    solution: list[tuple[str, float, list[str]]],
    experiment_name: str,
    input_csv_path: str,
    input_directory: str,
    paths_csv_column_name: str,
    species_names_csv_column_name: str,
    output_directory: str,
    ) -> tuple[str, float]:

    output_txt_path = os.path.join(output_directory, experiment_name + '.txt')
    output_library_path = os.path.join(output_directory, experiment_name + '_library.txt')
    output_csv_path = os.path.join(output_directory, experiment_name + '.csv')

    output_bowtie_temp_dir = os.path.join(output_directory, 'temp')
    output_align_temp_path = os.path.join(output_directory, 'temp', 'temp_library.txt')

    if not os.path.exists(output_bowtie_temp_dir):
        os.mkdir(output_bowtie_temp_dir)

    with open(output_txt_path, 'w') as f:

        for tup in solution:
            seq = tup[0]
            score = tup[1]
            hit_species = tup[2]

            f.write(f'Guide {seq} targets {len(hit_species)} locations.\n')

        f.write(f'We can cut the following {len(species_names)} species: {str(species_names)}.\n')
        f.write(f'Using the following {str(len(solution))} guides: {str(solution)}.\n')

    # guide_finder = GuideFinder()
    input_df = pandas.read_csv(input_csv_path)[[species_names_csv_column_name, paths_csv_column_name]]
    
    solution_dict = dict()

    for tup in solution:
        seq = tup[0]
        score = tup[1]
        hit_species = tup[2]

        if seq not in solution_dict:
            solution_dict[seq] = (score, hit_species)
        else:
            print('Dev error in utils/write_solution_to_file.py. Not your fault. Just bad coding.')

    library = list(solution_dict.keys())
    with open(output_align_temp_path, 'w') as f:
        for g in library:
            for n in ['AGG', 'CGG', 'TGG', 'GGG']:
                f.write(f'>{g}{n}\n{g}{n}\n')

    with open(output_library_path, 'w') as f:
        for g in library:
            f.write(f'{g}\n')
    print(f'{bcolors.BLUE}>{bcolors.RESET} Your library is ready. Find it in {output_library_path}.')

    total_time_elapsed = 0.0

    print(f'{bcolors.BLUE}>{bcolors.RESET} Aligning the library against the input species to generate the final report...')

    # paths: list[str] = list()
    # species_list: list[str] = list()
    # scores: list[float] = list()
    # strands: list[str] = list()
    # sequences: list[str] = list()
    # start_positions: list[int] = list()
    # chromosomes_or_genes: list[str] = list()
    # ortho_list: list[str] = list()

    # for idx, row in input_df.iterrows():
    #     name = row[species_names_csv_column_name]
    #     df_file_path = row[paths_csv_column_name]
    #     full_path = os.path.join(input_directory, df_file_path)

    #     (bowtie_sequences,
    #      bowtie_strands,
    #      reference_names,
    #      orthos,
    #      bowtie_start_positions,
    #      time_elapsed) = align_guides_to_seq_bowtie(
    #         name=name,
    #         output_align_temp_path=output_align_temp_path,
    #         file_path=full_path,
    #         output_directory=output_bowtie_temp_dir)
        
    #     total_time_elapsed += time_elapsed

    #     bowtie_scores = [solution_dict[s[:-len(pam)]][0] for s in bowtie_sequences]

    #     sequences.extend(bowtie_sequences)
    #     paths.extend([df_file_path] * len(bowtie_sequences))
    #     scores.extend(bowtie_scores)
    #     strands.extend(bowtie_strands)
    #     start_positions.extend(bowtie_start_positions)
    #     chromosomes_or_genes.extend(reference_names)
    #     species_list.extend([name] * len(bowtie_sequences))
    #     ortho_list.extend(orthos)

    #     print(f'{bcolors.BLUE}>{bcolors.RESET} Done with {idx + 1}/{len(input_df)} species...', end='\r')
    # print()

    (sequences,
    paths,
    scores,
    strands,
    start_positions,
    chromosomes_or_genes,
    species_list,
    ortho_list,
    time_elapsed) = process_dataframe(input_df,
                                      species_names_csv_column_name,
                                      paths_csv_column_name,
                                      input_directory,
                                      output_align_temp_path,
                                      output_bowtie_temp_dir,
                                      solution_dict,
                                      pam)
    
    total_time_elapsed += time_elapsed

    pandas.DataFrame(list(zip(
        sequences,
        species_list,
        chromosomes_or_genes,
        ortho_list,
        strands,
        start_positions,
        scores,
        paths,
        )),
    columns=['sequence', 'target', 'reference_name', 'orthologous_to',
    'strand', 'start_position', 'score', 'path']).to_csv(output_csv_path, index=False)
    
    print(f'{bcolors.BLUE}>{bcolors.RESET} Done. Check {output_csv_path} for the output.')

    shutil.rmtree(output_bowtie_temp_dir)

    return output_csv_path, total_time_elapsed

# DEPRECATED
# def write_solution_to_file(
#     pam: str,
#     species_names: list[str],
#     solution: list[tuple[str, float, str]],
#     experiment_name: str,
#     input_csv_path: str,
#     input_directory: str,
#     paths_csv_column_name: str,
#     species_names_csv_column_name: str,
#     output_directory: str,
#     ) -> str:

#     output_txt_path = os.path.join(output_directory, experiment_name + '.txt')
#     output_csv_path = os.path.join(output_directory, experiment_name + '.csv')

#     print(f'{bcolors.BLUE}>{bcolors.RESET} Writing to file: {output_txt_path}')
#     with open(output_txt_path, 'w') as f:

#         for tup in solution:
#             seq = tup[0]
#             score = tup[1]
#             hit_species = tup[2]

#             f.write(f'Guide {seq} targets {len(hit_species)} genes/chromosomes.\n')
#         f.write(f'We can cut the following {len(species_names)} species: {str(species_names)}.\n')
#         f.write(f'Using the following {str(len(solution))} guides: {str(solution)}.\n')

#     paths: list[str] = list()
#     misc_list: list[str] = list()
#     species_list: list[str] = list()
#     scores: list[float] = list()
#     strands: list[str] = list()
#     sequences: list[str] = list()
#     end_positions: list[int] = list()
#     start_positions: list[int] = list()
#     chromosomes_or_genes: list[str] = list()
#     ortho_list: list[str] = list()

#     guide_finder = GuideFinder()
#     input_df = pandas.read_csv(input_csv_path)[[species_names_csv_column_name, paths_csv_column_name]]

#     print(f'{bcolors.BLUE}>{bcolors.RESET} Locating the solution guides in the input species to generate the final report...')
    
#     for tup in solution:
#         seq = tup[0]
#         score = tup[1]
#         hit_species = tup[2]

#         for species in hit_species:
#             df_file_path = input_df[input_df[species_names_csv_column_name] == species][paths_csv_column_name].values[0]
#             full_path = os.path.join(input_directory, df_file_path)

#             list_of_tuples = guide_finder.locate_guides_in_sequence(pam=pam, sequence=seq, file_path=full_path, to_upper=True)

#             for tupe in list_of_tuples:
#                 container = tupe[0]
#                 strand = tupe[1]
#                 start_pos = tupe[2]
#                 end_pos = tupe[3]
#                 ortho_to = tupe[4]
#                 misc = tupe[5]

#                 sequences.append(seq)
#                 paths.append(df_file_path)
#                 scores.append(score)
#                 strands.append(strand)
#                 start_positions.append(start_pos)
#                 end_positions.append(end_pos)
#                 chromosomes_or_genes.append(container)
#                 species_list.append(species)
#                 ortho_list.append(ortho_to)
#                 misc_list.append(misc)
            
    
#     pandas.DataFrame(list(zip(
#         sequences,
#         species_list,
#         scores,
#         chromosomes_or_genes,
#         strands,
#         start_positions,
#         end_positions,
#         ortho_list,
#         misc_list,
#         paths,
#         )),
#     columns=['sequence', 'target', 'score', 'chromosome_or_gene',
#     'strand', 'start_position', 'end_position', 'orthologous_to', 'misc', 'path']).to_csv(output_csv_path, index=False)
    
#     print(f'{bcolors.BLUE}>{bcolors.RESET} Done. Check {output_csv_path} for the output.')
#     return output_csv_path


# ARTIFACT CODE BELOW

# def graph_size_dist(
#     beta: int,
#     exp_name: str, 
#     output_dir: str,
#     size_of_solutions_for_n_trials: list[int],
#     ) -> None:

#     print('Drawing size distribution graph...')
#     output_path = os.path.join(output_dir, exp_name + '_b' + str(beta) + '_size_hist.png')

#     y = numpy.asarray(size_of_solutions_for_n_trials)
#     average_size_over_all_trials = numpy.mean(y)

#     n = len(size_of_solutions_for_n_trials)
#     x = numpy.arange(n)

#     title_1 = 'Set size distribution for {n} trials'.format(n=len(size_of_solutions_for_n_trials))
#     title_2 = 'Average size over all trials: {n:.2f}'.format(n=average_size_over_all_trials)
#     title_3 = 'Beta: {b}'.format(b=beta)

#     mininum_size = min(size_of_solutions_for_n_trials)
#     red_label = 'Smallest size: {n}'.format(n=mininum_size)
#     blue_mask = y != mininum_size
#     red_mask = y == mininum_size

#     matplotlib.pyplot.figure(figsize=(20, 3))

#     matplotlib.pyplot.bar(x=x[blue_mask], height=y[blue_mask])
#     matplotlib.pyplot.bar(x=x[red_mask], height=y[red_mask], color='red', label=red_label)

#     matplotlib.pyplot.legend()
#     matplotlib.pyplot.title(title_1 + '\n' + title_2 + '\n' + title_3)
#     matplotlib.pyplot.xlabel('Trial')
#     matplotlib.pyplot.ylabel('Set size')
#     matplotlib.pyplot.margins()

#     matplotlib.pyplot.tight_layout()
#     matplotlib.pyplot.savefig(output_path)


# def graph_score_dist(
#     beta: int,
#     exp_name: str, 
#     output_dir: str,
#     average_scores_for_n_trials: list[float],
#     ) -> None:

#     print('Drawing score distribution graph...')
#     output_path = os.path.join(output_dir, exp_name + '_b' + str(beta) + '_avg_score_hist.png')

#     title = 'Solution guides\' average score distribution over {n} trials'.format(n=len(average_scores_for_n_trials))

#     matplotlib.pyplot.figure(figsize=(10, 3))

#     matplotlib.pyplot.hist(average_scores_for_n_trials, edgecolor='black', linewidth=1.2)
    
#     matplotlib.pyplot.title(title)
#     matplotlib.pyplot.xlabel('Average score')
#     matplotlib.pyplot.ylabel('Count')

#     matplotlib.pyplot.tight_layout()
#     matplotlib.pyplot.savefig(output_path)


# def output_csv(
#     beta: int,
#     experiment_name: str,
#     output_directory: str, 
#     solver,
#     ) -> str:

#     output_path = os.path.join(output_directory, experiment_name + '_metrics.csv')

#     print(f'Writing to file: {output_path}')

#     avg_num_while_iters_for_n_trials = 0
#     if len(solver.num_while_iters_for_each_trial) > 0:
#         avg_num_while_iters_for_n_trials = numpy.mean(solver.num_while_iters_for_each_trial).round(2)

#     df = pandas.DataFrame({
#         'budget': [beta],
#         'avg_score_over_all_trials': [solver.average_score_for_all_trials],
#         'avg_score_over_each_trial': [solver.average_score_for_each_trial],
#         'avg_num_while_iters_for_n_trials': [avg_num_while_iters_for_n_trials],
#         'num_while_iters_for_each_trial': [solver.num_while_iters_for_each_trial],
#         'fractional_glop_vals': [solver.fractional_glop_vals],
#         'solver_time': [solver.solver_time],
#         'exhaustive_time': [solver.exhaustive_time],
#         'random_rounding_time': [solver.randomized_rounding_time],
#         'num_non_zero_feasible': [solver.num_non_zero_feasible],
#         'num_non_zero_feasible_under_one': [solver.num_feasible_guides_with_prob_lt_one],
#         'num_exhausted_combos': [solver.num_exhausted_combos],
#         'solved_with_exhaustive': ['Yes' if solver.solved_with_exhaustive else 'No'],
#         'size_of_solutions_for_n_trials': [solver.set_size_for_each_trial],
#         'num_species_constraints': [len(solver.species)],
#         'num_guide_variables': [solver.k],
#         'num_rounding_trials': [solver.num_trials],
#     })

#     write_header_if_file_doesnt_exist = not os.path.exists(output_path)  # otherwise append w/ mode='a'

#     df.to_csv(output_path, mode='a', header=write_header_if_file_doesnt_exist, index=False)
#     print(f'Done. Check {output_path} for the output.')
#     return output_path