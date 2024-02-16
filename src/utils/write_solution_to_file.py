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


def output_solution_to_text(
    species_names: list[str],
    solution: list[tuple[str, float, list[str]]],
    experiment_name: str,
    output_directory: str,
    ) -> None:

    output_txt_path = os.path.join(output_directory, experiment_name + '.txt')
    with open(output_txt_path, 'w') as f:
        for tup in solution:
            seq = tup[0]
            score = tup[1]
            hit_species = tup[2]

            f.write(f'Guide {seq} targets {len(hit_species)} locations.\n')

        f.write(f'We can cut the following {len(species_names)} species: {str(species_names)}.\n')
        f.write(f'Using the following {str(len(solution))} guides: {str(solution)}.\n')


def output_solution_to_library(
    solution: list[tuple[str, float, list[str]]],
    experiment_name: str,
    output_directory: str,
    ) -> None:

    output_library_path = os.path.join(output_directory, experiment_name + '_library.txt')

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
    with open(output_library_path, 'w') as f:
        for g in library:
            f.write(f'{g}\n')
    print(f'{bcolors.BLUE}>{bcolors.RESET} Your library is ready. Find it in {output_library_path}.')


def cross_reference_solution_to_input(
    solution: list[tuple[str, float, list[str]]],
    experiment_name: str,
    input_csv_path: str,
    output_directory: str
    ) -> None:

    output_csv_path = os.path.join(output_directory, experiment_name + '.csv')

    input_df = pandas.read_csv(input_csv_path)

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
    
    input_df[input_df['sequence'].isin(library)].to_csv(output_csv_path, index=False)

    print(f'{bcolors.BLUE}>{bcolors.RESET} Done. Check {output_csv_path} for the output.')


def align_solution_to_input_bowtie(
    pam: str,
    solution: list[tuple[str, float, list[str]]],
    experiment_name: str,
    input_csv_path: str,
    input_directory: str,
    paths_csv_column_name: str,
    species_names_csv_column_name: str,
    output_directory: str,
    ) -> tuple[str, float]:

    output_csv_path = os.path.join(output_directory, experiment_name + '.csv')

    output_bowtie_temp_dir = os.path.join(output_directory, 'temp')
    output_align_temp_path = os.path.join(output_directory, 'temp', 'temp_library.txt')

    if not os.path.exists(output_bowtie_temp_dir):
        os.mkdir(output_bowtie_temp_dir)

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

    total_time_elapsed = 0.0
    print(f'{bcolors.BLUE}>{bcolors.RESET} Aligning the library against the input species to generate the final report...')

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
