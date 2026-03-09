import os
import re
import sys
import pandas

from allegro._kirschtorte_core import _KirschtorteEngine

from allegro.classes.species import Species
from allegro.utils.shell_colors import bcolors
from allegro.utils import records_count_finder
from allegro.classes.guide_container_factory import GuideContainerFactory

class KirschtorteCython:
    def __init__(
        self, 
        beta, 
        track,
        multiplicity,
        monophonic_threshold,
        preclustering,
        seed_length,
        mismatched_allowed_after_seed,
        early_stopping_patience,
        input_directory,
        output_directory, 
        input_species_csv_file_path,
        input_species_path_column,
        enable_solver_diagnostics):

        self.multiplicity = multiplicity
        self.input_directory = input_directory
        self.input_species_path_column = input_species_path_column
        self.species_df = pandas.read_csv(input_species_csv_file_path)
        self.clusters = self.species_df.shape[0] if track == 'track_a' else \
            records_count_finder.count_records(input_directory, self.species_df[input_species_path_column].tolist())
        
        # Instantiate the C++ engine
        self.engine = _KirschtorteEngine(
            preclustering, 
            beta, 
            seed_length,
            self.multiplicity,
            self.clusters, 
            monophonic_threshold, 
            early_stopping_patience,
            mismatched_allowed_after_seed, 
            enable_solver_diagnostics, 
            output_directory
        )

        # To translate indices back to legible names later.
        self.guide_origin: dict[int, str] = dict()

        # Instantiate the GuideContainerFactory to assign to each Species.
        # Each Species object tells the factory about its guide source.
        # This makes mix-and-matching the guide sources easier. You'd have to
        # only add another column to the species_df .csv file and read it in in the loop below.
        self.guide_container_factory = GuideContainerFactory()

        # Choose the appropriate track.
        if track == 'track_a':
            self.track_a()
        elif track == 'track_e':
            self.track_e()

    @property
    def species_names(self) -> list[str]:
        return list(self.guide_origin.values())

    def track_a(self) -> list[tuple[str, float, list[str]]]:
        total_number_of_guides: int = 0

        # Make an object for each species
        for idx, row in self.species_df.iterrows():
            self.guide_origin[idx] = row.species_name
            total_available_guides_for_this_species = 0

            records_path = os.path.join(self.input_directory, row[self.input_species_path_column])
            
            species_object = Species(
                name=row.species_name,
                records_path=records_path,
                guide_container_factory=self.guide_container_factory)

            guide_objects_list: list[Guide] = list()

            species_object.make_guide_containers()
            guide_objects_list: list[Guide] = species_object.get_guides_from_containers()

            if len(guide_objects_list) == 0:
                print(f'{bcolors.RED}> Warning{bcolors.RESET}: Species {row.species_name} contains no target sequences, or ' +
                    f'all of its targets have been marked as repetitive and thus removed in ' +
                    f'a preprocessing step. First, check the input file. Set the filter_repetitive and/or filter_by_gc option(s) to False in config.yaml ' +
                    f'to include them, or remove this species from your input file and try again.')
                continue

            total_available_guides_for_this_species += len(guide_objects_list)
            total_number_of_guides += len(guide_objects_list)

            for guide_object in guide_objects_list:
                guide_sequence = guide_object.sequence
                
                # interact with C++ -- encode and pass the sequence string, score, and index.
                self.engine.encode_and_save_dna(guide_sequence.encode('utf-8'), guide_object.score, idx)
                
            if total_available_guides_for_this_species < self.multiplicity:
                print(f'{bcolors.RED}> Error{bcolors.RESET}: The genes in species {row.species_name} ' +
                    f'contain fewer total targets ({total_available_guides_for_this_species}) than ' +
                    f'the requested multiplicity {self.multiplicity} for Track A. ' +
                    f'Either remove this species from your input file, or reduce your multiplicity ' +
                    f'to at most the total available targets for this species ({total_available_guides_for_this_species}) and try again. Exiting.')
                sys.exit(1)

            print(f'{bcolors.BLUE}>{bcolors.RESET} Done with {idx + 1}/{self.clusters} species...', end='\r')
        print(f'\n{bcolors.BLUE}>{bcolors.RESET} Created coversets for all species containing a total of {total_number_of_guides} targets.')
        print(f'{bcolors.BLUE}>{bcolors.RESET} Setting up and solving the linear program...')

        # Interface with the C++ functions.
        return self.setup_and_solve()
    
    def track_e(self) -> list[tuple[str, float, list[str]]]:
        container_idx: int = 0
        total_number_of_guides: int = 0
        
        for _, row in self.species_df.iterrows():  # Make an object for each species.
            total_available_guides_for_this_species = 0

            records_path = os.path.join(self.input_directory, row[self.input_species_path_column])
            
            species_object = Species(
                name=row.species_name,
                records_path=records_path,
                guide_container_factory=self.guide_container_factory)

            guide_containers_list: list[GuideContainer] = list()

            species_object.make_guide_containers()
            guide_containers_list = species_object.guide_containers_list

            if len(guide_containers_list) == 0:
                print(f'{bcolors.RED}> Error{bcolors.RESET}: Species {row.species_name} contains no target sequnces, or ' +
                    f'all of its targets have been marked as repetitive and thus removed in ' +
                    f'a preprocessing step. First, check the input file. Set the filter_repetitive and/or filter_by_gc option(s) to False in config.yaml ' +
                    f'to include them, or remove this species from your input file and try again. Exiting.')
                sys.exit(1)

            for guide_container in guide_containers_list:
                guide_objects_list: list[Guide] = guide_container.get_guides()

                if len(guide_objects_list) < self.multiplicity:
                    print(f'{bcolors.RED}> Warning{bcolors.RESET}: In species {row.species_name}, gene {guide_container.string_id} ' +
                        f'contains fewer targets ({len(guide_objects_list)}) than the requested multiplicity ({self.multiplicity}) for Track E. Discarding this gene.')
                    continue
                
                total_number_of_guides += len(guide_objects_list)
                total_available_guides_for_this_species += len(guide_objects_list)

                for guide_object in guide_objects_list:
                    guide_sequence = guide_object.sequence
                    
                    # Interact with C++ -- encode and pass the sequence string, score, and index.
                    status = self.engine.encode_and_save_dna(
                        guide_sequence.encode('utf-8'),
                        guide_object.score,
                        container_idx)

                    if status == 0:
                        self.guide_origin[container_idx] = row.species_name # + ', ' + container_target_name
                
                container_idx += 1
                
                print(f'{bcolors.BLUE}>{bcolors.RESET} Done with {container_idx}/{self.clusters} genes...', end='\r')

            if total_available_guides_for_this_species < self.multiplicity:
                print(f'{bcolors.RED}> Error{bcolors.RESET}: The genes in species {row.species_name} ' +
                f'contain fewer total target sequences ({total_available_guides_for_this_species}) than ' +
                f'the requested multiplicity {self.multiplicity} for Track E. ' +
                f'Either remove this species from your input file, or reduce your multiplicity ' +
                f'to at most the total available targets for this species ({total_available_guides_for_this_species}) and try again. Exiting.')
                sys.exit(1)
        print()
        print(f'{bcolors.BLUE}>{bcolors.RESET} Created coversets for all species containing a total of {total_number_of_guides} targets.')
        print(f'{bcolors.BLUE}>{bcolors.RESET} Setting up and solving the linear program...')

        # Interface with the C++ functions.
        return self.setup_and_solve()
        
    def setup_and_solve(self) -> list[tuple[str, float, list[str]]]:
        # Call C++ engine
        guide_struct_vector = self.engine.setup_and_solve()

        # Nichts zu tun
        if len(guide_struct_vector) == 0:
            sys.exit(1)

        self.solution: list[tuple[str, float, list[str]]] = list()
        for guide_struct in guide_struct_vector:
            seq = guide_struct.sequence
            score = guide_struct.score
            binary_hits = guide_struct.species_hit

            # Decode the bytes object and reverse the binary string.
            binary_hits = binary_hits[::-1] #.decode('utf-8')[::-1]  # e.g., '0111'
            # seq = seq#.decode('utf-8')  # e.g., 'ACCTGAG...'
            
            # e.g., ['saccharomyces_cerevisiae', 'yarrowia_lipolytica', 'kluyveromyces_marxianus', ...]
            names_hits: list[str] = list()
            
            # Find all the indices where you have a '1' and transform back to actual species names.
            for idx in [idx.start() for idx in re.finditer('1', binary_hits)]:
                names_hits.append(self.guide_origin[idx])
                # -- need to know exactly where the guide came from. not just its species
            
            # E.g., ('ACCTGAG...', 6, ['saccharomyces', 'yarrowia', 'kmarx']).
            self.solution.append((seq, score, names_hits))
            
        return self.solution
