# PYTHON LIBS.
import os
import re
import sys
import pandas
from queue import Queue
from threading import Thread, Semaphore, Lock

# ALLEGRO CYTHON CUSTOM LIBS.
# Cythonized custom C++ lib, AKA kirschtorte.so.
import kirschtorte  # type: ignore
from libcpp.string cimport string
from libcpp.vector cimport vector

# ---
sys.path.append('..')  # Required to import below.

# ALLEGRO PYTHON CUSTOM LIBS.
from classes.guide import Guide
from classes.species import Species
from utils.shell_colors import bcolors
from scorers.scorer_factory import ScorerFactory
from utils.offtarget_finder import OfftargetFinder
from classes.guide_container import GuideContainer
from classes.guide_container_factory import GuideContainerFactory
import utils.records_count_finder as records_count_finder


# Declare a C++ class with cdef.
cdef extern from "allegro/guide_struct.h":
    cdef cppclass GuideStruct:
        string sequence
        double score
        string species_hit


cdef extern from "allegro/kirschtorte.h" namespace "Kirschtorte":
    cdef cppclass Kirschtorte:
        Kirschtorte(
            size_t num_containers,
            size_t guide_length,
            size_t num_trials,
            string output_directory) except +
        
        # Google OR-Tools solver code
        # Returns a python-iterable list of tuples. tuple[0] is a plain-text sequence
        # such as 'ACTG...', and tuple[1] is the string representation of the bitset
        # that shows which species this sequence hits, e.g., '1110'.
        # The width of the bitset is equal to the number of input species found in
        # species_df .csv normally found in data/input/species.csv.
        vector[GuideStruct] setup_and_solve(
            size_t monophonic_threshold,
            size_t cut_multiplicity,
            size_t beta)

        # Encodes each seq into a boost::dynamic_bitset and saves it,
        # plus the species this seq hits.
        int encode_and_save_dna(
            string seq,
            double score,
            size_t container_id)


cdef class KirschtorteCython:
    cdef dict __dict__  # Enable cython self.attribute binding.
    cdef Kirschtorte *kirschtorte  # Pointer to C++ class instance.

    def __cinit__(
        self,
        beta: int,
        track: str,
        num_trials: int,
        max_threads: int,
        scorer_name: str,
        cas_variant: str,
        guide_length: int,
        input_directory: str,
        scorer_settings: dict,
        output_directory: str,
        file_column_name: str,
        cut_multiplicity: int,
        monophonic_threshold: int,
        # output_offtargets: bool,
        input_species_csv_file_path: str,
        # input_species_offtarget_dir: str,
        # input_species_offtarget_column: str,
        # report_up_to_n_mismatches: int,
        # seed_region_is_n_from_pam: int,
        ) -> None:

        self.beta = beta
        self.cas_variant = cas_variant
        self.input_directory = input_directory
        self.file_column_name = file_column_name
        self.cut_multiplicity = cut_multiplicity
        self.monophonic_threshold = monophonic_threshold
        # self.input_species_offtarget_dir = input_species_offtarget_dir
        # self.input_species_offtarget_column = input_species_offtarget_column
        # self.report_up_to_n_mismatches = report_up_to_n_mismatches
        # self.seed_region_is_n_from_pam = seed_region_is_n_from_pam
        # self.output_offtargets = output_offtargets
        # self.offtarget_finder = OfftargetFinder()
        self.species_df = pandas.read_csv(input_species_csv_file_path)

        self.lock = Lock()
        self.semaphore = Semaphore(max_threads)

        # Set how many clusters we need.
        # If track_a is selected, we have the same number of clusters as species. This is the same number for "Number of constraints" in solver_log.txt 
        # If track_e is selected, count the total number of fasta records in all species input files (uses multithreading). 
        #   This is the same number for "Number of constraints" in solver_log.txt aka total number of genes/chromosomes.
        clusters = self.species_df.shape[0] if track == 'track_a' else records_count_finder.count_records(self.input_directory)

        # Instantiate a C++ class. The linear programming part of ALLEGRO is done there.
        self.kirschtorte = new Kirschtorte(clusters, guide_length, num_trials, output_directory.encode('utf-8'))

        # To translate indices back to legible names later.
        self.guide_origin: dict[int, str] = dict()

        # Instantiate the appropriate scorer.
        scorer_factory = ScorerFactory()
        self.scorer = scorer_factory.make_scorer(scorer_name=scorer_name, scorer_settings=scorer_settings)

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

    
    # def create_and_init_species(self, row: pandas.Series, created_species_objects: list) -> None:
    #     records_path = os.path.join(self.input_directory, row[self.file_column_name])

    #     offtarget_background_path = ''
    #     if self.output_offtargets:
    #         offtarget_background_path = os.path.join(self.input_species_offtarget_dir, row[self.input_species_offtarget_column]) 
        
    #     with self.semaphore:
    #         species_object = Species(
    #             name=row.species_name,
    #             records_path=records_path,
    #             guide_scorer=self.scorer,
    #             guide_container_factory=self.guide_container_factory,
    #         )

    #         if self.output_offtargets:
    #             self.offtarget_finder.run_bowtie_build(species_name=species_object.name, path_to_background_fasta=offtarget_background_path)

    #         with self.lock:
    #             created_species_objects.append(species_object)


    # def offtarget_finder_run_bowtie_against_other(self, this_species_name: str, that_species_name: str, guide_seq_list: list[str], guide_objects_list: list[Guide]) -> None:
    #     with self.semaphore:
    #         genome_offtargets_indices = self.offtarget_finder.run_bowtie_against_other(this_species_name, that_species_name, guide_seq_list, self.report_up_to_n_mismatches, self.seed_region_is_n_from_pam)

    #         for idx, val in enumerate(genome_offtargets_indices):
    #             if val > 0:
    #                 with self.lock:
    #                     guide_objects_list[idx].has_offtargets = True


    # def track_a(self) -> list[tuple[str, float, list[str]]]:
    #     species_list = list[Species]
    #     threads: list[Thread] = list()

    #     # Limit to n concurrent threads
    #     created_species_objects: list[Species] = list()

    #     for idx, row in self.species_df.iterrows():
    #         self.guide_origin[idx] = row.species_name

    #         thread = Thread(target=self.create_and_init_species, args=(row, created_species_objects))
    #         thread.start()
    #         threads.append(thread)

    #     # Ensure all threads are completed
    #     for thread in threads:
    #         thread.join()

    #     threads = list()

    #     # Process objects from the list
    #     for species_object in created_species_objects:

    #         species_object.make_guide_containers()
    #         guide_objects_list: list[Guide] = species_object.get_guides_from_containers()
            
    #         if len(guide_objects_list) == 0:
    #             print(f'{bcolors.RED}> Warning{bcolors.RESET}: Species {species_object.name} contains no Cas9 guides, or ' +
    #             'all of its Cas9 guides have been marked as repetitive/blacklisted and thus removed in ' +
    #             'a preprocessing step. Set the filter_repetitive option to False in config.yaml ' +
    #             'to include them. Excluding', species_object.name, 'from further consideration.')
    #             continue

    #         if self.output_offtargets:
    #             guide_seq_list = self.offtarget_finder.create_and_return_reads(species_object.name, guide_objects_list)

    #             for other_species_object in created_species_objects:
    #                 # Outputs num_mismatches .sam files for each species
    #                 # Need to check each sam file for offtargets in non-seed region -- let's say, if we have at least 1
    #                 # exact offtarget in any, more than 20% of the species, discard the guide -- start with anywhere 
    #                 # -- check lib size
    #                 thread = Thread(target=self.offtarget_finder_run_bowtie_against_other, args=(species_object.name, other_species_object.name, guide_seq_list, guide_objects_list))
    #                 thread.start()
    #                 threads.append(thread)

    #             # Ensure all threads are completed
    #             for thread in threads:
    #                 thread.join()

    #             # for each guide obj, check if it has offtargets in any? (threshold? 30% of all?)

    #         for guide_object in guide_objects_list:
    #             if guide_object.has_offtargets:
    #                 continue
                
    #             # interact with C++ -- encode and pass the sequence string, score, and index.
    #             self.kirschtorte.encode_and_save_dna(guide_object.sequence.encode('utf-8'), guide_object.score, idx)
                
    #         print(f'{bcolors.BLUE}>{bcolors.RESET} Done with {idx + 1} species...', end='\r')
        
    #     print(f'\n{bcolors.BLUE}>{bcolors.RESET} Created coversets for all species.')
    #     print(f'{bcolors.BLUE}>{bcolors.RESET} Setting up and solving the linear program...')

    #     # Deallocate.
    #     del self.species_df

    #     # Interface with the C++ functions.
    #     guide_struct_vector = self.kirschtorte.setup_and_solve(self.monophonic_threshold, self.cut_multiplicity, self.beta)

    #     if guide_struct_vector.size() == 0:
    #         sys.exit(0)

    #     self.solution: list[tuple[str, float, list[str]]] = list()
    #     for guide_struct in guide_struct_vector:
    #         seq = guide_struct.sequence
    #         score = guide_struct.score
    #         binary_hits = guide_struct.species_hit

    #         # Decode the bytes object and reverse the binary string.
    #         binary_hits = binary_hits.decode('utf-8')[::-1]  # e.g., '0111'.
    #         seq = seq.decode('utf-8')  # e.g., 'ACCTGAG...'
            
    #         # e.g., ['saccharomyces', 'yarrowia', 'kmarx', ...].
    #         names_hits: list[str] = list()
            
    #         # Find all the indices where you have a '1' and transform back to actual species names.
    #         for idx in [idx.start() for idx in re.finditer('1', binary_hits)]:
    #             names_hits.append(self.guide_origin[idx])
            
    #         # E.g., ('ACCTGAG...', 6, ['saccharomyces', 'yarrowia', 'kmarx']).
    #         self.solution.append((seq, score, names_hits))

    #     return self.solution

    
    # def track_e(self) -> list[tuple[str, float, list[str]]]:
    #     container_idx: int = 0
        
    #     species_list = list[Species]
    #     threads: list[Thread] = list()

    #     # Limit to n concurrent threads
    #     created_species_objects: list[Species] = list()

    #     for idx, row in self.species_df.iterrows():
    #         self.guide_origin[idx] = row.species_name

    #         thread = Thread(target=self.create_and_init_species, args=(row, created_species_objects))
    #         thread.start()
    #         threads.append(thread)

    #     # Ensure all threads are completed
    #     for thread in threads:
    #         thread.join()

    #     threads = list()

    #     # Process objects from the list
    #     for species_object in created_species_objects:
    #         species_object.make_guide_containers()
    #         guide_containers_list = species_object.guide_containers_list

    #         # TODO fix this
    #         # guide_objects_list: list[Guide] = species_object.get_guides_from_containers()

    #         # if self.output_offtargets:
    #         #     guide_seq_list = self.offtarget_finder.create_and_return_reads(species_object.name, guide_objects_list)

    #         #     for other_species_object in created_species_objects:
    #         #         thread = Thread(target=self.offtarget_finder_run_bowtie_against_other, args=(species_object.name, other_species_object.name, guide_seq_list, guide_objects_list))
    #         #         thread.start()
    #         #         threads.append(thread)

    #         #     for thread in threads:
    #         #         thread.join()

    #         for guide_container in guide_containers_list:
    #             guide_attributes = guide_container.get_attributes_dict()

    #             record_ortho_to = guide_attributes['record_orthologous_to']
    #             record_string_id = guide_attributes['record_string_id']

    #             container_target_name = record_ortho_to if record_ortho_to != 'N/A' else record_string_id

    #             guide_objects_list: list[Guide] = guide_container.get_cas9_guides()
    #             for guide_object in guide_objects_list:
    #                 if guide_object.has_offtargets:
    #                     continue

    #                 # Interact with C++ -- encode and pass the sequence string, score, and index.
    #                 status = self.kirschtorte.encode_and_save_dna(
    #                     guide_sequence.encode('utf-8'),
    #                     guide_object.sequence.score,
    #                     container_idx,
    #                     )

    #                 if status == 0:
    #                     self.guide_origin[container_idx] = row.species_name # + ', ' + container_target_name
                
    #             container_idx += 1

    #             print(f'{bcolors.BLUE}>{bcolors.RESET} Done with {idx + 1} species...', end='\r')
    #     print(f'\n{bcolors.BLUE}>{bcolors.RESET} Created coversets for all species.')
    #     print(f'{bcolors.BLUE}>{bcolors.RESET} Setting up and solving the linear program...')

    #     #Deallocate.
    #     del self.species_df

    #     # Interface with the C++ functions.
    #     guide_struct_vector = self.kirschtorte.setup_and_solve(self.monophonic_threshold, self.cut_multiplicity, self.beta)

    #     if guide_struct_vector.size() == 0:
    #         sys.exit(0)

    #     self.solution: list[tuple[str, float, list[str]]] = list()
    #     for guide_struct in guide_struct_vector:
    #         seq = guide_struct.sequence
    #         score = guide_struct.score
    #         binary_hits = guide_struct.species_hit

    #         # Decode the bytes object and reverse the binary string.
    #         binary_hits = binary_hits.decode('utf-8')[::-1]  # e.g., '0111'
    #         seq = seq.decode('utf-8')  # e.g., 'ACCTGAG...'
            
    #         # e.g., ['saccharomyces, LYS2', 'yarrowia, URA3', 'kmarx, LYS3', ...]
    #         names_hits: list[str] = list()
            
    #         # Find all the indices where you have a '1' and transform back to actual species names.
    #         for idx in [idx.start() for idx in re.finditer('1', binary_hits)]:
    #             names_hits.append(self.guide_origin[idx])
            
    #         # E.g., ('ACCTGAG...', 6, ['saccharomyces, LYS2', 'yarrowia, URA3', 'kmarx, LYS3']).
    #         self.solution.append((seq, score, names_hits))
            
    #     return self.solution
  

    @property
    def species_names(self) -> list[str]:
        return list(self.guide_origin.values())


    def __dealloc__(self):
        del self.kirschtorte


    def track_a(self) -> list[tuple[str, float, list[str]]]:
        # Make an object for each species
        for idx, row in self.species_df.iterrows():
            self.guide_origin[idx] = row.species_name
            
            records_path = os.path.join(self.input_directory, row[self.file_column_name])
            
            species_object = Species(
                name=row.species_name,
                records_path=records_path,
                guide_scorer=self.scorer,
                guide_container_factory=self.guide_container_factory
            )

            guide_objects_list: list[Guide] = list()

            if self.cas_variant == 'cas9':
                species_object.make_guide_containers()
                guide_objects_list: list[Guide] = species_object.get_guides_from_containers()

                if len(guide_objects_list) == 0:
                    print(f'{bcolors.RED}> Warning{bcolors.RESET}: Species', row.species_name, 'contains no Cas9 guides, or ' +
                    'all of its Cas9 guides have been marked as repetitive and thus removed in ' +
                    'a preprocessing step. Set the filter_repetitive option to False in config.yaml ' +
                    'to include them. Excluding', row.species_name, 'from further consideration.')
            else:
                print(f'{bcolors.RED}> Warning{bcolors.RESET}: No such cas variant as {self.cas_variant}. Modify this value in config.yaml. Exiting.\n')
                sys.exit(1)

            for guide_object in guide_objects_list:
                guide_sequence = guide_object.sequence
                
                # interact with C++ -- encode and pass the sequence string, score, and index.
                self.kirschtorte.encode_and_save_dna(guide_sequence.encode('utf-8'), guide_object.score, idx)
                
            print(f'{bcolors.BLUE}>{bcolors.RESET} Done with {idx + 1} species...', end='\r')
        print(f'\n{bcolors.BLUE}>{bcolors.RESET} Created coversets for all species.')
        print(f'{bcolors.BLUE}>{bcolors.RESET} Setting up and solving the linear program...')

        # Deallocate.
        del self.species_df

        # Interface with the C++ functions.
        guide_struct_vector = self.kirschtorte.setup_and_solve(self.monophonic_threshold, self.cut_multiplicity, self.beta)

        if guide_struct_vector.size() == 0:
            sys.exit(0)

        self.solution: list[tuple[str, float, list[str]]] = list()
        for guide_struct in guide_struct_vector:
            seq = guide_struct.sequence
            score = guide_struct.score
            binary_hits = guide_struct.species_hit

            # Decode the bytes object and reverse the binary string.
            binary_hits = binary_hits.decode('utf-8')[::-1]  # e.g., '0111'.
            seq = seq.decode('utf-8')  # e.g., 'ACCTGAG...'
            
            # e.g., ['saccharomyces', 'yarrowia', 'kmarx', ...].
            names_hits: list[str] = list()
            
            # Find all the indices where you have a '1' and transform back to actual species names.
            for idx in [idx.start() for idx in re.finditer('1', binary_hits)]:
                names_hits.append(self.guide_origin[idx])
            
            # E.g., ('ACCTGAG...', 6, ['saccharomyces', 'yarrowia', 'kmarx']).
            self.solution.append((seq, score, names_hits))

        return self.solution
    

    def track_e(self) -> list[tuple[str, float, list[str]]]:
        container_idx: int = 0
        
        for _, row in self.species_df.iterrows():  # Make an object for each species.
            records_path = os.path.join(self.input_directory, row[self.file_column_name])
            
            species_object = Species(
                name=row.species_name,
                records_path=records_path,
                guide_scorer=self.scorer,
                guide_container_factory=self.guide_container_factory,
            )

            guide_containers_list: list[GuideContainer] = list()

            if self.cas_variant == 'cas9':
                species_object.make_guide_containers()
                guide_containers_list = species_object.guide_containers_list

                if len(guide_containers_list) == 0:
                    print(f'{bcolors.RED}> Warning{bcolors.RESET}: Species', row.species_name, 'contains no cas9 guides, or ' +
                    'all of its cas9 guides have been marked as repetitive and thus removed in ' +
                    'a preprocessing step. Set the filter_repetitive option to False in config.yaml ' +
                    'to include them. Excluding', row.species_name, 'from further consideration.')
            else:
                print(f'{bcolors.RED}> Warning{bcolors.RESET}: No such cas variant as {self.cas_variant}. Modify this value in config.yaml. Exiting.\n')
                sys.exit(1)

            for guide_container in guide_containers_list:
                # guide_attributes = guide_container.get_attributes_dict()

                # record_ortho_to = guide_attributes['record_orthologous_to']
                # record_string_id = guide_attributes['record_string_id']

                # container_target_name = record_ortho_to if record_ortho_to != 'N/A' else record_string_id

                guide_objects_list: list[Guide] = guide_container.get_cas9_guides()
                for guide_object in guide_objects_list:
                    guide_sequence = guide_object.sequence
                    
                    # Interact with C++ -- encode and pass the sequence string, score, and index.
                    status = self.kirschtorte.encode_and_save_dna(
                        guide_sequence.encode('utf-8'),
                        guide_object.score,
                        container_idx,
                        )

                    if status == 0:
                        self.guide_origin[container_idx] = row.species_name # + ', ' + container_target_name
                
                container_idx += 1
                
                print(f'{bcolors.BLUE}>{bcolors.RESET} Done with {container_idx} genes...', end='\r')
        print(f'\n{bcolors.BLUE}>{bcolors.RESET} Created coversets for all genes.')
        print(f'{bcolors.BLUE}>{bcolors.RESET} Setting up and solving the linear program...')

        # Deallocate.
        del self.species_df

        # Interface with the C++ functions.
        guide_struct_vector = self.kirschtorte.setup_and_solve(self.monophonic_threshold, self.cut_multiplicity, self.beta)

        if guide_struct_vector.size() == 0:
            sys.exit(0)

        self.solution: list[tuple[str, float, list[str]]] = list()
        for guide_struct in guide_struct_vector:
            seq = guide_struct.sequence
            score = guide_struct.score
            binary_hits = guide_struct.species_hit

            # Decode the bytes object and reverse the binary string.
            binary_hits = binary_hits.decode('utf-8')[::-1]  # e.g., '0111'
            seq = seq.decode('utf-8')  # e.g., 'ACCTGAG...'
            
            # e.g., ['saccharomyces, LYS2', 'yarrowia, URA3', 'kmarx, LYS3', ...]
            names_hits: list[str] = list()
            
            # Find all the indices where you have a '1' and transform back to actual species names.
            for idx in [idx.start() for idx in re.finditer('1', binary_hits)]:
                names_hits.append(self.guide_origin[idx])
            
            # E.g., ('ACCTGAG...', 6, ['saccharomyces, LYS2', 'yarrowia, URA3', 'kmarx, LYS3']).
            self.solution.append((seq, score, names_hits))
            
        return self.solution