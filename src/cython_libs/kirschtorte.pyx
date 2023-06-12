# PYTHON LIBS
import os
import re
import sys
import pandas

# ALLEGRO CYTHON CUSTOM LIBS
# C++ function setup_and_solve() returns type vector[tuple[string, int, string]]
import kirschtorte  # Cythonized custom C++ lib, AKA kirschtorte.so
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

# ---
sys.path.append('..')  # Required to import below

# ALLEGRO PYTHON CUSTOM LIBS
from classes.guide import Guide
from classes.species import Species
from scorers.scorer_factory import ScorerFactory
from classes.guide_container import GuideContainer
from classes.guide_container_factory import GuideContainerFactory
import utils.records_count_finder as records_count_finder

# ctypedef (string, int, string) tuple_string_int_string

# Declare the class with cdef
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
        # species_df .csv normally found in data/input/species.csv
        vector[pair[string, string]] setup_and_solve(
            size_t monophonic_threshold,
            size_t cut_multiplicity,
            size_t beta)

        # Encodes each seq into a boost::dynamic_bitset and saves it, plus the species 
        # this seq hits.
        int encode_and_save_dna(
            string seq,
            size_t score,
            size_t container_id)


cdef class KirschtorteCython:
    cdef dict __dict__  # Enable cython self.attribute binding.
    cdef Kirschtorte *kirschtorte  # Pointer to C++ class instance.

    def __cinit__(
        self,
        beta: int,
        cut_multiplicity: int,
        monophonic_threshold: int,
        track: str,
        num_trials: int,
        scorer_name: str,
        cas_variant: str,
        guide_length: int,
        scorer_settings: dict,
        output_directory: str,
        input_directory: str,
        file_column_name: str,
        input_species_csv_file_path: str,
        ) -> None:

        try:
            print('Reading species input file from {path}'.format(path=input_species_csv_file_path))
            self.species_df = pandas.read_csv(input_species_csv_file_path)
        except pandas.errors.EmptyDataError:
            print('EmptyDataError exception in kirschtorte.pyx: File', input_species_csv_file_path, 'is empty. Exiting.')
            sys.exit(1)
        except FileNotFoundError:
            print('FileNotFoundError exception in kirschtorte.pyx: Cannot find file', input_species_csv_file_path, 'Did you spell the path/file name correctly? Exiting.')
            sys.exit(1)

        self.beta = beta
        self.cas_variant = cas_variant
        self.cut_multiplicity = cut_multiplicity
        self.monophonic_threshold = monophonic_threshold
        self.input_directory = input_directory
        self.file_column_name = file_column_name

        clusters = self.species_df.shape[0] if track == 'track_a' else records_count_finder.count_records(self.input_directory)

        self.kirschtorte = new Kirschtorte(clusters, guide_length, num_trials, output_directory.encode('utf-8'))

        # To translate indices back to legible names later.
        self.guide_origin: dict[int, str] = dict()

        # Instantiate the appropriate scorer
        scorer_factory = ScorerFactory()
        self.scorer = scorer_factory.make_scorer(scorer_name=scorer_name, scorer_settings=scorer_settings)

        # Instantiate the GuideContainerFactory to assign to each Species.
        # Each Species object tells the factory about its guide source.
        # This makes mix-and-matching the guide sources easier. You'd have to
        # only add another column to the species_df .csv file and read it in in the loop below.
        self.guide_container_factory = GuideContainerFactory()

        if track == 'track_a':
            self.track_a()
        elif track == 'track_e':
            self.track_e()


    def track_a(self) -> list[tuple[str, list[str]]]:
        # Make an object for each species
        for idx, row in self.species_df.iterrows():
            self.guide_origin[idx] = row.species_name
            
            records_path = os.path.join(self.input_directory, row[self.file_column_name])
            
            species_object = Species(
                name=row.species_name,
                records_path=records_path,
                guide_scorer=self.scorer,
                guide_container_factory=self.guide_container_factory,
            )

            guide_objects_list: list[Guide] = list()

            if self.cas_variant == 'cas9':
                species_object.make_guide_containers()
                guide_objects_list: list[Guide] = species_object.get_guides_from_containers()

                if len(guide_objects_list) == 0:
                    print('* WARNING: Species', row.species_name, 'contains no cas9 guides, or ' +
                    'all of its cas9 guides have been marked as repetitive and thus removed in ' +
                    'a preprocessing step. Set the include_repetitive option to False in config.yaml ' +
                    'to include them. Excluding', row.species_name, 'from further consideration.')
            else:
                print('No such cas variant as', self.cas_variant, '. Modify this value in config.yaml. Exiting.\n')
                raise NotImplementedError

            for guide_object in guide_objects_list:
                guide_sequence = guide_object.sequence
                
                # interact with C++ -- encode and pass the sequence string, score, and index
                self.kirschtorte.encode_and_save_dna(guide_sequence.encode('utf-8'), guide_object.score, idx)
                
            print('Done with', idx + 1, 'species...', end='\r')
        print('\nCreated coversets for all species.')
        print('Setting up the linear program...')

        # Deallocate
        del self.species_df

        # Interface with the C++ functions.
        winners_bytes_pairs = self.kirschtorte.setup_and_solve(self.monophonic_threshold, self.cut_multiplicity, self.beta)

        self.solution: list[tuple[str, list[str]]] = list()
        for seq, binary_hits in winners_bytes_pairs:
            # Decode the bytes object and reverse the binary string
            binary_hits = binary_hits.decode('utf-8')[::-1]  # e.g., '0111'
            seq = seq.decode('utf-8')  # e.g., 'ACCTGAG...'
            
            # e.g., ['saccharomyces', 'yarrowia', 'kmarx', ...]
            names_hits: list[str] = list()
            
            # Find all the indices where you have a '1' and transform back to actual species names
            for idx in [idx.start() for idx in re.finditer('1', binary_hits)]:
                names_hits.append(self.guide_origin[idx])
            
            # E.g., ('ACCTGAG...', ['saccharomyces', 'yarrowia', 'kmarx'])
            self.solution.append((seq, names_hits))
    

    def track_e(self) -> list[tuple[str, list[str]]]:
        container_idx: int = 0
        # Make an object for each species
        for _, row in self.species_df.iterrows():
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
                    print('* WARNING: Species', row.species_name, 'contains no cas9 guides, or ' +
                    'all of its cas9 guides have been marked as repetitive and thus removed in ' +
                    'a preprocessing step. Set the include_repetitive option to False in config.yaml ' +
                    'to include them. Excluding', row.species_name, 'from further consideration.')
            else:
                print('No such cas variant as', self.cas_variant, '. Modify this value in config.yaml. Exiting.\n')
                raise NotImplementedError

            for guide_container in guide_containers_list:
                guide_attributes = guide_container.get_attributes_dict()

                record_ortho_to = guide_attributes['record_orthologous_to']
                record_string_id = guide_attributes['record_string_id']

                container_target_name = record_ortho_to if record_ortho_to != 'N/A' else record_string_id

                guide_objects_list: list[Guide] = guide_container.get_cas9_guides()

                for guide_object in guide_objects_list:
                    guide_sequence = guide_object.sequence
                    
                    # interact with C++ -- encode and pass the sequence string, score, and index
                    status = self.kirschtorte.encode_and_save_dna(
                        guide_sequence.encode('utf-8'),
                        guide_object.score,
                        container_idx,
                        )

                    if status == 0:
                        self.guide_origin[container_idx] = row.species_name + ', ' + container_target_name
                
                container_idx += 1
                
                print('Done with', container_idx, 'genes...', end='\r')
        print('\nCreated coversets for all genes.')
        print('Setting up the linear program...')

        # Deallocate
        del self.species_df

        # Interface with the C++ functions.
        winners_bytes_pairs = self.kirschtorte.setup_and_solve(self.monophonic_threshold, self.cut_multiplicity, self.beta)

        self.solution: list[tuple[str, list[str]]] = list()
        for seq, binary_hits in winners_bytes_pairs:
            # Decode the bytes object and reverse the binary string
            binary_hits = binary_hits.decode('utf-8')[::-1]  # e.g., '0111'
            seq = seq.decode('utf-8')  # e.g., 'ACCTGAG...'
            
            # e.g., ['saccharomyces, LYS2', 'yarrowia, URA3', 'kmarx, LYS3', ...]
            names_hits: list[str] = list()
            
            # Find all the indices where you have a '1' and transform back to actual species names
            for idx in [idx.start() for idx in re.finditer('1', binary_hits)]:
                names_hits.append(self.guide_origin[idx])
            
            # E.g., ('ACCTGAG...', ['saccharomyces, LYS2', 'yarrowia, URA3', 'kmarx, LYS3'])
            self.solution.append((seq, names_hits))


    @property
    def species_names(self) -> list[str]:
        return list(self.guide_origin.values())


    def __dealloc__(self):
        del self.kirschtorte
