# distutils: language = c++

# NATIVE IMPORTS
import os
import re
import gc
import sys
import pandas

# ALLEGRO CYTHON CUSTOM IMPORTS
# C++ function ortools_solver() returns type vector[pair[string, string]]
import coverset
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

sys.path.append('..')  # Required to import below

# ALLEGRO PYTHON CUSTOM IMPORTS
from classes.guide import Guide
from classes.species import Species
from scorers.scorer_factory import ScorerFactory
from classes.guide_container_factory import GuideContainerFactory


# Declare the class with cdef
cdef extern from "allegro/coverset.h" namespace "coversets":
    cdef cppclass CoversetCPP:
        CoversetCPP(size_t num_species, size_t guide_length, size_t num_trials) except +
        
        # Google OR-Tools solver code
        # Returns a python-iterable list of tuples. tuple[0] is a plain-text sequence
        # such as 'ACTG...', and tuple[1] is the string representation of the bitset
        # that shows which species this sequence hits, e.g., '1110'.
        # The width of the bitset is equal to the number of input species found in
        # species_df .csv normally found in data/input/species.csv
        vector[pair[string, string]] ortools_solver()

        # Encodes each seq into a boost::dynamic_bitset and saves it, plus the species 
        # this seq hits.
        void encode_and_save_dna(string &seq, unsigned char score, unsigned short species_id)


cdef class CoversetsCython:
    cdef dict __dict__  # Enable cython self.attribute binding.
    cdef CoversetCPP *coverset  # Pointer to C++ class instance

    def __cinit__(
        self,
        num_trials: int,
        scorer_name: str,
        cas_variant: str,
        guide_length: int,
        guide_source: str,
        scorer_settings: dict,
        input_cds_directory: str,
        input_genome_directory: str,
        input_species_csv_file_path: str,
        ) -> None:

        try:
            print('Reading species input file from {path}'.format(path=input_species_csv_file_path))
            species_df = pandas.read_csv(input_species_csv_file_path)
        except pandas.errors.EmptyDataError:
            print('Error in coverset.pyx: File', input_species_csv_file_path, 'is empty. Exiting.')
            sys.exit(1)

        self.num_species = species_df.shape[0]

        self.coverset = new CoversetCPP(self.num_species, guide_length, num_trials)

        # To translate indices back to legible names later.
        self.idx_to_species: dict[int, str] = dict()

        # Instantiate the proper scorer
        scorer_factory = ScorerFactory()
        scorer = scorer_factory.make_scorer(scorer_name=scorer_name, scorer_settings=scorer_settings)

        # Instantiate the GuideContainerFactory to assign to each Species.
        # Each Species object tells the factory about its guide source.
        # This makes mix-and-matching the guide sources easier. You'd have to
        # only add another column to species_df .csv file and read it in in the loop below.
        guide_container_factory = GuideContainerFactory()

        # Make an object for each species
        for row in species_df.itertuples():
            idx = row.Index
            self.idx_to_species[idx] = row.species_name

            cds_path = os.path.join(input_cds_directory, row.cds_file_name)
            genome_path = os.path.join(input_genome_directory, row.genome_file_name)
            
            species_object = Species(
                name=row.species_name,
                cds_path=cds_path,
                genome_path=genome_path,
                guide_scorer=scorer,
                guide_source=guide_source,
                guide_container_factory=guide_container_factory,
            )

            guide_objects_list: list[Guide] = list()

            if cas_variant == 'cas9':
                guide_objects_list = species_object.get_cas9_guides()

                if len(guide_objects_list) == 0:
                    print('* WARNING: Species', row.species_name, 'contains no cas9 guides, or ' +
                    'all of its cas9 guides have been marked as repetitive and thus removed in ' +
                    'a preprocessing step. Set the include_repetitive option to False in config.yaml ' +
                    'to include them. Excluding', row.species_name, 'from further consideration.')
            else:
                print('No such cas variant as', cas_variant, '. Modify this value in config.yaml. Exiting.\n')
                raise NotImplementedError

            for guide_object in guide_objects_list:
                guide_sequence = guide_object.sequence
                
                # interact with C++ -- encode and pass the sequence string, score, and index
                self.coverset.encode_and_save_dna(guide_sequence.encode('utf-8'), guide_object.score, idx)
                
            print('Done with', idx + 1, 'species...')
        print('Created coversets for all species.')

        # Free up memory
        del species_df
        gc.collect()

        # TODO: REMOVE. DEBUGGING
        scorer.guide_finder.write_removed_guides_to_dataframe()

        # Interface with C++ functions.
        winners_bytes_pairs = self.coverset.ortools_solver()

        self.solution: list[tuple[str, list[str]]] = list()
        for seq, binary_hits in winners_bytes_pairs:
            # Decode the bytes object and reverse the binary string
            binary_hits = binary_hits.decode('utf-8')[::-1]  # e.g., '0111'
            seq = seq.decode('utf-8')  # e.g., 'ACCTGAG...'
            
            # e.g., ['saccharomyces', 'yarrowia', 'kmarx', ...]
            names_hits: list[str] = list()
            
            # Find all the indices where you have a '1' and transform back to actual species names
            for idx in [idx.start() for idx in re.finditer('1', binary_hits)]:
                names_hits.append(self.idx_to_species[idx])
            
            # E.g., ('ACCTGAG...', ['saccharomyces', 'yarrowia', 'kmarx'])
            self.solution.append((seq, names_hits))


    @property
    def species_names(self) -> list[str]:
        return list(self.idx_to_species.values())


    def __dealloc__(self):
        del self.coverset
