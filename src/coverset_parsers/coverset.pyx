# distutils: language = c++

# NATIVE IMPORTS
import os
import sys
import pandas

# CYTHON CUSTOM IMPORTS
import coverset
from libcpp.string cimport string
from libcpp.unordered_set cimport unordered_set

sys.path.append("..")

# PYTHON CUSTOM IMPORTS
from classes.guide import Guide
from classes.species import Species
from scorers.scorer_factory import ScorerFactory
from classes.guide_container_factory import GuideContainerFactory


# Declare the class with cdef
cdef extern from "allegro/coverset.h" namespace "coversets":
    cdef cppclass CoversetCPP:
        CoversetCPP(size_t num_species, size_t guide_length, size_t num_trials) except +
        
        void encode_and_save_dna(string &seq, unsigned char score, unsigned short species_id)
        unordered_set[string] ortools_solver()
        string get_str()  # DEBUG only


cdef class CoversetsCython:
    cdef dict __dict__ 
    cdef CoversetCPP *coverset

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

        scorer_factory = ScorerFactory()
        scorer = scorer_factory.make_scorer(
            scorer_name=scorer_name,
            scorer_settings=scorer_settings,
        )

        guide_container_factory = GuideContainerFactory()

        print('Reading species input file from {path}'.format(
            path=input_species_csv_file_path)
            )
        species_df = pandas.read_csv(input_species_csv_file_path)

        num_species = species_df.shape[0]
        self.coverset = new CoversetCPP(num_species, guide_length, num_trials)

        # Make the species objects
        for row in species_df.itertuples():
            idx = row.Index
            cds_path = os.path.join(input_cds_directory, row.cds_file_name)
            genome_path = os.path.join(input_genome_directory, row.genome_file_name)

            species_object = Species(
                cds_path=cds_path,
                guide_scorer=scorer,
                name=row.species_name,
                genome_path=genome_path,
                guide_source=guide_source,
                guide_container_factory=guide_container_factory,
            )

            guide_objects_list: list[Guide] = list()
            guides_attributes_list: list[dict] = list()


            if cas_variant == 'cas9':
                guide_objects_list = species_object.get_cas9_guides()
            else:
                print('No such cas variant as', cas_variant, '. Modify this value in config.yaml. Exiting.\n')
                raise NotImplementedError


            for guide_object in guide_objects_list:
                guide_sequence = guide_object.sequence

                # interact with C++ -- encode and pass the sequence string
                python_bytes = guide_sequence.encode('utf-8')
                self.coverset.encode_and_save_dna(python_bytes, guide_object.score, idx)
                
            print('Done with', idx + 1, 'species...')
        print('Created coversets for all species.')

        # Call the C++ solver
        winners = self.ortools_solver()


    def ortools_solver(self):
        return self.coverset.ortools_solver()


    def __dealloc__(self):
        del self.coverset
