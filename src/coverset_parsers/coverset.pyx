# distutils: language = c++

# NATIVE IMPORTS
import os
import sys
import pandas

# CYTHON CUSTOM IMPORTS
from libcpp.string cimport string
import coverset

sys.path.append("..")

# PYTHON CUSTOM IMPORTS
from classes.guide import Guide
from classes.species import Species
from scorers.scorer_factory import ScorerFactory
from classes.guide_container_factory import GuideContainerFactory


# Declare the class with cdef
cdef extern from "allegro/coverset.h" namespace "coversets":
    cdef cppclass CoversetCPP:
        CoversetCPP(size_t num_species, size_t guide_length) except +
        
        void encode_and_save_dna(string &, unsigned char, unsigned short)
        void ortools_solver()
        string get_str()

# Create a Cython extension type which holds a C++ instance
#  as an attribute and create several forwarding methods
# This is the interface exposed to a Python script.
#cdef class PyCoverset:
#    cdef CoversetCPP *coverset  # Hold the C++ instance we're wrapping
#
#    def __cinit__(self):
#        self.coverset = new CoversetCPP()
#        
#    def __dealloc__(self):
#        del self.coverset
#
#    def OR_Tools_Solver(self):
#        self.coverset.ortools_solver()


cdef class CoversetsCython:
    cdef dict __dict__ 
    cdef CoversetCPP *coverset

    def __cinit__(
        self,
        scorer_name: str,
        cas_variant: str,
        guide_length: int,
        guide_source: str,
        scorer_settings: dict,
        input_cds_directory: str,
        input_genome_directory: str,
        input_species_csv_file_path: str,
        ) -> None:

        self.species_set = set()
        self.species_names = list()

        scorer_factory = ScorerFactory()
        scorer = scorer_factory.make_scorer(
            scorer_name=scorer_name,
            scorer_settings=scorer_settings,  # type: ignore
        )

        guide_container_factory = GuideContainerFactory()

        print('Reading species input file from {path}'.format(
            path=input_species_csv_file_path)
            )
        species_df = pandas.read_csv(input_species_csv_file_path)

        num_species = species_df.shape[0]
        self.coverset = new CoversetCPP(num_species, guide_length)

        # Make the species objects
        for row in species_df.itertuples():
            idx = row.Index
            self.species_set.add(idx)  # {0, 1, ..., num_species-1}
            self.species_names.append(row.species_name)  # ['kmarxianus', 'scerevisiae', ...]

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
                print('No such cas variant as', cas_variant, '.\n')
                raise NotImplementedError


            for guide_object in guide_objects_list:
                guide_sequence = guide_object.sequence

                print(guide_sequence)

                if 'N' in guide_sequence:
                    continue

                # interact with CPP -- pass in the sequence string
                python_bytes = guide_sequence.encode('utf-8')
                self.coverset.encode_and_save_dna(python_bytes, 1, idx)
                
            print('Done with', idx + 1, 'species...')

        print('Created coversets for all species.')

        # Call solver in CPP
        print(self.coverset.get_str())
        # self.coverset.ortools_solver()


    def __dealloc__(self):
        del self.coverset


    def ortools_solver(self):
        self.coverset.ortools_solver()