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


def count_kmers(sequence, k):
    kmers = dict()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]

        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1

    return kmers


# Declare the class with cdef
cdef extern from "allegro/coverset.h" namespace "coversets":
    cdef cppclass CoversetCPP:
        CoversetCPP(size_t num_species, size_t guide_length, size_t num_trials) except +
        
        void encode_and_save_dna(string &, unsigned char, unsigned short)
        void ortools_solver()
        string get_str()


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

        self.species_set = set()
        self.species_names = list()

        illegal_characters = ['N', 'W', 'R', 'Y', 'K']

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
        self.coverset = new CoversetCPP(num_species, guide_length, num_trials)

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


                # Skip guides with bad nucleotides.
                if any(c in guide_sequence for c in illegal_characters):
                    continue

                # Skip guides such as GGAGGAGGAGGAGGAGGAGG where GG is repeated
                # 7 times or GA is repeated 6 times.
                if max(count_kmers(guide_sequence, 2).values()) > 4:
                    continue

                # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0784-0#Abs1:~:text=Repetitive%20bases%20are%20defined%20as%20any%20of%20the%20following
                if any(substr in guide_sequence for substr in ['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT']):
                    continue
                
                # https://www.nature.com/articles/s41467-019-12281-8#Abs1:~:text=but%20not%20significant).-,The%20contribution%20of,-repetitive%20nucleotides%20to
                if any(substr in guide_sequence for substr in ['AAAA', 'CCCC', 'GGGG', 'TTTT']):
                    continue

                # interact with CPP -- pass in the sequence string
                python_bytes = guide_sequence.encode('utf-8')
                self.coverset.encode_and_save_dna(python_bytes, 1, idx)
                
            print('Done with', idx + 1, 'species...')

        print('Created coversets for all species.')

        # Call solver in CPP
        # print(self.coverset.get_str())
        self.coverset.ortools_solver()


    def __dealloc__(self):
        del self.coverset


    def ortools_solver(self):
        self.coverset.ortools_solver()