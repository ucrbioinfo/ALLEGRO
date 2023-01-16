import os
import pandas

from classes.guide import Guide
from classes.species import Species
from scorers.scorer_factory import ScorerFactory
from classes.guide_container_factory import GuideContainerFactory


class Coversets:
    def __init__(
        self, 
        guide_source: str, 
        scorer_name: str, 
        scorer_settings: dict[str, str], 
        input_species_csv_file_path: str, 
        input_genome_directory: str, 
        input_cds_directory: str,
        ) -> None:

        self.species_set: set[int] = set()
        self.species_names: list[str] = list()
        self.cover_sets: dict[str, tuple[float, set[int]]] = dict()  # To be returned for the solver to use
        self.int_to_species_dict: dict[int, Species] = dict()
        self.seq_to_guides_dict: dict[str, list[Guide]] = dict()

        scorer_factory = ScorerFactory()
        scorer_obj = scorer_factory.make_scorer(
            scorer_name=scorer_name, 
            scorer_settings=scorer_settings,
            )

        guide_container_factory = GuideContainerFactory(guide_source=guide_source)

        print('Reading species input file from {path}'.format(path=input_species_csv_file_path))
        species_df = pandas.read_csv(input_species_csv_file_path)

        # Make the species objects
        # TODO: Add support for another cas variant -- move this section to a function specific to cas9
        for row in species_df.itertuples():
            idx = row.Index
            self.species_set.add(idx)  # {0, 1, ..., num_species-1}
            self.species_names.append(row.species_name)  # ['kmarxianus', 'scerevisiae', ...]
            genome_path = os.path.join(input_genome_directory, row.genome_file_name)
            cds_path = os.path.join(input_cds_directory, row.cds_file_name)

            # {0: Species(name='kmarxianus'), 1: Species(name='scerevisiae'), ...}
            self.int_to_species_dict[idx] = Species(
                guide_scorer=scorer_obj,
                guide_container_factory=guide_container_factory,
                name=row.species_name,
                id=idx,
                genome_path=genome_path,
                cds_path=cds_path
                )

        print('Creating coversets for each species.')
        # For each species, get the cas9 guide objects, and create coversets from their sequences
        for species_id, species_object in self.int_to_species_dict.items():
            guide_objects_list = species_object.get_cas9_guides()  # TODO: Add support for another cas variant

            for guide_object in guide_objects_list:
                guide_score = guide_object.score
                sequence = guide_object.sequence

                # TODO: A little hacky whacky
                if sequence in self.cover_sets:
                    self.cover_sets[sequence][1].add(species_id)
                else:
                    self.cover_sets[sequence] = tuple((guide_score, set({species_id})))  # type: ignore

                self.seq_to_guides_dict.setdefault(sequence, list()).append(guide_object)
                

    def get_coversets(self) -> dict[str, tuple[float, set[int]]]:
        return self.cover_sets

    
    def get_species_set(self) -> set[int]:
        return self.species_set

    
    def get_species_names(self) -> list[str]:
        return self.species_names


    def get_species_from_int(self, i: int) -> Species:
        return self.int_to_species_dict[i]

    
    def get_guides_from_seq(self, seq: str) -> list[Guide]:
        return self.seq_to_guides_dict[seq]