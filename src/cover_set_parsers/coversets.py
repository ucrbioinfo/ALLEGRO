import os
import pandas

from classes.guide import Guide
from classes.species import Species
from scorers.scorer_factory import ScorerFactory
from classes.guide_container_factory import GuideContainerFactory


class Coversets:
    def __init__(
        self, 
        scorer_name: str, 
        cas_variant: str,
        guide_source: str,
        input_cds_directory: str,
        input_genome_directory: str, 
        scorer_settings: dict[str, str], 
        input_species_csv_file_path: str, 
        ) -> None:

        self.species_set: set[int] = set()
        self.species_names: list[str] = list()
        self.int_to_species_dict: dict[int, Species] = dict()
        self.seq_to_guides_dict: dict[str, list[Guide]] = dict()

        # To be returned for the solver to use
        self.cover_sets: dict[str, tuple[float, set[int]]] = dict()

        scorer_factory = ScorerFactory()
        scorer = scorer_factory.make_scorer(
            scorer_name=scorer_name, 
            scorer_settings=scorer_settings,
            )

        guide_container_factory = GuideContainerFactory(guide_source=guide_source)

        print('Reading species input file from {path}'.format(path=input_species_csv_file_path))
        species_df = pandas.read_csv(input_species_csv_file_path)

        # Make the species objects
        for row in species_df.itertuples():
            idx = row.Index
            self.species_set.add(idx)  # {0, 1, ..., num_species-1}
            self.species_names.append(row.species_name)  # ['kmarxianus', 'scerevisiae', ...]

            cds_path = os.path.join(input_cds_directory, row.cds_file_name)
            genome_path = os.path.join(input_genome_directory, row.genome_file_name)

            # {0: Species(name='kmarxianus'), 1: Species(name='scerevisiae'), ...}
            self.int_to_species_dict[idx] = Species(
                id=idx,
                cds_path=cds_path,
                name=row.species_name,
                guide_scorer=scorer,
                genome_path=genome_path,
                guide_container_factory=guide_container_factory,
                )

        print('Creating coversets for each species.')

        # For each species, get the cas9 guide objects, and create coversets from their sequences
        for species_id, species_object in self.int_to_species_dict.items():

            guide_objects_list: list[Guide] = list()

            match cas_variant:
                case 'cas9':
                    guide_objects_list = species_object.get_cas9_guides()
                case 'cas12a' | 'cpf1' | _:
                    print('No such cas variant as', cas_variant, 'implemented in coversets.py')
                    raise NotImplementedError

            for guide_object in guide_objects_list:
                guide_score = guide_object.score
                sequence = guide_object.sequence

                # TODO: A little hacky whacky -- cover_sets should be an object?
                if sequence in self.cover_sets:
                    self.cover_sets[sequence][1].add(species_id)
                else:
                    self.cover_sets[sequence] = tuple((guide_score, set({species_id})))  # type: ignore

                self.seq_to_guides_dict.setdefault(sequence, list()).append(guide_object)
                

    def get_species_from_int(self, i: int) -> Species:
        return self.int_to_species_dict[i]

    
    def get_guides_from_seq(self, seq: str) -> list[Guide]:
        return self.seq_to_guides_dict[seq]