import os
import pandas

from classes.guide import Guide
from classes.species import Species
from scorers.scorer_factory import ScorerFactory
from coverset_parsers.coversets_base import CoversetsBase
from classes.guide_container_factory import GuideContainerFactory


class CoversetsRAM(CoversetsBase):
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
        self.coversets: dict[str, tuple[float, set[int]]] = dict() # To be returned for the solver to use
        self.seq_to_guides_attributes_dict: dict[str, list[dict]] = dict()

        scorer_factory = ScorerFactory()
        scorer = scorer_factory.make_scorer(
            scorer_name=scorer_name,
            scorer_settings=scorer_settings,
        )

        guide_container_factory = GuideContainerFactory()

        print('Reading species input file from {path}'.format(path=input_species_csv_file_path))
        species_df = pandas.read_csv(input_species_csv_file_path)

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

            match cas_variant:
                case 'cas9':
                    guide_objects_list = species_object.get_cas9_guides()
                
                case 'cas12a' | 'cpf1' | _:
                    print('No such cas variant as', cas_variant, 'implemented in coversets.py')
                    raise NotImplementedError

            for guide_object in guide_objects_list:
                guide_sequence = guide_object.sequence

                guides_attributes_list.append(guide_object.get_attributes_dict())
                self.seq_to_guides_attributes_dict.setdefault(guide_sequence, list()).append(guide_object.get_attributes_dict())

                # If we have multiple guides with the same sequence from different species,
                # assign their average score to a unique representative guide
                # that is input to the solver.
                if guide_sequence in self.coversets:
                    average_score = self.get_average_score(
                        self.seq_to_guides_attributes_dict[guide_sequence],
                    )
                    self.coversets[guide_sequence][1].add(idx)
                    old_set = self.coversets[guide_sequence][1]
                    self.coversets[guide_sequence] = (average_score, old_set)  # Update the old avg
                else:
                    self.coversets[guide_sequence] = tuple((guide_object.score, set({idx})))  # type: ignore
            
            print('Done with', idx + 1, 'species...')

        print('Created coversets for all species.')


    def get_average_score(
        self,
        guides_list: list[dict],
        ) -> float:
        '''
        ## Args:
            list of Guide attribute dictionaries.

        ## Returns:
            A float value. The average guide_score from all dicts.
        '''
        avg = 0.0

        for d in guides_list:
            avg += d['guide_score']

        return avg / len(guides_list)


    def get_guide_attributes_dicts_from_seq(self, seqs_iterable: list[str] | set[str]) -> list[dict]:
        results: list[dict] = list()
        for seq in seqs_iterable:
            results.extend(self.seq_to_guides_attributes_dict[seq])

        return results