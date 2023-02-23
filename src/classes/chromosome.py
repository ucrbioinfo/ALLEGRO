from classes.guide import Guide
from scorers.scorer_base import Scorer
from classes.guide_container import GuideContainer


class Chromosome(GuideContainer):
    __slots__ = ['sequence', 'string_id', 'species_name', 'guide_scorer']

    sequence: str
    string_id: str
    species_name: str
    guide_scorer: Scorer

    def __init__(
        self,
        sequence: str,
        string_id: str,
        species_name: str,
        guide_scorer: Scorer,
        ) -> None:

        self.sequence = sequence
        self.string_id = string_id
        self.species_name = species_name
        self.guide_scorer = guide_scorer


    def get_cas9_guides(self) -> list[Guide]:
        self.cas9_guide_objects: list[Guide] = list()

        (guides_list,
        guides_context_list,
        strands_list, 
        locations_list,
        scores_list) = self.guide_scorer.score_sequence(self)

        for i in range(len(guides_list)):
            self.cas9_guide_objects.append(Guide(
                score=scores_list[i],
                strand=strands_list[i],
                sequence=guides_list[i],
                genomic_location=locations_list[i],
                guide_container_metadata_dict=self.get_attributes_dict(),
                )
            )

        self.sequence = ''

        return self.cas9_guide_objects


    def get_attributes_dict(self) -> dict[str, str]:
        return dict({
            'genome_string_id': self.string_id,
            'genome_species_name': self.species_name,
        })