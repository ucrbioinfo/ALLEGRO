from classes.guide import Guide
from scorers.scorer_base import Scorer
from classes.guide_container import GuideContainer


class Chromosome(GuideContainer):
    __slots__ = ['sequence', 'string_id', 'integer_id', 'species_name',
    'guide_scorer', 'cas9_guide_objects']

    sequence: str
    string_id: str
    integer_id: int
    species_name: str
    guide_scorer: Scorer
    cas9_guide_objects: list[Guide]

    def __init__(
        self,
        sequence: str,
        string_id: str,
        integer_id: int,
        species_name: str,
        guide_scorer: Scorer,
        ) -> None:

        self.sequence = sequence
        self.string_id = string_id
        self.integer_id = integer_id
        self.species_name = species_name
        self.guide_scorer = guide_scorer

        self.cas9_guide_objects = list()
        

    def get_cas9_guides(self) -> list[Guide]:
        if len(self.cas9_guide_objects) == 0:

            (guides_list,
            guides_context_list,
            strands_list, 
            locations_list,
            scores_list) = self.guide_scorer.score_sequence(self)

            for i in range(len(guides_list)):
                self.cas9_guide_objects.append(Guide(
                    pam='GG',
                    container=self,
                    endonuclease='cas9',
                    score=scores_list[i],
                    strand=strands_list[i],
                    sequence=guides_list[i],
                    genomic_location=locations_list[i],
                    sequence_with_context=guides_context_list[i]
                    )
                )

        self.sequence = ''

        return self.cas9_guide_objects


    def get_attributes_dict(self) -> dict:
        return dict({
            'genome_string_id': self.string_id,
            'genome_integer_id': self.integer_id,
            'genome_species_name': self.species_name,
        })