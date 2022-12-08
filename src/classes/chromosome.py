from __future__ import annotations
import typing

if typing.TYPE_CHECKING:
    from classes.species import Species

from classes.guide import Guide
from classes.guide_container import GuideContainer
from scorers.scorer_base import Scorer


class Chromosome(GuideContainer):
    def __init__(
        self,
        species: Species,
        guide_scorer: Scorer,
        sequence: str,
        string_id: str,
        integer_id: int,
        ) -> None:

        self.species = species
        self.guide_scorer = guide_scorer
        self.sequence = sequence
        self.string_id = string_id
        self.integer_id = integer_id

        self.guide_objects: list[Guide] = list()


    def get_cas9_guides(self) -> list[Guide]:
        guide_strand_score_tuple_list = self.guide_scorer.score_sequence(self)

        for guide_strand_score_tuple in guide_strand_score_tuple_list:
            self.guide_objects.append(Guide(
                pam='GG',
                endonuclease='cas9',
                score=guide_strand_score_tuple[2],
                strand=guide_strand_score_tuple[1],
                sequence=guide_strand_score_tuple[0],
                container=self,
                )
            )
        
        return self.guide_objects


    def get_species(self) -> Species:
        return self.species


    def get_sequence(self) -> str:
        return self.sequence


    def get_string_id(self) -> str:
        return self.string_id

    
    def get_integer_id(self) -> int:
        return self.integer_id


    def get_guides(self) -> list[Guide]:
        return self.guide_objects


    def get_attributes_dict(self) -> dict:
        return dict({
            'species': self.species,
            'sequence': self.sequence,
            'string_id': self.string_id,
            'integer_id': self.integer_id
        })


    def get_gene_name(self) -> str:
        return ''