from __future__ import annotations
import typing

if typing.TYPE_CHECKING:
    from classes.species import Species

from classes.guide import Guide
from scorers.scorer_base import Scorer
from classes.guide_container import GuideContainer


class Gene(GuideContainer):
    def __init__(
        self,
        sequence: str,
        gene_name: str,
        locus_tag: str,
        string_id: str,
        integer_id: int,
        protein_id: str,
        species: Species,
        ref_species: str,
        orthologous_to: str,
        guide_scorer: Scorer,
        ) -> None:

        self.sequence = sequence
        self.gene_name = gene_name
        self.locus_tag = locus_tag
        self.string_id = string_id
        self.integer_id = integer_id
        self.protein_id = protein_id
        self.species = species
        self.ref_species = ref_species
        self.orthologous_to = orthologous_to
        self.guide_scorer = guide_scorer    

        self.cas9_guide_objects: list[Guide] = list()


    @property
    def species_name(self) -> str:
        return self.species.name


    def get_cas9_guides(self) -> list[Guide]:
        if len(self.cas9_guide_objects) == 0:
            guide_strand_score_tuple_list = self.guide_scorer.score_sequence(self)

            for guide_strand_score_tuple in guide_strand_score_tuple_list:
                self.cas9_guide_objects.append(Guide(
                    pam='GG',
                    endonuclease='cas9',
                    score=guide_strand_score_tuple[2],
                    strand=guide_strand_score_tuple[1],
                    sequence=guide_strand_score_tuple[0],
                    container=self,
                    )
                )
            
        return self.cas9_guide_objects


    def get_attributes_dict(self) -> dict:
        return dict({
            'gene_name': self.gene_name,
            'protein_id': self.protein_id,
            'gene_string_id': self.string_id,
            'gene_locus_tag': self.locus_tag,
            'gene_integer_id': self.integer_id,
            'gene_ref_species': self.ref_species,
            'gene_species_name': self.species_name,
            'gene_orthologous_to': self.orthologous_to,
        })