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
        guide_scorer: Scorer,
        orthologous_to_gene: str,
        orthologous_to_prot: str,
        ) -> None:

        self.sequence = sequence
        self.gene_name = gene_name
        self.locus_tag = locus_tag
        self.string_id = string_id
        self.integer_id = integer_id
        self.protein_id = protein_id
        self.species = species
        self.ref_species = ref_species
        self.guide_scorer = guide_scorer
        self.orthologous_to_gene = orthologous_to_gene
        self.orthologous_to_prot = orthologous_to_prot

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
                    container=self,
                    endonuclease='cas9',
                    score=guide_strand_score_tuple[2],
                    strand=guide_strand_score_tuple[1],
                    sequence=guide_strand_score_tuple[0],
                    )
                )
            
        return self.cas9_guide_objects


    def get_attributes_dict(self) -> dict:
        return dict({
            'cds_name': self.gene_name,
            'protein_id': self.protein_id,
            'cds_string_id': self.string_id,
            'cds_locus_tag': self.locus_tag,
            'cds_integer_id': self.integer_id,
            'cds_ref_species': self.ref_species,
            'cds_species_name': self.species_name,
            'cds_orthologous_to': self.orthologous_to_gene,
            'cds_orthologous_to_ref_prot_id': self.orthologous_to_prot,
        })