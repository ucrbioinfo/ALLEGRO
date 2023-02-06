from classes.guide import Guide
from scorers.scorer_base import Scorer
from classes.guide_container import GuideContainer


class Gene(GuideContainer):
    __slots__ = ['sequence', 'gene_name', 'locus_tag', 'string_id', 'integer_id',
    'protein_id', 'species_name', 'ref_species', 'guide_scorer', 'orthologous_to_gene',
    'orthologous_to_prot', 'cas9_guide_objects']
    
    sequence: str
    gene_name: str
    locus_tag: str
    string_id: str
    integer_id: int
    protein_id: str
    species_name: str
    ref_species: str
    guide_scorer: Scorer
    orthologous_to_gene: str
    orthologous_to_prot: str
    cas9_guide_objects: list[Guide]

    def __init__(
        self,
        sequence: str,
        gene_name: str,
        locus_tag: str,
        string_id: str,
        integer_id: int,
        protein_id: str,
        species_name: str,
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
        self.ref_species = ref_species  # Which reference species are we using? 
                                        # This is found in utils/find_orthogroup_config
                                        # under the species attribute.
        self.species_name = species_name  # Which species does this gene belong to?
        self.guide_scorer = guide_scorer  # Any assigned scorer in scorers/
        self.orthologous_to_gene = orthologous_to_gene  # Which input gene is this gene orthologous to?
        self.orthologous_to_prot = orthologous_to_prot  # Which input protein_id is this gene orthologous to?
                                                        # This is the protein_id of the orthologous_to_gene

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