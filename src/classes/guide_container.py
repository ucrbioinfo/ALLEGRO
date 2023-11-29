from classes.guide import Guide
from scorers.scorer_base import Scorer


class GuideContainer:
    __slots__ = [
        'sequence',
        'gene_name',
        'locus_tag',
        'string_id',
        'protein_id',
        'ref_species',
        'species_name',
        'guide_scorer',
        'orthologous_to_gene',
        'orthologous_to_prot'
        ]
    
    sequence: str
    gene_name: str
    locus_tag: str
    string_id: str
    protein_id: str
    ref_species: str
    species_name: str
    guide_scorer: Scorer
    orthologous_to_gene: str
    orthologous_to_prot: str

    def __init__(
        self,
        sequence: str,
        gene_name: str,
        locus_tag: str,
        string_id: str,
        protein_id: str,
        ref_species: str,
        species_name: str,
        guide_scorer: Scorer,
        orthologous_to_gene: str,
        orthologous_to_prot: str,
        ) -> None:

        self.sequence = sequence
        self.gene_name = gene_name
        self.locus_tag = locus_tag
        self.string_id = string_id
        self.protein_id = protein_id
        self.ref_species = ref_species  # Which reference species are we using?
        self.species_name = species_name  # Which species does this gene belong to?
        self.guide_scorer = guide_scorer  # Any assigned scorer in scorers/
        self.orthologous_to_gene = orthologous_to_gene  # Which input gene is this gene orthologous to?
        self.orthologous_to_prot = orthologous_to_prot  # Which input protein_id is this gene orthologous to?
                                                        # This is the protein_id of orthologous_to_gene


    def get_cas9_guides(self) -> list[Guide]:
        cas9_guide_objects: list[Guide] = list()

        (guides_list,
        guides_context_list,
        strands_list, 
        locations_list,
        scores_list) = self.guide_scorer.score_sequence(self)

        for i in range(len(guides_list)):
            cas9_guide_objects.append(Guide(
                score=scores_list[i],
                strand=strands_list[i],
                sequence=guides_list[i],
                genomic_location=locations_list[i],
                # sequence_with_context=guides_context_list[i],
                guide_container_metadata_dict=self.get_attributes_dict(),
                )
            )
        
        return cas9_guide_objects
        

    def get_attributes_dict(self) -> dict:
        return dict({
            'protein_id': self.protein_id,
            'record_name': self.gene_name,
            'record_string_id': self.string_id,
            'record_locus_tag': self.locus_tag,
            'record_ref_species': self.ref_species,
            'record_species_name': self.species_name,
            'record_orthologous_to': self.orthologous_to_gene,
            'record_orthologous_to_ref_prot_id': self.orthologous_to_prot,
        })