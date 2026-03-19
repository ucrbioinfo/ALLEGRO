from allegro.classes.guide import Guide
from allegro.utils.guide_finder import GuideFinder

class GuideContainer:
    __slots__ = [
        'sequence',
        'gene_name',
        'locus_tag',
        'string_id',
        'protein_id',
        'ref_species',
        'species_name',
        'orthologous_to_gene',
        'orthologous_to_prot',
        'guide_finder']
    
    sequence: str
    gene_name: str
    locus_tag: str
    string_id: str
    protein_id: str
    ref_species: str
    species_name: str
    orthologous_to_gene: str
    orthologous_to_prot: str
    guide_finder: GuideFinder
    
    def __init__(
        self,
        sequence: str,
        gene_name: str,
        locus_tag: str,
        string_id: str,
        protein_id: str,
        ref_species: str,
        species_name: str,
        orthologous_to_gene: str,
        orthologous_to_prot: str) -> None:
        
        self.guide_finder = GuideFinder()
        self.sequence = str(sequence).upper()
        self.gene_name = gene_name
        self.locus_tag = locus_tag
        self.string_id = string_id
        self.protein_id = protein_id
        self.ref_species = ref_species  # Which reference species are we using?
        self.species_name = species_name  # Which species does this gene belong to?
        self.orthologous_to_gene = orthologous_to_gene  # Which input gene is this gene orthologous to?
        self.orthologous_to_prot = orthologous_to_prot  # Which input protein_id is this gene orthologous to?
                                                        # This is the protein_id of orthologous_to_gene
    def get_guides(self) -> list[Guide]:
        guide_objects: list[Guide] = list()

        (guides_list,
        guides_context_list,
        strands_list,
        locations_list,
        scores_list) = self.guide_finder.identify_guides_and_indicate_strand(sequence=self.sequence)

        # Create Guide objects
        for i in range(len(guides_list)):
            guide_objects.append(Guide(
                score=scores_list[i],
                strand=strands_list[i],
                sequence=guides_list[i],
                genomic_location=locations_list[i],
                guide_container_metadata_dict=self.get_attributes_dict()))
        
        return guide_objects

    def get_attributes_dict(self) -> dict:
        return dict({
            'protein_id': self.protein_id,
            'record_name': self.gene_name,
            'record_string_id': self.string_id,
            'record_locus_tag': self.locus_tag,
            'record_ref_species': self.ref_species,
            'record_species_name': self.species_name,
            'record_orthologous_to': self.orthologous_to_gene,
            'record_orthologous_to_ref_prot_id': self.orthologous_to_prot})
    