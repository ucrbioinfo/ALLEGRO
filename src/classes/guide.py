from __future__ import annotations
import typing

if typing.TYPE_CHECKING:
    from classes.species import Species
    from classes.guide_container import GuideContainer


class Guide:
    __slots__ = ['pam', 'strand', 'score', 'sequence', 'endonuclease', 
    'genomic_location', 'container', 'sequence_with_context']
    
    pam: str
    strand: str
    score: float
    sequence: str
    endonuclease: str
    genomic_location: int
    container: GuideContainer
    sequence_with_context: str

    def __init__(
        self,
        pam: str,
        strand: str,
        score: float,
        sequence: str,
        endonuclease: str,
        genomic_location: int,
        container: GuideContainer,
        sequence_with_context: str = '',
        ) -> None:

        self.pam = pam
        self.score = score
        self.strand = strand
        self.sequence = sequence
        self.container = container
        self.endonuclease = endonuclease
        self.genomic_location = genomic_location
        self.sequence_with_context = sequence_with_context


    @property
    def species_name(self) -> str:
        return self.container.species.name


    @property
    def location(self) -> str:
        return self.container.string_id


    @property
    def guide_in_gene_orthologous_to(self) -> str:
        return self.container.orthologous_to

    
    def get_attributes_dict(self) -> dict:
        guide_attributes = dict({
            'guide_pam': self.pam,
            'guide_score': self.score,
            'guide_strand': self.strand,
            'guide_sequence': self.sequence,
            'guide_endonuclease': self.endonuclease,
            'guide_genomic_location': self.genomic_location,
            'guide_sequence_with_context': self.sequence_with_context,  # Chopchop scorer does not return context
        })

        # merge dictionaries
        return guide_attributes | self.container.get_attributes_dict()


    def get_species(self) -> Species:
        '''
        ## Returns
            A Species object. To get the species name as a str, use the class'
             `species_name`` property.
        '''
        return self.container.get_species()