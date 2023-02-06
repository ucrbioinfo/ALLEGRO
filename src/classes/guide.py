class Guide:
    __slots__ = ['pam', 'strand', 'score', 'sequence', 'endonuclease', 
    'genomic_location', 'sequence_with_context', 'guide_container_metadata_dict']
    
    pam: str
    strand: str
    score: float
    sequence: str
    endonuclease: str
    genomic_location: int
    sequence_with_context: str
    guide_container_metadata_dict: dict[str, str | int]

    def __init__(
        self,
        pam: str,
        strand: str,
        score: float,
        sequence: str,
        endonuclease: str,
        genomic_location: int,
        guide_container_metadata_dict: dict[str, str | int],
        sequence_with_context: str = '',
        ) -> None:

        self.pam = pam
        self.score = score
        self.strand = strand
        self.sequence = sequence
        self.endonuclease = endonuclease
        self.genomic_location = genomic_location
        self.sequence_with_context = sequence_with_context
        self.guide_container_metadata_dict = guide_container_metadata_dict

    
    def get_attributes_dict(self) -> dict[str, str | int]:
        guide_attributes = dict({
            'guide_pam': self.pam,
            'guide_score': self.score,
            'guide_strand': self.strand,
            'guide_sequence': self.sequence,
            'guide_endonuclease': self.endonuclease,
            'guide_genomic_location': self.genomic_location,
            'guide_sequence_with_context': self.sequence_with_context,  # Chopchop scorer does not return context
        })

        # merge dictionaries with the | operator
        return guide_attributes | self.guide_container_metadata_dict