class Guide:
    __slots__ = ['strand', 'score', 'sequence',
    'genomic_location', 'guide_container_metadata_dict']
    
    strand: str
    score: float
    sequence: str
    genomic_location: int
    guide_container_metadata_dict: dict[str, str]

    def __init__(
        self, # HASH TABLE instead of dict # see if you can remove irrelevant guides
        strand: str,
        score: float,
        sequence: str, # 2 bits per symbol
        genomic_location: int,
        guide_container_metadata_dict: dict[str, str],
        ) -> None:

        self.score = score
        self.strand = strand
        self.sequence = sequence
        self.genomic_location = genomic_location
        self.guide_container_metadata_dict = guide_container_metadata_dict

    
    def get_attributes_dict(self) -> dict[str, str | float | int]:
        guide_attributes = dict({
            'guide_score': self.score,
            'guide_strand': self.strand,
            'guide_sequence': self.sequence,
            'guide_genomic_location': self.genomic_location,
        })

        # merge dictionaries with the | operator
        return guide_attributes | self.guide_container_metadata_dict