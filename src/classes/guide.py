class Guide:
    __slots__ = [
        'strand',
        'score',
        'sequence',
        'has_offtargets',
        'genomic_location',
        'guide_container_metadata_dict'
        ]
    
    strand: str
    score: float
    sequence: str
    has_offtargets: bool
    genomic_location: int
    guide_container_metadata_dict: dict

    def __init__(
        self,
        strand: str,
        score: float,
        sequence: str,
        genomic_location: int,
        guide_container_metadata_dict: dict,
        has_offtargets: bool = False,
        ) -> None:

        self.score = score
        self.strand = strand
        self.sequence = sequence
        self.has_offtargets = has_offtargets
        self.genomic_location = genomic_location
        self.guide_container_metadata_dict = guide_container_metadata_dict

    
    def get_attributes_dict(self) -> dict:
        guide_attributes = dict({
            'guide_score': self.score,
            'guide_strand': self.strand,
            'guide_sequence': self.sequence,
            'guide_has_offtargets': self.has_offtargets,
            'guide_genomic_location': self.genomic_location,
        })

        # merge dictionaries with the | operator
        return guide_attributes | self.guide_container_metadata_dict