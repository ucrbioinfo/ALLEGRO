from abc import ABC, abstractmethod

class CoversetsBase(ABC):
    @abstractmethod
    def get_guide_attributes_dicts_from_seq(self, seqsseqs_iterable: list[str] | set[str]) -> list[dict]:
        '''
        ## Args:
            * seqs_iterable: a list or set object containing string sequences.
        
        ## Returns:
            * a list containing dictionaries of Guide attributes.
        '''