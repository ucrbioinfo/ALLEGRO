from __future__ import annotations

from abc import ABC, abstractmethod

class Scorer(ABC):
    @abstractmethod
    def score_sequence(self, guides_context_list: list[str]) -> list[float]:
        '''
        ## Args:
            * guides_context_list: a list of guide strings with extra base context on their 5' and 3'.
        
        ## Returns:
            * list[float] indicates the efficiency scores of each guide.
        '''
        