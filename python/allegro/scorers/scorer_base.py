from __future__ import annotations
import typing

if typing.TYPE_CHECKING:
    from allegro.classes.guide_container import GuideContainer

from abc import ABC, abstractmethod

class Scorer(ABC):
    @abstractmethod
    def score_sequence(self, guides_context_list: list[str]) -> list[float]:
        '''
        ## Args:
            * guide_container: Either a Gene or a Chromosome type guide container.
        
        ## Returns:
            * list[float] indicates the efficiency scores of each guide.
        '''
        