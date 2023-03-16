from __future__ import annotations
import typing

if typing.TYPE_CHECKING:
    from classes.guide_container import GuideContainer

from abc import ABC, abstractmethod


class Scorer(ABC):
    @abstractmethod
    def score_sequence(
        self,
        guide_container: GuideContainer,
        ) -> tuple[list[str], list[str], list[str], list[int], list[int]]:
        '''
        ## Args:
            * guide_container: Either a Gene or a Chromosome type guide container.
        
        ## Returns:
            A tuple of four lists:
            * The first list[str] is a list of the guides found in `sequence`.
            * The second list[str] is a list of guides in `sequence` with context around them.
            * The third list[str] is a list of '0.0's and '1.0's indicating on
              which strand, forward (0.0) or reverse complement (1.0), each respective guide resides.
            * The fourth list[int] shows the locations of each guides in `sequence`.
            * The fifth list[int] indicates the efficiency scores of each guide.
        '''