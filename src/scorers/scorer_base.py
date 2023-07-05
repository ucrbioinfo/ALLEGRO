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
        ) -> tuple[list[str], list[str], list[str], list[int], list[float]]:
        '''
        ## Args:
            * guide_container: Either a Gene or a Chromosome type guide container.
        
        ## Returns:
            A tuple of four lists:
            * The first list[str] is a list of the guides found in `guide_container.sequence`.
            * The second list[str] is a list of guides in `guide_container.sequence` with context around them.
            * The third list[str] is a list of 'F's and 'RC's indicating on
              which strand, forward (F) or reverse complement (RC), each respective guide resides.
            * The fourth list[int] shows the locations of each guides in `guide_container.sequence`.
            * The fifth list[float] indicates the efficiency scores of each guide.
        '''

    @abstractmethod
    def save_guides_to_disk(self) -> None:
        '''
        if self.use_secondary_memory == True, saves scored guides to data/cache/saved_guides.pickle
        '''