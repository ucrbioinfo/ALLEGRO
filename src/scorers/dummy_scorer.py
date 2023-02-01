from classes.guide_container import GuideContainer
from scorers.scorer_base import Scorer
from utils.find_cas9_guides_in_seq import GuideFinder


class DummyScorer(Scorer):
    def __init__(self) -> None:
        super().__init__()
        

    def score_sequence(
        self,
        guide_container: GuideContainer
        ) -> tuple[list[str], list[str], list[int], list[float]]:
        '''
        Assigns all guides in param `guide_container` a score of 1.0
        '''

        sequence = guide_container.sequence
        gf = GuideFinder()

        (guides_list, 
        strands_list, 
        locations_list) = gf.find_guides_and_indicate_strand(sequence=sequence)
        
        one_scores = [1.0] * len(guides_list)
        
        return guides_list, strands_list, locations_list, one_scores
