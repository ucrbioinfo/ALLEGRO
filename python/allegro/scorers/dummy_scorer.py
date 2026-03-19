from allegro.classes.guide_container import GuideContainer
from allegro.scorers.scorer_base import Scorer
from allegro.utils.guide_finder import GuideFinder

class DummyScorer(Scorer):
    def __init__(self) -> None:
        self.context_toward_five_prime = 0
        self.context_toward_three_prime = 0

    def score_sequence(self, guides_context_list: list[str]) -> list[float]:
        '''
        Assigns all guides in param sequence a score of 1.
        '''
        
        scores = [1.0] * len(guides_context_list)
    
        return scores
