from classes.guide_container import GuideContainer
from scorers.scorer_base import Scorer
from utils.guide_finder import GuideFinder


class DummyScorer(Scorer):
    def __init__(self, settings: dict) -> None:
        self.pam: str = settings['pam']
        self.gc_min: float = settings['gc_min']
        self.gc_max: float = settings['gc_max']
        self.filter_by_gc: bool = settings['filter_by_gc']
        self.patterns_to_exclude: list[str] = settings['patterns_to_exclude']
        self.protospacer_length: int = settings['protospacer_length']
        self.context_toward_five_prime: int = settings['context_toward_five_prime']
        self.context_toward_three_prime: int = settings['context_toward_three_prime']
        self.guide_finder = GuideFinder()
        

    def score_sequence(
        self,
        guide_container: GuideContainer
        ) -> tuple[list[str], list[str], list[str], list[int], list[float]]:
        '''
        Assigns all guides in param sequence a score of 1
        '''

        (guides_list,
        guides_context_list,
        strands_list,
        locations_list) = self.guide_finder.identify_guides_and_indicate_strand(
            pam=self.pam,
            sequence=guide_container.sequence,
            protospacer_length=self.protospacer_length,
            context_toward_five_prime=self.context_toward_five_prime,
            context_toward_three_prime=self.context_toward_three_prime,
            gc_min=self.gc_min,
            gc_max=self.gc_max,
            filter_by_gc=self.filter_by_gc,
            patterns_to_exclude=self.patterns_to_exclude
        )
        
        scores = [1.0] * len(guides_list)
        
        return guides_list, guides_context_list, strands_list, locations_list, scores
