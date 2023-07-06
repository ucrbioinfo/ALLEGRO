from classes.guide_container import GuideContainer
from scorers.scorer_base import Scorer
from utils.guide_finder import GuideFinder


class DummyScorer(Scorer):
    def __init__(self, settings: dict) -> None:
        self.pam = settings['pam']
        self.filter_repetitive = settings['filter_repetitive']
        self.protospacer_length = settings['protospacer_length']
        self.context_toward_five_prime = settings['context_toward_five_prime']
        self.context_toward_three_prime = settings['context_toward_three_prime']
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
            filter_repetitive=self.filter_repetitive
        )
        
        scores = [1.0] * len(guides_list)
        
        return guides_list, guides_context_list, strands_list, locations_list, scores


    def save_guides_to_disk(self) -> None:
        return