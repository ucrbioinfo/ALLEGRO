from classes.guide_container import GuideContainer
from scorers.scorer_base import Scorer
from utils.guide_finder import GuideFinder


class DummyScorer(Scorer):
    __slots__ = [
        'pam',
        'protospacer_length',
        'include_repetitive',
        'context_toward_five_prime',
        'context_toward_three_prime',
        'guide_finder'
        ]

    pam: str
    protospacer_length: int
    include_repetitive: bool
    context_toward_five_prime: int
    context_toward_three_prime: int
    guide_finder: GuideFinder

    def __init__(self, settings: dict) -> None:
        self.pam = settings['pam']
        self.include_repetitive = settings['include_repetitive']
        self.protospacer_length = settings['protospacer_length']
        self.context_toward_five_prime = settings['context_toward_five_prime']
        self.context_toward_three_prime = settings['context_toward_three_prime']
        self.guide_finder = GuideFinder()
        

    def score_sequence(
        self,
        guide_container: GuideContainer
        ) -> tuple[list[str], list[str], list[str], list[int], list[int]]:
        '''
        Assigns all guides in param sequence a score of 1
        '''

        (guides_list,
        guides_context_list,
        strands_list,
        locations_list) = self.guide_finder.find_guides_and_indicate_strand(
            pam=self.pam,
            sequence=guide_container.sequence,
            protospacer_length=self.protospacer_length,
            context_toward_five_prime=self.context_toward_five_prime,
            context_toward_three_prime=self.context_toward_three_prime,
            include_repetitive=self.include_repetitive
        )
        
        scores = [1] * len(guides_list)
        
        return guides_list, guides_context_list, strands_list, locations_list, scores
