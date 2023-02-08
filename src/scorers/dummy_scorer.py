from classes.guide_container import GuideContainer
from scorers.scorer_base import Scorer
from utils.find_guides_in_seq import find_guides_and_indicate_strand


class DummyScorer(Scorer):
    __slots__ = ['pam', 'protospacer_length',
    'context_toward_five_prime', 'context_toward_three_prime']

    pam: str
    protospacer_length: str
    context_toward_five_prime: str
    context_toward_three_prime: str

    def __init__(self, settings: dict[str, str]) -> None:
        super().__init__(settings=settings)

        self.pam = settings['pam']
        self.protospacer_length = settings['protospacer_length']
        self.context_toward_five_prime = settings['context_toward_five_prime']
        self.context_toward_three_prime = settings['context_toward_three_prime']
        

    def score_sequence(
        self,
        guide_container: GuideContainer
        ) -> tuple[list[str], list[str], list[str], list[str], list[int], list[float]]:
        '''
        Assigns all guides in param `guide_container` a score of 1.0
        '''

        (guides_list, 
        guides_context_list,
        strands_list, 
        locations_list) = find_guides_and_indicate_strand(
            pam=self.pam,
            sequence=guide_container.sequence,
            protospacer_length=self.protospacer_length,
            context_toward_five_prime=self.context_toward_five_prime,
            context_toward_three_prime=self.context_toward_three_prime
        )
        
        one_scores = [1.0] * len(guides_list)
        
        return guides_list, guides_context_list, strands_list, locations_list, one_scores
