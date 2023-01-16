from scorers.scorer_base import Scorer
from scorers.dummy_scorer import DummyScorer
from scorers.chopchop_wrapper import ChopChopWrapper
from scorers.deepguide_wrapper import DeepGuideWrapper
from scorers.random_scorer import RandomScorer


class ScorerFactory:
    def __init__(self) -> None:
        pass


    def make_scorer(self, scorer_name: str, scorer_settings: dict[str, str]) -> Scorer:
        print('Selected scorer:', scorer_name)

        match scorer_name:
            case 'chopchop':
                return ChopChopWrapper(scorer_settings)
            case 'random':
                return RandomScorer()  # Assigns all guides a random score between [0, 10)
            case _:
                print('Unknown/dummy scorer selected. Assigning all guides a default score of 1.0')
                return DummyScorer()