from scorers.scorer_base import Scorer
from scorers.dummy_scorer import DummyScorer
from scorers.chopchop_wrapper import ChopChopWrapper


class ScorerFactory:
    def __init__(self) -> None:
        pass


    def make_scorer(
        self,
        scorer_name: str,
        scorer_settings: dict[str, str] = None,
        ) -> Scorer:
        print('Selected scorer:', scorer_name)

        match scorer_name:
            case 'chopchop':
                if scorer_settings == None:
                    print('Error: CHOPCHOP requires scorer_settings.')
                    raise ValueError

                return ChopChopWrapper(scorer_settings)
            
            case 'dummy':
                return DummyScorer()
            
            case _:
                print('Unknown/dummy scorer selected. Assigning all guides a default score of 1.0')
                return DummyScorer()