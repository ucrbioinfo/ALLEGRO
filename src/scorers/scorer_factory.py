from scorers.scorer_base import Scorer
from scorers.dummy_scorer import DummyScorer
from scorers.chopchop_wrapper import ChopChopWrapper


class ScorerFactory:
    def __init__(self) -> None:
        pass


    def make_scorer(
        self,
        scorer_name: str,
        scorer_settings: dict,
        ) -> Scorer:
        print('Selected scorer:', scorer_name)

        match scorer_name:
            case 'chopchop':
                if scorer_settings == None:
                    print('Error: ChopChopWrapper requires scorer_settings.')
                    raise ValueError

                return ChopChopWrapper(scorer_settings)
            
            case 'dummy':
                if scorer_settings == None:
                    print('Error: DummyScorer requires scorer_settings.')
                    raise ValueError
                
                return DummyScorer(scorer_settings)
            
            case _:
                print('Unknown scorer selected. Aborting.')
                raise ValueError