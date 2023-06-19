from scorers.scorer_base import Scorer
from scorers.dummy_scorer import DummyScorer
from scorers.chopchop_wrapper import ChopChopWrapper
from scorers.ucrispr_scorer import uCRISPR_scorer


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
                    print('Dev Error: ChopChopWrapper requires scorer_settings.')
                    raise ValueError

                return ChopChopWrapper(scorer_settings)
            
            case 'dummy':
                if scorer_settings == None:
                    print('Dev Error: DummyScorer requires scorer_settings.')
                    raise ValueError
                
                return DummyScorer(scorer_settings)
            
            case 'uCRISPR' | 'ucrispr':
                if scorer_settings == None:
                    print('Dev Error: uCRISPR_scorer requires scorer_settings.')
                    raise ValueError
                
                return uCRISPR_scorer(scorer_settings)
            
            case _:
                print('Unknown scorer selected. Aborting.')
                raise ValueError