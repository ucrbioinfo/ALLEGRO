from utils.shell_colors import bcolors
from scorers.scorer_base import Scorer
from scorers.dummy_scorer import DummyScorer
from scorers.ucrispr_scorer import uCRISPR_scorer
from scorers.chopchop_wrapper import ChopChopWrapper

scorer_names = {
    'dummy': 'dummy',
    'ucrispr': 'uCRISPR',
}

class ScorerFactory:
    def __init__(self) -> None:
        pass


    def make_scorer(
        self,
        scorer_name: str,
        scorer_settings: dict,
        ) -> Scorer:
        print(f'{bcolors.BLUE}>{bcolors.RESET} Selected scorer: {scorer_names[scorer_name]}.')

        match scorer_name.lower():
            case 'chopchop':
                if scorer_settings == None:
                    print(f'{bcolors.RED}> Dev Error{bcolors.RESET}: ChopChopWrapper requires scorer_settings.')
                    raise ValueError

                return ChopChopWrapper(scorer_settings)
            
            case 'dummy':
                if scorer_settings == None:
                    print(f'{bcolors.RED}> Dev Error{bcolors.RESET}: DummyScorer requires scorer_settings.')
                    raise ValueError
                
                return DummyScorer(scorer_settings)
            
            case 'ucrispr':
                if scorer_settings == None:
                    print(f'{bcolors.RED}> Dev Error{bcolors.RESET}: uCRISPR_scorer requires scorer_settings.')
                    raise ValueError
                
                return uCRISPR_scorer(scorer_settings)
            
            case _:
                print(f'{bcolors.RED}> Error{bcolors.RESET}: Unknown scorer selected. Exiting.')
                raise ValueError