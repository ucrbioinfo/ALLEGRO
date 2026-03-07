import sys

from allegro.utils.shell_colors import bcolors
from allegro.scorers.scorer_base import Scorer
from allegro.scorers.dummy_scorer import DummyScorer
from allegro.scorers.ucrispr_scorer import uCRISPR_scorer

scorer_names = {
    'dummy': 'dummy',
    'ucrispr': 'uCRISPR',
}

class ScorerFactory:
    def __init__(self) -> None:
        pass

    def make_scorer(self, scorer_name: str) -> Scorer:
        print(f'{bcolors.BLUE}>{bcolors.RESET} Selected scorer: {scorer_names[scorer_name]}.')

        match scorer_name.lower():
            case 'dummy':
                return DummyScorer(guide_settings)
            
            case 'ucrispr':
                return uCRISPR_scorer(guide_settings)
            
            case _:
                print(f'{bcolors.RED}> Error{bcolors.RESET}: Unknown scorer selected. Exiting.')
                sys.exit(1)
                