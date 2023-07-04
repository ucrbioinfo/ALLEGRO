import subprocess
import scipy.stats

from classes.guide_container import GuideContainer
from scorers.scorer_base import Scorer
from utils.guide_finder import GuideFinder


class uCRISPR_scorer(Scorer):
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
        Identifies guides in guide_container.sequence and pipes the guides to the
        uCRISPR_scorer C++ program. Receives the scores back from the uCRISPR 
        program and returns them.
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

        cpp_program_path = 'src/scorers/uCRISPR/uCRISPR_scorer'  # Relative path from root to 
                                                                 # uCRISPR C++ program executable
        process = subprocess.Popen(cpp_program_path, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        
        for guide in guides_context_list:
            process.stdin.write(guide.encode() + b'\n')  # Send each encoded binary string from 
                                                         # the list to the standard input of the C++ program
        
        process.stdin.close()  # Close the standard input of the C++ program

        output = process.stdout.read().decode()  # Read the output strings from the C++ program
        output = output.strip().split('\n')  # Split the output into a list of strings

        # Caclulation method taken from CHOPCHOP
        scores = [scipy.stats.norm.cdf(float(s.split(' ')[1]), loc=11.92658, scale=0.2803797) * 100 for s in output]

        return guides_list, guides_context_list, strands_list, locations_list, scores