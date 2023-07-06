import os
import pickle
import subprocess
import scipy.stats

from classes.guide_container import GuideContainer
from scorers.scorer_base import Scorer
from utils.guide_finder import GuideFinder
from utils.shell_colors import bcolors


class uCRISPR_scorer(Scorer):
    def __init__(self, settings: dict) -> None:
        self.pam: str = settings['pam']
        self.filter_repetitive: bool = settings['filter_repetitive']
        self.protospacer_length: int = settings['protospacer_length']
        self.use_secondary_memory: bool = settings['use_secondary_memory']
        self.context_toward_five_prime: int = settings['context_toward_five_prime']
        self.context_toward_three_prime: int = settings['context_toward_three_prime']

        self.guide_finder: GuideFinder = GuideFinder()
        self.saved_guides: dict[str, float] = dict()

        if self.use_secondary_memory == True:
            if not os.path.exists('data/cache/'):
                os.mkdir('data/cache/')

            elif os.path.exists('data/cache/saved_guides.pickle'):
                with open('data/cache/saved_guides.pickle', 'rb') as f:
                    self.saved_guides = pickle.load(f)
                    print(f'{bcolors.BLUE}>{bcolors.RESET} Loaded cached guide scores.')


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

        scores: list = list()
        indices_of_uncached_guides: list[int] = list()

        # If there are guides which we do not have the scores for, save their indices.
        # Score them with the C++ uCRISPR after this block.
        for idx, guide in enumerate(guides_context_list):
            score = self.saved_guides.get(guide)
            scores.append(score)
            
            if score == None:
                indices_of_uncached_guides.append(idx)

        if len(indices_of_uncached_guides) > 0:
            cpp_program_path = 'src/scorers/uCRISPR/uCRISPR_scorer'  # Relative path from root to 
                                                                     # uCRISPR C++ program executable.
            process = subprocess.Popen(cpp_program_path, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            
            for i in indices_of_uncached_guides:
                process.stdin.write(guides_context_list[i].encode() + b'\n')  # Send each encoded binary string from 
                                                                              # the list to the standard input of the C++ program.
            
            process.stdin.close()  # Close the standard input of the C++ program.

            output = process.stdout.read().decode()  # Read the output strings from the C++ program.
            output = output.strip().split('\n')      # Split the output into a list of strings.

            # Caclulation method taken from CHOPCHOP's code.
            uncached_scores = [scipy.stats.norm.cdf(float(s.split(' ')[1]), loc=11.92658, scale=0.2803797) * 100 for s in output]

            # Update the saved guides cache. Call save to disk later when done with all scorings.
            for idx, s in enumerate(uncached_scores):
                scores[indices_of_uncached_guides[idx]] = s
                self.saved_guides[guides_context_list[indices_of_uncached_guides[idx]]] = s

        return guides_list, guides_context_list, strands_list, locations_list, scores
    

    def save_guides_to_disk(self) -> None:
        if self.use_secondary_memory == True:
            with open('data/cache/saved_guides.pickle', 'wb') as f:
                pickle.dump(self.saved_guides, f)