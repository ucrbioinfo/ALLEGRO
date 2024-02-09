import os
import sys
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
        self.patterns_to_exclude: bool = settings['patterns_to_exclude']
        self.protospacer_length: int = settings['protospacer_length']
        self.use_secondary_memory: bool = settings['use_secondary_memory']
        self.context_toward_five_prime: int = settings['context_toward_five_prime']
        self.context_toward_three_prime: int = settings['context_toward_three_prime']

        self.guide_finder: GuideFinder = GuideFinder()

        self.cpp_program_path = 'src/scorers/uCRISPR/uCRISPR_scorer'  # Relative path from root to 
                                                                      # uCRISPR C++ program executable.

        if self.use_secondary_memory == True:
            if not os.path.exists('data/cache/'):
                os.mkdir('data/cache/')


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
            patterns_to_exclude=self.patterns_to_exclude
        )

        scores: list = list()
        saved_guides: dict[str, float] = dict()
        indices_of_uncached_guides: list[int] = list()

        if self.use_secondary_memory == True:
            if os.path.exists(f'data/cache/{guide_container.species_name}.pickle'):
                with open(f'data/cache/{guide_container.species_name}.pickle', 'rb') as f:
                    saved_guides = pickle.load(f)

        # If there are guides which we do not have the scores for, save their indices.
        # Score them with the C++ uCRISPR after this block.
        for idx, guide in enumerate(guides_context_list):
            score = saved_guides.get(guide)
            scores.append(score)
            
            if score == None:
                indices_of_uncached_guides.append(idx)

        if len(indices_of_uncached_guides) > 0:
            for i in range(0, len(indices_of_uncached_guides), 1000):
                payload = b''
                for j in indices_of_uncached_guides[i:i+1000]:
                    payload += guides_context_list[j].encode() + b'\n'

                process = None
                try:
                    process = subprocess.Popen(self.cpp_program_path, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                except PermissionError:
                    print(f'{bcolors.RED}> Error{bcolors.RESET}: {bcolors.ORANGE}ALLEGRO{bcolors.RESET} could not initiate the uCRISPR executable for guide scoring. Please manually navigate to src/scorers/uCRISPR/ and give the appropriate execution permission to this file: uCRISPR_scorer ($ chmod +x uCRISPR_scorer) and try again. Exiting.')
                    sys.exit(1)

                output, stderr = process.communicate(payload)

                # Check for any errors
                if stderr:
                    print(f'{bcolors.RED}>{bcolors.RESET} ucrispr_scorer.py: ucrispr_scorer encountered an error: {stderr.decode()}')
                    print('Exiting.')
                    sys.exit(1)

                output = output.decode().strip().split('\n')      # Split the output into a list of strings.

                # Calculation method taken from CHOPCHOP's code.
                uncached_scores = [scipy.stats.norm.cdf(float(s), loc=11.92658, scale=0.2803797) * 100 for s in output]

                # Update the saved guides cache. Call save to disk later when done with all scorings.
                for idx, s in enumerate(uncached_scores):
                    scores[indices_of_uncached_guides[idx+i]] = s
                    saved_guides[guides_context_list[indices_of_uncached_guides[idx+i]]] = s

            if self.use_secondary_memory == True:
                with open(f'data/cache/{guide_container.species_name}.pickle', 'wb') as f:
                    pickle.dump(saved_guides, f)

        return guides_list, guides_context_list, strands_list, locations_list, scores
    