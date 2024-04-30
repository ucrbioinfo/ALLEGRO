import os
import sys
import pickle
import subprocess
import scipy.stats

from utils.shell_colors import bcolors
from scorers.scorer_base import Scorer
from utils.guide_finder import GuideFinder
from classes.guide_container import GuideContainer


class uCRISPR_scorer(Scorer):
    def __init__(self, settings: dict) -> None:
        self.pam: str = settings['pam']
        self.gc_min: float = settings['gc_min']
        self.gc_max: float = settings['gc_max']
        self.filter_by_gc: bool = settings['filter_by_gc']
        self.score_threshold: int = settings['guide_score_threshold']
        self.patterns_to_exclude: list[str] = settings['patterns_to_exclude']
        self.protospacer_length: int = settings['protospacer_length']
        self.use_secondary_memory: bool = settings['use_secondary_memory']
        self.context_toward_five_prime: int = settings['context_toward_five_prime']
        self.context_toward_three_prime: int = settings['context_toward_three_prime']

        self.guide_finder: GuideFinder = GuideFinder()

        # Relative path from root to uCRISPR C++ program executable.
        self.cpp_program_path = 'src/scorers/uCRISPR/uCRISPR_scorer'  

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
            gc_min=self.gc_min,
            gc_max=self.gc_max,
            filter_by_gc=self.filter_by_gc,
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

                # Check for any errors.
                if stderr:
                    print(f'{bcolors.RED}>{bcolors.RESET} ucrispr_scorer.py: ucrispr_scorer encountered an error: {stderr.decode()}')
                    print('Exiting.')
                    sys.exit(1)

                # Split the output into a list of strings.
                output = output.decode().strip().split('\n')

                # Calculation method taken from CHOPCHOP's code.
                uncached_scores = [scipy.stats.norm.cdf(float(s), loc=11.92658, scale=0.2803797) * 100 for s in output]

                # Update the saved guides cache. Call save to disk later when done with all scorings.
                for idx, s in enumerate(uncached_scores):
                    scores[indices_of_uncached_guides[idx+i]] = s
                    saved_guides[guides_context_list[indices_of_uncached_guides[idx+i]]] = s

            if self.use_secondary_memory == True:
                with open(f'data/cache/{guide_container.species_name}.pickle', 'wb') as f:
                    pickle.dump(saved_guides, f)
        
        # Remove guides with score < threshold.
        marked_for_removal = list()
        for idx, score in enumerate(scores):
            if score < self.score_threshold:
                marked_for_removal.append(idx)

        guides_list = [guide for idx, guide in enumerate(guides_list) if idx not in marked_for_removal]
        guides_context_list = [guide for idx, guide in enumerate(guides_context_list) if idx not in marked_for_removal]
        strands_list = [strand for idx, strand in enumerate(strands_list) if idx not in marked_for_removal]
        locations_list = [loc for idx, loc in enumerate(locations_list) if idx not in marked_for_removal]
        # Truncate decimals. Don't want 24.423871381236. Just take the 24.
        scores = [int(s) for idx, s in enumerate(scores) if idx not in marked_for_removal]
        
        return guides_list, guides_context_list, strands_list, locations_list, scores
    