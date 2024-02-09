# PYTHON LIBS
import sys

# ALLEGRO CYTHON CUSTOM LIBS
import ucrispr_scorer
from libcpp.string cimport string
from libcpp.vector cimport vector

# ---
sys.path.append('..')  # Required to import below

# ALLEGRO PYTHON CUSTOM LIBS
from classes.guide_container import GuideContainer
from scorers.scorer_base import Scorer
from utils.guide_finder import GuideFinder


# Declare the class with cdef
cdef extern from "include/uCRISPR_scorer.h" namespace "uCRISPR_scorer":
    cdef cppclass uCRISPR_scorer:
        uCRISPR_scorer() except +
        
        void hello()

        vector[size_t] score_guides(vector[string] sequences)
        size_t score_guide(string sequence)


cdef class uCRISPR_scorer_cython(Scorer):
    cdef dict __dict__  # Enable cython self.attribute binding.
    cdef uCRISPR_scorer *ucrispr_scorer  # Pointer to C++ class instance.

    def __cinit__(self, settings: dict) -> None:
        self.pam = settings['pam']
        self.filter_repetitive = settings['filter_repetitive']
        self.protospacer_length = settings['protospacer_length']
        self.context_toward_five_prime = settings['context_toward_five_prime']
        self.context_toward_three_prime = settings['context_toward_three_prime']
        self.guide_finder = GuideFinder()
        

    def score_sequence(
        self,
        guide_container: GuideContainer
        ) -> tuple[list[str], list[str], list[str], list[int], list[int]]:

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

        # need 23-mers aka guides_context_list with 0 contexts on both sides
        self.ucrispr_scorer = new uCRISPR_scorer()

        ucrispr_scorer.hello()
        
        return guides_list, guides_context_list, strands_list, locations_list, scores


# CHOPCHOP does this for uCRISPR. TODO: figure out why its doing the extra calculations

# prog = Popen("%s/uCRISPR/uCRISPR -on %s" % (f_p, zhangInputFile), stdout=PIPE, stderr=PIPE, shell=True)
#output = prog.communicate()
#output = output[0].splitlines()
#output = output[1:]
# distribution calculated on 100k random guides
#output = [ss.norm.cdf(float(x.split()[1]), loc=11.92658, scale=0.2803797) for x in output]
#for i, guide in enumerate(results):
#    guide.CoefficientsScore["ZHANG_2019"] = output[i] * 100
#    if args.scoringMethod == "ZHANG_2019":
#        guide.score -= (guide.CoefficientsScore["ZHANG_2019"] / 100) * SCORE['COEFFICIENTS']