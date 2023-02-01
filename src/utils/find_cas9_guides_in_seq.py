# Script used by ALLEGRO. No need to run it manually.

import re
import numpy
from Bio.Seq import Seq



class GuideFinder:
    def __init__(self) -> None:
        # guide_to_guide_with_context_dict: dict[str, str] = dict()
        pass

    def find_guides_on_both_strands(
        self,
        sequence: str,
        pam_regex: str = r'(?=(GG))',
        toward_5_prime_from_pam: int = 20,
        context_toward_5_prime_from_protospacer: int = 0,
        context_toward_3_prime_from_pam: int = 0,
        ) -> list[str]:

        guides_list: list[str] = list()

        # Store the reverse complement
        sequence_rev_comp = str(Seq(sequence).reverse_complement())

        # Find PAMs on the forward strand
        matches = re.finditer(pam_regex, sequence)
        pam_positions = [match.start() for match in matches]

        for position in pam_positions:
            if (position - toward_5_prime_from_pam >= 0):
                # -20 -1 toward 5' and position-1 because it's NGG and we find GG,
                # we need to discard N and take the next 20 letters.
                guide = sequence[position-toward_5_prime_from_pam-1:position-1]
                guide_with_context = sequence[position-toward_5_prime_from_pam-1-context_toward_5_prime_from_protospacer:position-1+context_toward_3_prime_from_pam]
                
                guides_list.append(guide)
                self.guide_to_guide_with_context_dict[guide] = guide_with_context
                

        # Find PAMs on the reverse comp. strand
        matches = re.finditer(pam_regex, sequence_rev_comp)
        pam_positions = [match.start() for match in matches]

        for position in pam_positions:
            if (position - toward_5_prime_from_pam >= 0):
                guide = sequence_rev_comp[position-toward_5_prime_from_pam:position-1]
                guides_list.append(guide)

        return guides_list


    def find_guides_on_forward_strand(
        self,
        sequence: str,
        pam_regex: str = r'(?=(GG))',
        toward_5_prime: int = 21,
        ) -> list[str]:

        guides_list: list[str] = list()

        # Find PAMs on the forward strand
        matches = re.finditer(pam_regex, sequence)
        pam_positions = [match.start() for match in matches]

        for position in pam_positions:
            if (position - toward_5_prime >= 0):
                # -21 and -1 because it's NGG and we find GG,
                # we need to discard N and take the next 20 letters.
                guide = sequence[position-toward_5_prime:position-1]  
                guides_list.append(guide)

        return guides_list


    def find_guides_on_reverse_strand(
        self,
        sequence: str,
        pam_regex: str = r'(?=(GG))', 
        toward_5_prime: int = 21,
        ) -> list[str]:

        guides_list: list[str] = list()

        # Store the reverse complement
        sequence_rev_comp = str(Seq(sequence).reverse_complement())

        # Find PAMs on the reverse comp. strand
        matches = re.finditer(pam_regex, sequence_rev_comp)
        pam_positions = [match.start() for match in matches]

        for position in pam_positions:
            if (position - toward_5_prime >= 0):
                guide = sequence_rev_comp[position-toward_5_prime:position-1]
                guides_list.append(guide)

        return guides_list


    def find_guides_and_indicate_strand(
        self,
        sequence: str,
        pam_regex: str = r'(?=(GG))',
        toward_5_prime: int = 21,
        ) -> tuple(list[str], list[str], list[int]):
        '''
        ## Args:
            sequence: A nucleotide/DNA sequence to find substrings of guide RNA in.
            pam_regex (optional): The PAM sequence to look for. Currently, only GG is
              supported.
            toward_5_prime (optional): When a PAM is found, how many nucleotides
              to its left do you want to extract?

        ## Returns:
            A tuple of three lists:
            * The first list[str] is a list of the guides found in `sequence`.
            * The seconds list[str] is a list of 'F's and 'R's indicating on
              which strand, forward or reverse, each respective guide resides.
            * The third list[int] shows the location of each guides in `sequence`.
        '''

        guides_list: list[str] = list()
        strands_list: list[str] = list()
        locations_list: list[int] = list()

        # Store the reverse complement
        sequence_rev_comp = str(Seq(sequence).reverse_complement())

        # Find PAMs on the forward strand
        matches = re.finditer(pam_regex, sequence)
        pam_positions = [match.start() for match in matches]

        for position in pam_positions:
            if (position - toward_5_prime >= 0):
                # -21 and -1 because it's NGG and we find GG,
                # we need to discard N and take the next 20 letters.
                guide = sequence[position-toward_5_prime:position-1]  
                
                guides_list.append(guide)
                strands_list.append('F')
                locations_list.append(position-toward_5_prime)

        # Find PAMs on the reverse comp. strand
        matches = re.finditer(pam_regex, sequence_rev_comp)
        pam_positions = [match.start() for match in matches]

        for position in pam_positions:
            if (position - toward_5_prime >= 0):
                guide = sequence_rev_comp[position-toward_5_prime:position-1]

                guides_list.append(guide)
                strands_list.append('R')
                locations_list.append(position-toward_5_prime)

        return guides_list, strands_list, locations_list