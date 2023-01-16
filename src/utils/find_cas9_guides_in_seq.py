# Script used by ALLEGRO. No need to run it manually.

import re
import numpy
from Bio.Seq import Seq


def find_guides_on_both_strands(sequence: str, pam_regex: str = r'(?=(GG))', toward_5_prime: int = 21) -> list[str]:
    guides_list: list[str] = list()

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

    # Find PAMs on the reverse comp. strand
    matches = re.finditer(pam_regex, sequence_rev_comp)
    pam_positions = [match.start() for match in matches]

    for position in pam_positions:
        if (position - toward_5_prime >= 0):
            guide = sequence_rev_comp[position-toward_5_prime:position-1]
            guides_list.append(guide)

    return guides_list


def find_guides_on_forward_strand(sequence: str, pam_regex: str = r'(?=(GG))', toward_5_prime: int = 21) -> list[str]:
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


def find_guides_on_reverse_strand(sequence: str, pam_regex: str = r'(?=(GG))', toward_5_prime: int = 21) -> list[str]:
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


def find_guides_and_indicate_strand(sequence: str, pam_regex: str = r'(?=(GG))', toward_5_prime: int = 21) -> list[tuple[str, str, float]]:
    guides_list: list[tuple[str, str, float]] = list()

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
            guides_list.append((guide, 'F', 1.0))

    # Find PAMs on the reverse comp. strand
    matches = re.finditer(pam_regex, sequence_rev_comp)
    pam_positions = [match.start() for match in matches]

    for position in pam_positions:
        if (position - toward_5_prime >= 0):
            guide = sequence_rev_comp[position-toward_5_prime:position-1]
            guides_list.append((guide, 'R', 1.0))

    return guides_list



def find_guides_and_indicate_strand_random_scores(sequence: str, pam_regex: str = r'(?=(GG))', toward_5_prime: int = 21) -> list[tuple[str, str, float]]:
    guides_list: list[tuple[str, str, float]] = list()

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

            guides_list.append((guide, 'F', numpy.random.rand(1)[0] * 10))

    # Find PAMs on the reverse comp. strand
    matches = re.finditer(pam_regex, sequence_rev_comp)
    pam_positions = [match.start() for match in matches]

    for position in pam_positions:
        if (position - toward_5_prime >= 0):
            guide = sequence_rev_comp[position-toward_5_prime:position-1]
            guides_list.append((guide, 'R', numpy.random.rand(1)[0] * 10))

    return guides_list