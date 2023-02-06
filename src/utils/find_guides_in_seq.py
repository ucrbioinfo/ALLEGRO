# Function imported by ALLEGRO. No need to run it manually.
import re
from Bio.Seq import Seq


def find_guides_and_indicate_strand(
    pam: str,
    sequence: str,
    protospacer_length: int,
    context_toward_five_prime: int,
    context_toward_three_prime: int,
    ) -> tuple[list[str], list[str], list[str], list[int]]:
    '''
    ## Args:
        pam: The protospacer adjacent motif to look for, e.g., 'NGG'.
        sequence: A nucleotide/DNA sequence to find substrings of guide RNA in.
        protospacer_length: The length of the string toward 5-prime of the PAM.
        context_toward_five_prime: The number of nucleotides to extract toward
            the 5-prime after the protospacer.
        context_toward_three_prime: The number of nucleotides to extract toward
            the 3-prime after (and excluding) the PAM.

    ## Returns:
        A tuple of three lists:
        * The first list[str] is a list of the guides found in `sequence`.
        * The second list[str] is a list of the guides with their context around them.
        * The third list[str] is a list of 'F's and 'R's indicating on
            which strand, forward or reverse, each respective guide resides.
        * The third list[int] shows the location of each guides in `sequence`.

    ## Example
        Input: find_guides_and_indicate_strand(
            pam='NGG',
            sequence='AAAAAACCAAAAAGGTTTTTT',
            protospacer_length=3,
            context_toward_five_prime=2,
            context_toward_three_prime=2,
        )
        Output:
        (['AAA', 'TTT'], ['CAAAAAGGTT', 'CTTTTTGGTT'], ['F', 'R'], [12, 12])
    '''

    # Store the reverse complement
    sequence_rev_comp = str(Seq(sequence).reverse_complement())

    guides_list: list[str] = list()
    strands_list: list[str] = list()
    locations_list: list[int] = list()
    guides_context_list: list[str] = list()

    pam_dict = {
        'NGG': r'(?=(AGG))|(?=(CGG))|(?=(TGG))|(?=(GGG))',
    }

    pam_regex = r''
    if pam in pam_dict:
        pam_regex = pam_dict[pam]
    else:
        pam_regex = r'(?=({PAM}))'.format(PAM=pam)


    def find_matches(seq: str, strand: str) -> None:
        matches = re.finditer(pam_regex, seq)  # Find PAMs on seq
        pam_positions = [match.start() for match in matches]

        for position in pam_positions:
            if (position - protospacer_length >= 0):
                guide = seq[position-protospacer_length:position]
                guide_with_context = ''

                # There is enough context on both sides of the PAM
                if position-protospacer_length-context_toward_five_prime >= 0 and position+len(pam)+context_toward_three_prime < len(seq) + 1:
                    guide_with_context = seq[position-protospacer_length-context_toward_five_prime:position+len(pam)+context_toward_three_prime]
                
                # There is not enough context on the right side -- extract as many chars as possible
                elif position-protospacer_length-context_toward_five_prime >= 0:
                    guide_with_context = seq[position-protospacer_length-context_toward_five_prime:]
                
                # There is not enough context on the left side
                #  extract from the start of the string to the specified context on the right side
                elif position+len(pam)+context_toward_three_prime < len(seq) + 1:
                    guide_with_context = seq[:position+len(pam)+context_toward_three_prime]
                
                # There is not enough context on either side
                else:
                    guide_with_context = seq
                
                guides_list.append(guide)
                strands_list.append(strand)
                guides_context_list.append(guide_with_context)
                locations_list.append(position)

    find_matches(seq=sequence, strand='F')
    find_matches(seq=sequence_rev_comp, strand='R')

    return guides_list, guides_context_list, strands_list, locations_list