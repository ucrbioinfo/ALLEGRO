# Functions imported by ALLEGRO. No need to run it manually.
import re
from Bio import SeqIO
from Bio.Seq import Seq


def count_kmers(sequence, k):
    kmers = dict()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers[kmer] = kmers.setdefault(kmer, 0) + 1

    return kmers


class GuideFinder:
    __slots__ = [
        'pam_dict',
        'num_containing_fourmers_removed',
        'num_containing_fivemers_removed',
        'num_more_than_four_twomers_removed'
        ]

    pam_dict: dict
    num_containing_fourmers_removed: int
    num_containing_fivemers_removed: int
    num_more_than_four_twomers_removed: int

    def __init__(self) -> None:
        self.pam_dict = {
            'NGG': r'(?=(AGG))|(?=(CGG))|(?=(TGG))|(?=(GGG))',
        }
        self.num_containing_fourmers_removed = 0
        self.num_containing_fivemers_removed = 0
        self.num_more_than_four_twomers_removed = 0


    def identify_guides_and_indicate_strand(
        self,
        pam: str,
        sequence: str,
        protospacer_length: int,
        context_toward_five_prime: int,
        context_toward_three_prime: int,
        include_repetitive: bool
        ) -> tuple[list[str], list[str], list[str], list[int]]:
        '''
        ## Args:
            * pam: The protospacer adjacent motif to look for, e.g., 'NGG'.
            * sequence: A nucleotide/DNA sequence to find substrings of guide RNA in.
            * protospacer_length: The length of the string toward 5-prime of the PAM.
            * context_toward_five_prime: The number of nucleotides to extract toward
                the 5-prime after the protospacer.
            * context_toward_three_prime: The number of nucleotides to extract toward
                the 3-prime after (and excluding) the PAM.
            * include_repetitive: Discards a guide if the protospace contains 2 or 
                fewer unique bases.

        ## Returns:
            A tuple of three lists:
            * The first list[str] is a list of the guides found in `sequence`.
            * The second list[str] is a list of the guides with their context around them.
            * The third list[str] is a list of 'F's and 'RC's indicating on
                which strand, forward or reverse complement, each respective guide resides.
            * The fourth list[int] shows the location of each guides in `sequence`.

        ## Example
            Input: find_guides_and_indicate_strand(
                pam='NGG',
                sequence='AAAAAACCAAAAAGGTTTTTT',
                protospacer_length=3,
                context_toward_five_prime=2,
                context_toward_three_prime=2,
            )
            Output:
            (['AAA', 'TTT'], ['CAAAAAGGTT', 'CTTTTTGGTT'], ['F', 'RC'], [12, 12])
        '''

        guides_list: list[str] = list()
        guides_context_list: list[str] = list()
        strands_list: list[str] = list()
        locations_list: list[int] = list()

        alphabet = ['A', 'C', 'G', 'T']

        # Store the reverse complement
        sequence_rev_comp = str(Seq(sequence).reverse_complement())

        pam_regex = r''
        if pam in self.pam_dict:
            pam_regex = self.pam_dict[pam]
        else:
            print(pam, 'not found in dictionary. Treating as literal.')
            pam_regex = r'(?=({PAM}))'.format(PAM=pam)
            

        def find_matches(seq: str, strand: str) -> None:        
            matches = re.finditer(pam_regex, seq)  # Find PAMs on seq
            pam_positions = [match.start() for match in matches]

            for position in pam_positions:
                if (position - protospacer_length >= 0):
                    guide = seq[position-protospacer_length:position]

                    if include_repetitive == False:
                        # Skip guides with bad nucleotides.
                        if any(c not in alphabet for c in guide):
                            continue

                        # Skip guides such as GGAGGAGGAGGAGGAGGAGG where GG is repeated
                        # 7 times or GA is repeated 6 times.
                        if max(count_kmers(guide, 2).values()) >= 5:
                            self.num_more_than_four_twomers_removed += 1
                            # print(guide, 'contains more than 5 2-mers')
                            continue

                        # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0784-0#Abs1:~:text=Repetitive%20bases%20are%20defined%20as%20any%20of%20the%20following
                        if any(substr in guide for substr in ['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT']):
                            self.num_containing_fivemers_removed += 1
                            # print(guide, 'contains 5-mers')
                            continue
                        
                        # https://www.nature.com/articles/s41467-019-12281-8#Abs1:~:text=but%20not%20significant).-,The%20contribution%20of,-repetitive%20nucleotides%20to
                        if any(substr in guide for substr in ['AAAA', 'CCCC', 'GGGG', 'TTTT']):
                            self.num_containing_fourmers_removed += 1
                            # print(guide, 'contains 4-mers')
                            continue

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
                    guides_context_list.append(guide_with_context)
                    strands_list.append(strand)
                    locations_list.append(position)

        find_matches(seq=sequence, strand='F')  # F means forward strand
        find_matches(seq=sequence_rev_comp, strand='RC')  # RC is the reverse complement strand

        return guides_list, guides_context_list, strands_list, locations_list


    def locate_guides_in_sequence(self, sequence, file_path, to_upper: bool = True):
        '''
        ## Args:
            * sequence: A guide RNA nucleotide sequence.
            * file_path: The file path to the DNA sequence (genome or CDS genes) fasta file. 
                Must be relative to the root directory of ALLEGRO.
            * to_upper: (Default: True) convert the fasta sequence to upper case or not?

        ## Returns:
            A list of tuples of 4 items list[tuple[str, str, int, int]]:
            * 1. (str) Chromosome name.
            * 2. (str) Strand. 'F' indicates the forward strand as read in the fasta file,
                'RC' indicates the reverse complement of the string in the same file.
            * 3. (int) Start position of `sequence` in `file_path` fasta.
            * 4. (int) End position of `sequence` in `file_path` fasta.
        '''

        sequence_rev_comp = str(Seq(sequence).reverse_complement())
        chrom_strand_start_end: list[tuple[str, str, int, int]] = list()

        def find_matches(seq: str, strand: str):
            for record in SeqIO.parse(open(file_path, 'r'), 'fasta'):
                chromosome = record.id
                dna = str(record.seq).upper() if to_upper else str(record.seq)

                for match in re.finditer(seq, dna):
                    chrom_strand_start_end.append((chromosome, strand, match.start(), match.end()))

        find_matches(seq=sequence, strand='F')  # F means forward strand
        find_matches(seq=sequence_rev_comp, strand='RC')  # RC is the reverse complement strand

        return chrom_strand_start_end