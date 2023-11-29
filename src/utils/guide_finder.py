# Functions imported by ALLEGRO. No need to run it manually.
# You can import this .py file separately, instantiate GuideFinder,
# and check for guides in your custom sequence by calling identify_guides_and_indicate_strand(...).
import re
import gc
import os
import pandas
from Bio import SeqIO
from Bio.Seq import Seq

from utils.shell_colors import bcolors


def count_kmers(sequence, k):
    kmers = dict()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers[kmer] = kmers.setdefault(kmer, 0) + 1

    return kmers


class GuideFinder:
    _self = None
    pam_dict = None
    exclusion_list = None

    def __new__(self):
        if self._self is None:
            self._self = super().__new__(self)
        return self._self
    
    
    def __init__(self) -> None:
        if self.pam_dict is None:
            self.pam_dict = {
                'NGG': r'(?=([ACTG]GG))',
            }

        if self.exclusion_list is None:
            # Whether there are guides to exclude
            if os.path.exists('data/input/_the_blacklist_.txt'):
                with open('data/input/_the_blacklist_.txt', 'r') as file:
                    self.exclusion_list = file.read().splitlines()
                    if len(self.exclusion_list) != 0:
                        print(f'{bcolors.BLUE}>{bcolors.RESET} Excluding the gRNA(s) in data/input/_the_blacklist_.csv')


    def identify_guides_and_indicate_strand(
        self,
        pam: str,
        sequence: str,
        protospacer_length: int,
        context_toward_five_prime: int,
        context_toward_three_prime: int,
        filter_repetitive: bool
        ) -> tuple[list[str], list[str], list[str], list[int]]:
        '''
        ## Args:
            * pam: The protospacer adjacent motif to look for, e.g., 'NGG'.
            * sequence: A nucleotide/DNA sequence to find substrings of guide RNA in.
            * protospacer_length: The length of the string toward 5-prime of the PAM.
            * context_toward_five_prime: The number of nucleotides to extract toward
                the 5-prime after the protospacer (to the left of the sequence).
            * context_toward_three_prime: The number of nucleotides to extract toward
                the 3-prime after (and excluding) the PAM (to the right side of the sequence).
            * filter_repetitive: Discards a guide if the protospace contains 2-mers
                repeated 5 or more times.

        ## Returns:
            A tuple of four lists:
            * The first list[str] is a list of the protospacers (WITHOUT PAM) found in `sequence`.
            * The second list[str] is a list of the protospacers with their context around them.
            This includes the PAM by default.
            * The third list[str] is a list of 'F's and 'RC's indicating on
                which strand, Forward or Reverse Complement, each respective guide resides.
            * The fourth list[int] shows the location of the start of the PAM of each guide in `sequence`.

        ## Example
            Input: find_guides_and_indicate_strand(
                pam='NGG',
                sequence='AGCGTACCCCCAGGTCTTGCAGG',
                protospacer_length=20,
                context_toward_five_prime=0,
                context_toward_three_prime=0,
            )
            Output:
            (['AGCGTACCCCCAGGTCTTGC'], ['AGCGTACCCCCAGGTCTTGCAGG'], ['F'], [20])
        '''

        guides_list: list[str] = list()
        guides_context_list: list[str] = list()
        strands_list: list[str] = list()
        locations_list: list[int] = list()

        pam_regex = r''
        if pam in self.pam_dict:
            pam_regex = self.pam_dict[pam]
        else:
            print('PAM', pam, 'not recognized. Treating as literal and searching for instances of', pam, '.')
            pam_regex = r'(?=({PAM}))'.format(PAM=pam)
            

        def find_matches(seq: str, strand: str) -> None:        
            matches = re.finditer(pam_regex, seq)  # Find PAMs on seq
            pam_positions = [match.start() for match in matches]

            for position in pam_positions:
                if (position - protospacer_length >= 0):
                    guide = seq[position-protospacer_length:position]

                    # Skip guides with non-standard nucleotides.
                    if any(c not in ['A', 'C', 'G', 'T'] for c in guide):
                        continue

                    if guide in self.exclusion_list:
                        continue

                    if filter_repetitive == True:
                        # Skip guides such as GGAGGAGGAGGAGGAGGAGG where GG is repeated
                        # 7 times or GA is repeated 6 times.
                        if max(count_kmers(guide, 2).values()) >= 5:
                            # self.num_more_than_four_twomers_removed += 1
                            # print(guide, 'contains more than 5 of the same 2-mer')
                            continue

                        # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0784-0#Abs1:~:text=Repetitive%20bases%20are%20defined%20as%20any%20of%20the%20following
                        if any(substr in guide for substr in ['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT']):
                            # self.num_containing_fivemers_removed += 1
                            # print(guide, 'contains 5-mers')
                            continue
                        
                        # https://www.nature.com/articles/s41467-019-12281-8#Abs1:~:text=but%20not%20significant).-,The%20contribution%20of,-repetitive%20nucleotides%20to
                        if any(substr in guide for substr in ['AAAA', 'CCCC', 'GGGG', 'TTTT']):
                            # self.num_containing_fourmers_removed += 1
                            # print(guide, 'contains 4-mers')
                            continue

                    guide_with_context = ''

                    # There is enough context on both sides of the PAM
                    if position-protospacer_length-context_toward_five_prime >= 0 and position+len(pam)+context_toward_three_prime < len(seq) + 1:
                        guide_with_context = seq[position-protospacer_length-context_toward_five_prime:position+len(pam)+context_toward_three_prime]
                    
                    else:
                        continue

                    guides_list.append(guide)
                    guides_context_list.append(guide_with_context)
                    strands_list.append(strand)
                    locations_list.append(position)

        sequence_rev_comp = str(Seq(sequence).reverse_complement())
        find_matches(seq=sequence, strand='+')  # F means forward strand
        find_matches(seq=sequence_rev_comp, strand='-')  # RC is the reverse complement strand

        return guides_list, guides_context_list, strands_list, locations_list


    def locate_guides_in_sequence(
        self,
        pam: str,
        sequence: str,
        file_path: str,
        to_upper: bool=True
        ) -> list[tuple[str, str, int, int, str, str]]:
        '''
        ## Args:
            * sequence: A guide RNA nucleotide sequence.
            * file_path: The file path to the DNA sequence (genome or CDS genes) fasta file. 
                Must be relative to the root directory of ALLEGRO.
            * to_upper: (Default: True) convert the fasta sequence to upper case?

        ## Returns:
            A list of tuples of 5 items: list[tuple[str, str, int, int, str]]:
            * 1. (str) Container (gene/chromosome) name.
            * 2. (str) Strand. 'F' indicates the forward strand as read in the fasta file,
                'RC' indicates the reverse complement of the string in the same file.
            * 3. (int) Start position of `sequence` in `file_path` fasta.
            * 4. (int) End position of `sequence` in `file_path` fasta.
            * 5. (str) Misc info about the guide and its location.
        '''

        chrom_strand_start_end_misc: list[tuple[str, str, int, int, str, str]] = list()


        def find_matches(seq: str, strand: str):
            for record in SeqIO.parse(open(file_path, 'r'), 'fasta'):

                dna: str = ''

                if strand == '+':
                    dna = str(record.seq).upper() if to_upper else str(record.seq)
                elif strand == '-':
                    dna = str(record.seq.reverse_complement()).upper() if to_upper else str(record.seq.reverse_complement())

                ortho_to = 'N/A'
                misc_list: list[str] = list()
                # match = re.search(r'\[Gene=(.*?)\]', record.description)
                # if match:
                #     misc_list.append(f'Gene: {match.group(1)}')

                match = re.search(r'\[orthologous_to_gene=(.*?)\]', record.description)
                if match:
                    ortho_to = match.group(1)

                # match = re.search(r'\[protein_id=(.*?)\]', record.description)
                # if match:
                #     misc_list.append(f'Protein: {match.group(1)}')

                misc_list.append(record.description.split(" ")[1])

                misc = ';'.join(misc_list)

                seq_and_pam = re.compile(seq + self.pam_dict[pam])
                for match in re.finditer(seq_and_pam, dna):
                    chrom_strand_start_end_misc.append((record.id, strand, match.start(), match.end(), ortho_to, misc))

        find_matches(seq=sequence, strand='+')  # F means forward strand
        find_matches(seq=sequence, strand='-')  # RC is the reverse complement strand

        return chrom_strand_start_end_misc


# UNUSED IN RELEASE CODE
class GuideFinderDebug:
    def __init__(self) -> None:
        self.pam_dict = {
            'NGG': r'(?=([ACTG]GG))',
        }

        self.removed_guide: list[str] = list()
        self.guide_species: list[str] = list()
        self.strand: list[str] = list()
        self.seq_ids: list[str] = list()
        self.location: list[str] = list()
        self.seq_length: list[int] = list()
        self.is_non_standard: list[int] = list()
        self.contains_fourmers: list[int] = list()
        self.contains_fivemers: list[int] = list()
        self.proportion_of_total_length: list[float] = list()
        self.contains_five_or_more_rep_twomers: list[int] = list()


    def identify_guides_and_indicate_strand(
        self,
        pam: str,
        sequence: str,
        protospacer_length: int,
        context_toward_five_prime: int,
        context_toward_three_prime: int,
        filter_repetitive: bool,
        sequence_id: str = '',
        species: str = '',
        ) -> tuple[list[str], list[str], list[str], list[int]]:

        guides_list: list[str] = list()
        guides_context_list: list[str] = list()
        strands_list: list[str] = list()
        locations_list: list[int] = list()

        pam_regex = r''
        if pam in self.pam_dict:
            pam_regex = self.pam_dict[pam]
        else:
            print('PAM', pam, 'not recognized. Treating as literal and searching for instances of', pam, '.')
            pam_regex = r'(?=({PAM}))'.format(PAM=pam)
            

        def find_matches(seq: str, strand: str) -> None:        
            matches = re.finditer(pam_regex, seq)  # Find PAMs on seq
            pam_positions = [match.start() for match in matches]

            for position in pam_positions:
                if (position - protospacer_length >= 0):
                    guide = seq[position-protospacer_length:position]

                    # Skip guides with non-standard nucleotides.
                    if any(c not in ['A', 'C', 'G', 'T'] for c in guide):
                        self.removed_guide.append(guide)
                        self.guide_species.append(species)
                        self.strand.append(strand)
                        self.seq_ids.append(sequence_id)
                        self.location.append(str(position))
                        self.seq_length.append(len(seq))
                        self.proportion_of_total_length.append(position / len(seq))
                        self.is_non_standard.append(1)
                        self.contains_five_or_more_rep_twomers.append(0)
                        self.contains_fivemers.append(0)
                        self.contains_fourmers.append(0)
                        continue

                    if filter_repetitive == True:
                        # Skip guides such as GGAGGAGGAGGAGGAGGAGG where GG is repeated
                        # 7 times or GA is repeated 6 times.
                        if max(count_kmers(guide, 2).values()) >= 5:
                            self.removed_guide.append(guide)
                            self.guide_species.append(species)
                            self.strand.append(strand)
                            self.seq_ids.append(sequence_id)
                            self.location.append(str(position))
                            self.seq_length.append(len(seq))
                            self.proportion_of_total_length.append(position / len(seq))
                            self.is_non_standard.append(0)
                            self.contains_five_or_more_rep_twomers.append(1)
                            self.contains_fivemers.append(0)
                            self.contains_fourmers.append(0)
                            continue

                        # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0784-0#Abs1:~:text=Repetitive%20bases%20are%20defined%20as%20any%20of%20the%20following
                        if any(substr in guide for substr in ['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT']):
                            self.removed_guide.append(guide)
                            self.guide_species.append(species)
                            self.strand.append(strand)
                            self.seq_ids.append(sequence_id)
                            self.location.append(str(position))
                            self.seq_length.append(len(seq))
                            self.proportion_of_total_length.append(position / len(seq))
                            self.is_non_standard.append(0)
                            self.contains_five_or_more_rep_twomers.append(0)
                            self.contains_fivemers.append(1)
                            self.contains_fourmers.append(0)
                            continue
                        
                        # https://www.nature.com/articles/s41467-019-12281-8#Abs1:~:text=but%20not%20significant).-,The%20contribution%20of,-repetitive%20nucleotides%20to
                        if any(substr in guide for substr in ['AAAA', 'CCCC', 'GGGG', 'TTTT']):
                            self.removed_guide.append(guide)
                            self.guide_species.append(species)
                            self.strand.append(strand)
                            self.seq_ids.append(sequence_id)
                            self.location.append(str(position))
                            self.seq_length.append(len(seq))
                            self.proportion_of_total_length.append(position / len(seq))
                            self.is_non_standard.append(0)
                            self.contains_five_or_more_rep_twomers.append(0)
                            self.contains_fivemers.append(0)
                            self.contains_fourmers.append(1)
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

        find_matches(seq=sequence, strand='+')  # F means forward strand

        sequence_rev_comp = str(Seq(sequence).reverse_complement())
        find_matches(seq=sequence_rev_comp, strand='-')  # RC is the reverse complement strand

        return guides_list, guides_context_list, strands_list, locations_list


    def locate_guides_in_sequence(self, sequence: str, file_path: str, to_upper: bool = True):
        '''
        ## Args:
            * sequence: A guide RNA nucleotide sequence.
            * file_path: The file path to the DNA sequence (genome or CDS genes) fasta file. 
                Must be relative to the root directory of ALLEGRO.
            * to_upper: (Default: True) convert the fasta sequence to upper case or not?

        ## Returns:
            A list of tuples of 4 items: list[tuple[str, str, int, int]]:
            * 1. (str) Chromosome/gene name.
            * 2. (str) Strand. 'F' indicates the forward strand as read in the fasta file,
                'RC' indicates the reverse complement of the string in the same file.
            * 3. (int) Start position of `sequence` in `file_path` fasta.
            * 4. (int) End position of `sequence` in `file_path` fasta.
        '''

        container_strand_start_end: list[tuple[str, str, int, int]] = list()

        def find_matches(seq: str, strand: str):
            for record in SeqIO.parse(open(file_path, 'r'), 'fasta'):

                dna = str(record.seq).upper() if to_upper else str(record.seq)

                for match in re.finditer(seq, dna):
                    container_strand_start_end.append((record.id, strand, match.start(), match.end()))

        find_matches(seq=sequence, strand='+')  # F means forward strand
        
        sequence_rev_comp = str(Seq(sequence).reverse_complement())
        find_matches(seq=sequence_rev_comp, strand='-')  # RC is the reverse complement strand

        return container_strand_start_end


    def write_removed_guides_to_dataframe(self) -> None:
        pandas.DataFrame({
            'removed_guide': self.removed_guide,
            'species': self.guide_species,
            'strand': self.strand,
            'seq_id': self.seq_ids,
            'location': self.location,
            'total_seq_length': self.seq_length,
            'proportion_of_total_location': self.proportion_of_total_length,
            'contains_nonstandard': self.is_non_standard,
            'contains_fourmers': self.contains_fourmers,
            'contains_fivemers': self.contains_fivemers,
            'contains_five_or_more_rep_twomers': self.contains_five_or_more_rep_twomers
        }).to_csv('removed_guides.csv', index=False)

        del self.strand
        del self.location
        del self.seq_ids
        del self.removed_guide
        del self.guide_species
        del self.is_non_standard
        del self.seq_length
        del self.proportion_of_total_length
        del self.contains_fivemers
        del self.contains_fourmers
        del self.contains_five_or_more_rep_twomers

        gc.collect()