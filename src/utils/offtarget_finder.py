import os
import re
import pandas as pd

from utils.shell_colors import bcolors

class OfftargetFinder:
    _self = None
    _base_path = None
    _cache_path = None

    def __new__(self):
        if self._self is None:
            self._self = super().__new__(self)
        return self._self
    
    def __init__(self) -> None:
        if not os.path.exists('data/cache/bowtie/index'):
            os.makedirs('data/cache/bowtie/index')
        if not os.path.exists('data/cache/bowtie/reads'):
            os.makedirs('data/cache/bowtie/reads')
        if not os.path.exists('data/cache/bowtie/alignments'):
            os.makedirs('data/cache/bowtie/alignments')
        
        if self._base_path is None:
            self._base_path = os.getcwd()
        
        if self._cache_path is None:
            self._cache_path = 'data/cache/bowtie'

    
    def write_guides_as_reads(self, species_name: str, guide_seqs_list: list[str]) -> None:
        guides_w_pam = dict()

        for seq in guide_seqs_list:
            with_pam = [seq + pam for pam in ['AGG', 'CGG', 'TGG', 'GGG']]

            for wp in with_pam:
                guides_w_pam[wp] = seq

        reads_path = f'{self._cache_path}/reads/{species_name}_reads.fq'
        with open(reads_path, 'w') as f:
            for idx, guide in enumerate(guides_w_pam.keys()):
                f.write(f'@READ_{idx+1}\n{guide}\n+\nIIIIIIIIIIIIIIIIIIIIIII\n')


    def run_bowtie_build(self, species_name: str, path_to_background_fasta: str) -> str:
        one_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_idx.1.ebwt')
        two_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_idx.2.ebwt')
        three_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_idx.3.ebwt')
        four_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_idx.4.ebwt')
        rev_one_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_idx.rev.1.ebwt')
        rev_two_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_idx.rev.2.ebwt')

        if one_exists and two_exists and three_exists and four_exists and rev_one_exists and rev_two_exists:
            return
        else:
            print(f'{bcolors.BLUE}>{bcolors.RESET} Creating Bowtie index for {species_name}')
            os.system(f'bowtie-build -q {self._base_path}/{path_to_background_fasta} {self._base_path}/{self._cache_path}/index/{species_name}_idx')

    def run_bowtie_against_other(self, this_species_name: str, that_species_name: str, guide_seq_list: list[str], seed_region_is_n_from_pam: int, num_of_mismatches: int) -> None:
        print(f'{bcolors.BLUE}>{bcolors.RESET} Running Bowtie with gRNA reads from {this_species_name} against {that_species_name}')

        def extract_complete_digits(s):
            return [int(num) for num in re.findall(r'\d+', s)]

        def find_indices_of_string(target, string_list):
            return [index for index, string in enumerate(string_list) if string == target]
        
        def determine_ot_in_nonseed(strand, digits):
            if strand == '+':
                return any(digit >= 23 - seed_region_is_n_from_pam for digit in digits)
            elif strand == '-':
                return any(23 - digit - 1 >= 23 - seed_region_is_n_from_pam for digit in digits)
            else:
                return False


        genome_offtargets = [0] * len(guide_seq_list)
        mm_allowed = [0] if (num_of_mismatches == 0) else [0, num_of_mismatches]

        for v in mm_allowed:
            os.system(f'bowtie -v {v} -a --quiet {self._base_path}/{self._cache_path}/index/{that_species_name}_idx {self._base_path}/{self._cache_path}/reads/{this_species_name}_reads.fq {self._base_path}/{self._cache_path}/alignments/{this_species_name}_against_{that_species_name}_{v}mm_alignment.sam')

            df_mm_genomic = pd.read_csv(f'{self._base_path}/{self._cache_path}/alignments/{this_species_name}_against_{that_species_name}_{v}mm_alignment.sam', sep='\t',
                    names=['query_name', 'strand', 'reference_name', 'mapping_position', 'aligned_seq', 'mapping_quality', 'idk', 'mismatch'])

            if len(df_mm_genomic) == 0:
                continue
            
            guide_w_pam = list()

            for idx, row in df_mm_genomic.iterrows():
                if row['strand'] == '-':
                    guide_w_pam.append(rc(row['aligned_seq']))
                else:
                    guide_w_pam.append(row['aligned_seq'])

            df_mm_genomic['guide_w_pam'] = guide_w_pam
            df_mm_genomic['guide'] = df_mm_genomic['guide_w_pam'].str[:-3]

            # Create a new column called mismatch_position and extract the mismatch position digits from
            # the mismatch column. Then fill NaN values with 0 because NaN means no mismatch.
            # Then drop rows where the mismatch occurs after the 10th base.
            df_mm_genomic['mismatch'] = df_mm_genomic['mismatch'].astype(str)
            df_mm_genomic['mismatch'].replace('nan', 'Exact Match', inplace=True)
            df_mm_genomic = df_mm_genomic[~df_mm_genomic['mismatch'].str.contains('20:')]
            df_mm_genomic['num_mutations'] = df_mm_genomic['mismatch'].apply(lambda x: len(x.split(',')))
            df_mm_genomic.loc[df_mm_genomic['mismatch'] == 'Exact Match', 'num_mutations'] = 0

            df_mm_genomic['num_mutations'] = df_mm_genomic['num_mutations'].astype(int)
            df_mm_genomic.drop(columns=['idk', 'mapping_quality'], inplace=True)

            if num_of_mismatches == 0:
                
                return df_mm_genomic

                # duplicates = df_mm_genomic.duplicated(subset='guide', keep='first')

                # # To see if there are any duplicates
                # if duplicates.any():
                #     for _, row in df_mm_genomic[duplicates].iterrows():
                #         idxs = find_indices_of_string(row['guide'], guide_seq_list)
                #         for idx in idxs:
                #             genome_offtargets[idx] += 1

            if num_of_mismatches > 0:
                df_mm_genomic = df_mm_genomic[df_mm_genomic.duplicated(subset='guide', keep=False)]

                # remove rows with nan in mismatch. these should be caught in the previous step
                df_mm_genomic = df_mm_genomic[df_mm_genomic['mismatch'] != 'Exact Match']

                # make a list out of the mismatch locations
                df_mm_genomic['digits'] = df_mm_genomic['mismatch'].apply(extract_complete_digits)

                # if positive strand and if any of the elements in the list are mm loc >= 23 - seed region, not an offtarget
                # if negative strand and if any of the elements in the list are 23 - mismatch location - 1 >= 23 - the seed region, not an offtarget
                df_mm_genomic['ot_in_seed'] = df_mm_genomic.apply(lambda row: determine_ot_in_nonseed(row['strand'], row['digits']), axis=1)

                # False are bad news. False in ot_in_seed means offtarget.
                if not all(df_mm_genomic['ot_in_seed']):
                    for _, row in df_mm_genomic.iterrows():
                        if row['ot_in_seed'] == False:
                            idxs = find_indices_of_string(row['guide'], guide_seq_list)
                            for idx in idxs:
                                genome_offtargets[idx] += 1

            # Clean up
            # os.remove(f'{self._base_path}/{self._cache_path}/alignments/{this_species_name}_against_{that_species_name}_{v}mm_alignment.sam')
        
        # return genome_offtargets


def rc(string):
    s = ''
    for c in string:
        if c == 'A': s += 'T'
        elif c == 'C': s += 'G'
        elif c == 'G': s += 'C'
        elif c == 'T': s += 'A'
    return s[::-1]
