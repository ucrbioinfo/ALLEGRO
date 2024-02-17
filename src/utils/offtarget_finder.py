import os
import sys
import subprocess
import pandas as pd

from utils.shell_colors import bcolors


def reverse_complement(string):
    s = ''
    for c in string:
        if c == 'A': s += 'T'
        elif c == 'C': s += 'G'
        elif c == 'G': s += 'C'
        elif c == 'T': s += 'A'
    return s[::-1]


class OfftargetFinder:
    _self = None
    _base_path = None
    _cache_path = None

    # Singleton
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


    def run_bowtie_build(self, species_name: str, path_to_background_fasta: str, background_source: str) -> str:
        one_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_{background_source}_idx.1.ebwt')
        two_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_{background_source}_idx.2.ebwt')
        three_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_{background_source}_idx.3.ebwt')
        four_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_{background_source}_idx.4.ebwt')
        rev_one_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_{background_source}_idx.rev.1.ebwt')
        rev_two_exists = os.path.exists(f'{self._base_path}/{self._cache_path}/index/{species_name}_{background_source}_idx.rev.2.ebwt')

        if one_exists and two_exists and three_exists and four_exists and rev_one_exists and rev_two_exists:
            return
        else:
            # print(f'{bcolors.BLUE}>{bcolors.RESET} Creating Bowtie index for {species_name}')
            bowtie_build_command = ['bowtie-build', '-q', f'{self._base_path}/{path_to_background_fasta}', f'{self._base_path}/{self._cache_path}/index/{species_name}_{background_source}_idx']
            process = subprocess.Popen(bowtie_build_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            _, stderr = process.communicate()
            
            if stderr:
                    print(f'{bcolors.RED}>{bcolors.RESET} offtarget_finder.py: bowtie-build encountered an error: {stderr.decode()}')
                    print(f'{bcolors.RED}> Exiting. {bcolors.RESET}')
                    sys.exit(1)
                    

    def run_bowtie_against_other(self, this_species_name: str, that_species_name: str, background_source: str, num_of_mismatches: int) -> None:
        print(f'{bcolors.BLUE}>{bcolors.RESET} Running Bowtie with gRNA reads from {this_species_name} against {that_species_name}')

        for v in [num_of_mismatches]:
            bowtie_command = ['bowtie', '-v', str(v), '-a', '--quiet', f'{self._base_path}/{self._cache_path}/index/{that_species_name}_{background_source}_idx', f'{self._base_path}/{self._cache_path}/reads/{this_species_name}_reads.fq', f'{self._base_path}/{self._cache_path}/alignments/{this_species_name}_against_{that_species_name}_{v}mm_alignment.sam']
            process = subprocess.Popen(bowtie_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            _, stderr = process.communicate()
            
            if stderr:
                    print(f'{bcolors.RED}>{bcolors.RESET} offtarget_finder.py: bowtie encountered an error: {stderr.decode()}')
                    print('Exiting.')
                    sys.exit(1)

            # os.system(f'bowtie -v {v} -a --quiet {self._base_path}/{self._cache_path}/index/{that_species_name}_{background_source}_idx {self._base_path}/{self._cache_path}/reads/{this_species_name}_reads.fq {self._base_path}/{self._cache_path}/alignments/{this_species_name}_against_{that_species_name}_{v}mm_alignment.sam')

            df_mm_genomic = pd.read_csv(f'{self._base_path}/{self._cache_path}/alignments/{this_species_name}_against_{that_species_name}_{v}mm_alignment.sam', sep='\t',
                    names=['query_name', 'strand', 'reference_name', 'start_position', 'aligned_seq', 'mapping_quality', 'idk', 'mismatch'])

            if len(df_mm_genomic) == 0:
                continue
            
            guide_w_pam = list()

            for _, row in df_mm_genomic.iterrows():
                if row['strand'] == '-':
                    guide_w_pam.append(reverse_complement(row['aligned_seq']))
                else:
                    guide_w_pam.append(row['aligned_seq'])

            df_mm_genomic['pam'] = [s[-3:] for s in guide_w_pam]
            df_mm_genomic['sequence'] = [s[:-3] for s in guide_w_pam]
            
            # Remove mismatches occuring at the N base of the NGG PAM
            df_mm_genomic['mismatch'] = df_mm_genomic['mismatch'].astype(str)
            df_mm_genomic['mismatch'].replace('nan', 'N/A', inplace=True)
            df_mm_genomic = df_mm_genomic[~df_mm_genomic['mismatch'].str.contains('20:')]

            # df_mm_genomic['num_mutations'] = df_mm_genomic['num_mutations'].astype(int)
            df_mm_genomic.drop(columns=['idk', 'mapping_quality'], inplace=True)

            os.remove(f'{self._base_path}/{self._cache_path}/alignments/{this_species_name}_against_{that_species_name}_{v}mm_alignment.sam')
            
            return df_mm_genomic.reset_index(drop=True)
        