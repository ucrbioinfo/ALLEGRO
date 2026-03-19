from itertools import product

iupac_dict: dict[str, list[str]] = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'],
}

def expand_iupac_sequence(seq: str) -> list[str]:
    seq = seq.strip().upper()

    for c in seq:
        if c not in iupac_dict:
            raise ValueError(f'{bcolors.RED}>{bcolors.RESET} Dev error. Invalid IUPAC code "{c}" in sequence "{seq}"')

    return [''.join(chars) for chars in product(*(iupac_dict[c] for c in seq))]

def expand_pam(pam: str) -> list[str]:
    return expand_iupac_sequence(pam)

def build_pam_regex(pam: str) -> str:
    pattern = ''

    for c in pam.upper():
        if c not in iupac_dict:
            raise ValueError(f'{bcolors.RED}>{bcolors.RESET} Dev error. Unknown IUPAC code: {c}')

        bases = ''.join(iupac_dict[c])
        pattern += f'[{bases}]'

    return rf'(?=({pattern}))'