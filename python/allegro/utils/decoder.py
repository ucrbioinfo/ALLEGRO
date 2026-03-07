nucleotide_map = {
    'A': 'AAA', 'B': 'AAT', 'C': 'AAC', 'D': 'AAG',
    'E': 'ATA', 'F': 'ATC', 'G': 'ATG', 'H': 'ATT',
    'I': 'ACA', 'J': 'ACC', 'K': 'ACG', 'L': 'ACT',
    'M': 'AGA', 'N': 'AGC', 'O': 'AGG', 'P': 'AGT',
    'Q': 'GAA', 'R': 'GAC', 'S': 'GAG', 'T': 'GAT',
    'U': 'GCA', 'V': 'GCC', 'W': 'GCG', 'X': 'GCT',
    'Y': 'GTA', 'Z': 'GTC', ' ': 'TTT', '!': 'TAA'
}

inverse_nucleotide_map = {v: k for k, v in nucleotide_map.items()}
def decode_sequence(nucleotide_sequence):
    triplets = [nucleotide_sequence[i:i+3] for i in range(0, len(nucleotide_sequence), 3)]
    return ''.join(inverse_nucleotide_map[triplet] for triplet in triplets if triplet in inverse_nucleotide_map)

encoded_message = ''
decoded_message = decode_sequence(encoded_message)

if encoded_message == '':
    print('The secret is in one of the figures.')
else:
    print(decoded_message)