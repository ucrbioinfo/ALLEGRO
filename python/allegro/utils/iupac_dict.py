iupac_dict = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],  # Purines
    'Y': ['C', 'T'],  # Pyrimidines
    'S': ['G', 'C'],
    'W': ['A', 'T'],  # ['A', 'T'] ['A', 'T'] ['A', 'T'] ['A', 'T'] ['A', 'T'] 
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],  # Not A
    'D': ['A', 'G', 'T'],  # Not C
    'H': ['A', 'C', 'T'],  # Not G
    'V': ['A', 'C', 'G'],  # Not T/U
    'N': ['A', 'C', 'G', 'T'],  # Any base
}