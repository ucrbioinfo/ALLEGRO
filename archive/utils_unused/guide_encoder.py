

class DNAEncoderDecoder:
    def __init__(self) -> None:
        self.encoding_map = {'A': 0, 'G': 1, 'C': 2, 'T': 3, 'N': 4}
        self.alphabet = ['A', 'G', 'C', 'T', 'N']

        self.bits_to_shift = (len(self.alphabet) - 1).bit_length()

        # For AGCT (alphabet of length 4) this is 4
        # For AGCTN this is 8
        # Cast to float to save some bytes... thanks python.
        self.integer_multiplier = float(2 ** self.bits_to_shift)

        # For AGCT (length 4) we need 2 bits to right shift -- this is 0b11 or 3 in decimal
        # For AGCTN (length 5) we need 3 bits to right shift -- this is 0b111 or 7 in decimal
        self.int_bits_to_shift = (len(self.alphabet) - 1).bit_length() | (len(self.alphabet) - 1)


    def encode(self, seq: str) -> float:
        code = 0.0
        for nc in seq:
            code *= self.integer_multiplier  # 4 for AGCT, 8 for AGCTN
            code += self.encoding_map[nc]
        return code


    def decode(self, encoding: float, length: int) -> str:
        decoded_str = ''
        code = int(encoding)
        
        for _ in range(length):
            index = code & self.int_bits_to_shift  # 3 for AGCT, 7 for AGCTN
            code >>= self.bits_to_shift  # 2 for AGCT, 3 for AGCTN
            decoded_str = self.alphabet[index] + decoded_str
            
        return decoded_str