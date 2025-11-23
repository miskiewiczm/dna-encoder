"""
DNA decoder class.

This module contains the DNADecoder class that provides high-level API
for decoding DNA sequences back to original text and binary data.
"""

import logging
from typing import Optional

from .utils import (
    bits_to_string,
    normalize_dna_sequence,
)
from .exceptions import (
    DecodingError,
    InputError
)

logger = logging.getLogger(__name__)

# Bit to nucleotide mapping - shared with encoder
BIT_TO_NUCLEOTIDES = {
    '0': ['TA', 'TT', 'GC', 'CC', 'AC', 'AG', 'GT', 'CT'],
    '1': ['AT', 'AA', 'CG', 'GG', 'TC', 'TG', 'GA', 'CA']
}


class DNADecoder:
    """
    Class for decoding DNA sequences back to original data.

    Implements decoding algorithm that reverses the encoding process
    to recover original text and binary data from DNA sequences.
    """

    def __init__(self):
        """Initialize DNA decoder."""
        logger.info("Initialized DNADecoder")

    def decode_sequence(self, sequence: str, initial_sequence: str) -> str:
        """
        Decode DNA sequence back to original text.

        Args:
            sequence: DNA sequence to decode
            initial_sequence: Original starting sequence

        Returns:
            Decoded text string

        Raises:
            DecodingError: If decoding fails
        """
        logger.info(f"Starting DNA sequence decoding (length: {len(sequence)})")

        try:
            sequence, initial_sequence = self._normalize_decode_inputs(sequence, initial_sequence)
            encoded_part = self._extract_encoded_part(sequence, initial_sequence)
            bit_string = self._decode_pairs_to_bits(encoded_part)
            return self._convert_bits_to_text(bit_string)

        except DecodingError:
            raise
        except Exception as e:
            raise DecodingError(f"Unexpected decoding error: {e}", sequence=sequence)

    def _normalize_decode_inputs(self, sequence: str, initial_sequence: str) -> tuple[str, str]:
        """Normalize and validate decoding inputs."""
        sequence = normalize_dna_sequence(sequence)
        initial_sequence = normalize_dna_sequence(initial_sequence)

        if not sequence.startswith(initial_sequence):
            raise DecodingError(
                "Sequence does not start with expected initial sequence",
                sequence=sequence[:50] + "..." if len(sequence) > 50 else sequence,
                initial_sequence=initial_sequence
            )

        return sequence, initial_sequence

    def _extract_encoded_part(self, sequence: str, initial_sequence: str) -> str:
        """Extract and validate encoded part of sequence."""
        encoded_part = sequence[len(initial_sequence):]

        if len(encoded_part) % 2 != 0:
            logger.warning(f"Encoded part has odd length: {len(encoded_part)}")
            encoded_part = encoded_part[:-1]  # Trim last nucleotide

        return encoded_part

    def _decode_pairs_to_bits(self, encoded_part: str) -> str:
        """Decode nucleotide pairs to bit string."""
        bit_list = []
        unknown_pairs = []

        for i in range(0, len(encoded_part), 2):
            pair = encoded_part[i:i+2]
            bit_found = False

            for bit, nucleotide_options in BIT_TO_NUCLEOTIDES.items():
                if pair in nucleotide_options:
                    bit_list.append(bit)
                    bit_found = True
                    break

            if not bit_found:
                unknown_pairs.append(pair)
                logger.warning(f"Unknown nucleotide pair: {pair}")

        if unknown_pairs:
            logger.warning(f"Found {len(unknown_pairs)} unknown pairs during decoding")

        bit_string = ''.join(bit_list)
        logger.info(f"Decoded {len(bit_list)} bits")
        return bit_string

    def _convert_bits_to_text(self, bit_string: str) -> str:
        """Convert bit string to text."""
        if not bit_string:
            return ""

        try:
            return bits_to_string(bit_string)
        except Exception as e:
            raise DecodingError(f"Failed to convert bits to text: {e}", bit_string=bit_string)


# Export main classes
__all__ = ['DNADecoder', 'BIT_TO_NUCLEOTIDES']