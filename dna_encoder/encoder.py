"""
Main DNA encoder class.

This module contains the main DNAEncoder class that provides high-level API
for encoding text and binary data into DNA sequences with biochemical quality control.
"""

import logging
from typing import Optional, Dict, Any

from .config import DNAEncoderConfig
from dna_commons import DNAValidator, Primer3Adapter
from .results import DNAResult
from .backtracking_engine import BacktrackingEngine
from dna_commons import (
    DeterministicRandom,
    normalize_dna_sequence,
    SequenceAnalyzer
)
from .utils import string_to_bits
from .exceptions import (
    DNAEncodingError,
    BacktrackingError,
    EncodingError,
    InputError
)

logger = logging.getLogger(__name__)

# Bit to nucleotide mapping - shared with backtracking engine
BIT_TO_NUCLEOTIDES = {
    '0': ['TA', 'TT', 'GC', 'CC', 'AC', 'AG', 'GT', 'CT'],
    '1': ['AT', 'AA', 'CG', 'GG', 'TC', 'TG', 'GA', 'CA']
}

# Helper functions for primer3 compatibility
def calculate_tm(sequence: str, config: Optional[DNAEncoderConfig] = None) -> float:
    """Calculate melting temperature of sequence."""
    from dna_commons import Primer3Adapter
    adapter = Primer3Adapter()
    return adapter.calculate_tm(sequence)

def calculate_hairpin_tm(sequence: str) -> float:
    """Calculate Tm for potential hairpin. Limits sequence to 60 bp."""
    from dna_commons import Primer3Adapter
    adapter = Primer3Adapter()
    return adapter.calculate_hairpin_tm(sequence)

def calculate_homodimer_tm(sequence: str) -> float:
    """Calculate Tm for potential homodimer. Limits sequence to 60 bp."""
    from dna_commons import Primer3Adapter
    adapter = Primer3Adapter()
    return adapter.calculate_homodimer_tm(sequence)

def check_3_prime_stability(sequence: str) -> bool:
    """Check if at most 3 GC nucleotides in last 5 nt."""
    last_5 = sequence[-5:]
    return (last_5.count('G') + last_5.count('C')) <= 3

def check_homopolymer_runs(sequence: str, max_length: int = 4) -> bool:
    """Check if sequence contains runs of repeating nucleotides."""
    return all(base * max_length not in sequence for base in "ATGC")

def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content."""
    return (sequence.count('G') + sequence.count('C')) / len(sequence)

def check_repeats(sequence: str) -> bool:
    """Check if sequence contains repeating dinucleotides."""
    return not any(sequence[i:i+2] == sequence[i+2:i+4] for i in range(len(sequence) - 3))


class DNAEncoder:
    """
    Main class for encoding data into DNA sequences.

    Implements backtracking algorithm with biochemical quality control
    compatible with original dna_encoder_backtrack.py
    """

    def __init__(self, config: Optional[DNAEncoderConfig] = None):
        """
        Initialize DNA encoder.

        Args:
            config: Encoder configuration (None = default)
        """
        self.config = config or DNAEncoderConfig()

        # Initialize core components with dna_commons adapter
        from dna_commons import ValidationRules, ThermodynamicParams

        # Convert DNAEncoderConfig to ValidationRules
        # Use validation_rules from config (loaded from profile JSON)
        val_rules = self.config.validation_rules
        rules = ValidationRules(
            gc_content=val_rules.get('gc_content', True),
            melting_temperature=val_rules.get('melting_temperature', True),
            homopolymer_runs=val_rules.get('homopolymer_runs', True),
            dinucleotide_repeats=val_rules.get('dinucleotide_repeats', True),
            three_prime_stability=val_rules.get('three_prime_stability', True),
            hairpin_structures=val_rules.get('hairpin_structures', True),
            homodimer_structures=val_rules.get('homodimer_structures', True),
            min_gc=self.config.min_gc,
            max_gc=self.config.max_gc,
            min_tm=self.config.min_tm,
            max_tm=self.config.max_tm,
            max_hairpin_tm=self.config.max_hairpin_tm,
            max_homodimer_tm=self.config.max_homodimer_tm,
            max_homopolymer_length=self.config.max_homopolymer_length,
            max_dinucleotide_repeats=self.config.max_dinucleotide_repeats
        )

        # Create thermodynamic parameters - load from default config
        thermoparams = ThermodynamicParams.load_default()

        self.validator = DNAValidator(rules, thermoparams)
        self.analyzer = SequenceAnalyzer()
        self.primer3_adapter = Primer3Adapter()

        # Initialize backtracking engine
        self.backtracking_engine = BacktrackingEngine(
            self.config, self.validator, self.analyzer
        )

        logger.info("Initialized DNAEncoder with modular architecture")

    def encode_text(self, initial_sequence: str, text: str, max_attempts: int = 5) -> DNAResult:
        """
        Encode text into DNA sequence.

        Args:
            initial_sequence: Starting DNA sequence
            text: Text to encode
            max_attempts: Maximum number of encoding attempts

        Returns:
            DNAResult with encoding outcome
        """
        logger.info(f"Starting text encoding: '{text}'")

        # Normalize input sequence and convert text to bits
        initial_sequence = normalize_dna_sequence(initial_sequence)
        bit_sequence = string_to_bits(text)

        return self.encode_bits(initial_sequence, bit_sequence, max_attempts)

    def encode_bits(self, initial_sequence: str, bit_sequence: str, max_attempts: int = 5) -> DNAResult:
        """
        Encode bit sequence into DNA.

        Args:
            initial_sequence: Starting DNA sequence
            bit_sequence: Bit string (e.g., "01001100")
            max_attempts: Maximum number of encoding attempts

        Returns:
            DNAResult with encoding outcome
        """
        logger.info(f"Starting encoding of {len(bit_sequence)} bits to DNA")
        initial_sequence = normalize_dna_sequence(initial_sequence)
        last_stats: Dict[str, Any] = {}

        for attempt in range(1, max_attempts + 1):
            logger.info(f"Encoding attempt {attempt}/{max_attempts}")

            try:
                result, attempt_stats = self.backtracking_engine.encode_bits_single_attempt(
                    initial_sequence, bit_sequence
                )
                last_stats = attempt_stats

                # Check if all bits were encoded (each bit = 2 nucleotides)
                expected_length = len(initial_sequence) + (len(bit_sequence) * 2)
                if result and len(result) >= expected_length:
                    return self._create_success_result(result, attempt, attempt_stats)
                else:
                    if result:
                        logger.info(f"Attempt {attempt} completed: FAILURE (partial: {len(result)}/{expected_length} nt)")
                    else:
                        logger.info(f"Attempt {attempt} completed: FAILURE")

            except Exception as e:
                logger.error(f"Error in attempt {attempt}: {e}")
                last_stats = {'error': str(e), 'attempt': attempt}

        return self._create_failure_result(max_attempts, last_stats)

    def _create_success_result(self, result: str, attempt: int, attempt_stats: Dict[str, Any]) -> DNAResult:
        """Create successful encoding result with validation."""
        logger.info(f"Success! Encoded in attempt {attempt}")

        encoding_stats = {
            'attempts': attempt,
            'successful_attempt': attempt,
            **attempt_stats
        }

        validation_details = self._validate_final_sequence(result)

        return DNAResult(
            success=True,
            sequence=result,
            encoding_stats=encoding_stats,
            validation_details=validation_details
        )

    def _validate_final_sequence(self, sequence: str) -> Optional[Dict[str, Any]]:
        """Validate final sequence and return validation details."""
        try:
            validation_result = self.validator.validate_sequence(sequence)
            return {
                'final_validation': validation_result,
                'is_valid': validation_result.is_valid
            }
        except Exception as e:
            logger.warning(f"Final validation failed: {e}")
            return {'final_validation_error': str(e)}

    def _create_failure_result(self, max_attempts: int, last_stats: Dict[str, Any]) -> DNAResult:
        """Create failure result when all attempts exhausted."""
        logger.error(f"Failed to encode in {max_attempts} attempts")

        error_stats = {
            'attempts': max_attempts,
            'successful_attempt': None,
            **last_stats
        }

        return DNAResult(
            success=False,
            error_message=f"Encoding failed after {max_attempts} attempts",
            encoding_stats=error_stats
        )


    def get_config(self) -> DNAEncoderConfig:
        """Get current encoder configuration."""
        return self.config

    def set_config(self, config: DNAEncoderConfig) -> None:
        """
        Update encoder configuration.

        Args:
            config: New configuration

        Note:
            This requires reinitializing internal components.
        """
        self.config = config

        # Convert DNAEncoderConfig to ValidationRules for dna_commons
        from dna_commons import ValidationRules, ThermodynamicParams

        # Use validation_rules from config (loaded from profile JSON)
        val_rules = self.config.validation_rules
        rules = ValidationRules(
            gc_content=val_rules.get('gc_content', True),
            melting_temperature=val_rules.get('melting_temperature', True),
            homopolymer_runs=val_rules.get('homopolymer_runs', True),
            dinucleotide_repeats=val_rules.get('dinucleotide_repeats', True),
            three_prime_stability=val_rules.get('three_prime_stability', True),
            hairpin_structures=val_rules.get('hairpin_structures', True),
            homodimer_structures=val_rules.get('homodimer_structures', True),
            min_gc=self.config.min_gc,
            max_gc=self.config.max_gc,
            min_tm=self.config.min_tm,
            max_tm=self.config.max_tm,
            max_hairpin_tm=self.config.max_hairpin_tm,
            max_homodimer_tm=self.config.max_homodimer_tm,
            max_homopolymer_length=self.config.max_homopolymer_length,
            max_dinucleotide_repeats=self.config.max_dinucleotide_repeats
        )

        # Create thermodynamic parameters - load from default config
        thermoparams = ThermodynamicParams.load_default()

        self.validator = DNAValidator(rules, thermoparams)
        self.primer3_adapter = Primer3Adapter()
        self.backtracking_engine = BacktrackingEngine(
            self.config, self.validator, self.analyzer
        )
        logger.info("Updated DNAEncoder configuration")

    def validate_sequence(self, sequence: str) -> bool:
        """
        Validate DNA sequence against current configuration.

        Args:
            sequence: DNA sequence to validate

        Returns:
            True if sequence meets quality criteria
        """
        try:
            metrics = self.validator.validate_sequence(sequence)
            return metrics.is_valid
        except Exception as e:
            logger.error(f"Validation error: {e}")
            return False


# Export main classes and functions
__all__ = [
    'DNAEncoder',
    'BIT_TO_NUCLEOTIDES',
    'calculate_tm',
    'calculate_hairpin_tm',
    'calculate_homodimer_tm',
    'check_3_prime_stability',
    'check_homopolymer_runs',
    'calculate_gc_content',
    'check_repeats'
]