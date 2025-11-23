"""
DNA Encoder - Refactored package for encoding data into DNA sequences.

This package enables encoding text and binary data into DNA sequences
using a backtracking algorithm and biochemical quality control.

The encoder is now fully compatible with the original dna_encoder_backtrack.py

Main components:
    DNAEncoder: Main class for encoding data into DNA
    DNADecoder: Class for decoding DNA sequences back to text
    DNAEncoderConfig: Configuration for encoding and validation parameters
    DNAValidator: DNA sequence quality validator

Example usage:
    from dna_encoder import DNAEncoder, DNAEncoderConfig

    # Use default configuration
    encoder = DNAEncoder()
    result = encoder.encode_text("ATGCATGC", "Hello World")

    # Or with custom configuration
    config = DNAEncoderConfig(min_gc=0.4, max_gc=0.6, window_size=15)
    encoder = DNAEncoder(config)
    result = encoder.encode_text("ATGCATGC", "Hello World")

    if result.success:
        print(f"DNA sequence: {result.sequence}")
        decoded = encoder.decode_sequence(result.sequence, "ATGCATGC")
        print(f"Decoded text: {decoded}")
"""

from .config import DNAEncoderConfig, DEFAULT_CONFIG
from .encoder import DNAEncoder
from .decoder import DNADecoder
from .results import DNAResult
from dna_commons import DNAValidator, QualityMetrics, Primer3Adapter, PRIMER3_AVAILABLE
from .backtracking_engine import BacktrackingEngine
from .exceptions import (
    DNAEncodingError,
    ConfigurationError,
    ValidationError,
    BacktrackingError,
    EncodingError,
    DecodingError,
    Primer3Error,
    InputError
)
from dna_commons import (
    DeterministicRandom,
    generate_seed_from_string,
    SequenceAnalyzer,
    normalize_dna_sequence
)
from .utils import (
    string_to_bits,
    bits_to_string,
    validate_dna_sequence,
    calculate_sequence_statistics,
    format_sequence_for_display
)

# Version information
__version__ = "0.1.0"
__author__ = "Marek Miskiewicz"
__description__ = "DNA data encoder with biochemical quality control"

# Main exported classes and functions
__all__ = [
    # Main classes
    'DNAEncoder',
    'DNADecoder',
    'DNAValidator',
    'QualityMetrics',
    'DNAEncoderConfig',
    'DNAResult',
    'DeterministicRandom',
    'SequenceAnalyzer',
    'Primer3Adapter',
    'BacktrackingEngine',

    # Configuration
    'DEFAULT_CONFIG',
    'PRIMER3_AVAILABLE',

    # Exceptions
    'DNAEncodingError',
    'ConfigurationError',
    'ValidationError',
    'BacktrackingError',
    'EncodingError',
    'DecodingError',
    'Primer3Error',
    'InputError',

    # Utility functions
    'string_to_bits',
    'bits_to_string',
    'generate_seed_from_string',
    'validate_dna_sequence',
    'normalize_dna_sequence',
    'calculate_sequence_statistics',
    'format_sequence_for_display',

    # Metadata
    '__version__',
    '__author__',
    '__description__'
]


def create_encoder(min_gc: float = 0.45, max_gc: float = 0.55,
                  min_tm: float = 55.0, max_tm: float = 65.0,
                  window_size: int = 20, deterministic: bool = True,
                  **kwargs) -> DNAEncoder:
    """
    Helper function to create an encoder with commonly used parameters.

    Args:
        min_gc: Minimum GC content (0.0-1.0)
        max_gc: Maximum GC content (0.0-1.0)
        min_tm: Minimum melting temperature in °C
        max_tm: Maximum melting temperature in °C
        window_size: Analysis window size
        deterministic: Whether to enable deterministic mode
        **kwargs: Additional parameters for DNAEncoderConfig

    Returns:
        Configured DNAEncoder instance

    Example:
        encoder = create_encoder(min_gc=0.4, max_gc=0.6, window_size=15)
        result = encoder.encode_text("ATGC", "Hello")
    """
    config = DNAEncoderConfig(
        min_gc=min_gc,
        max_gc=max_gc,
        min_tm=min_tm,
        max_tm=max_tm,
        window_size=window_size,
        enable_deterministic=deterministic,
        **kwargs
    )
    return DNAEncoder(config)


def quick_encode(text: str, initial_sequence: str = "ATGCATGCATGCATGC",
                **kwargs) -> DNAResult:
    """
    Quick function to encode text with default settings.

    Args:
        text: Text to encode
        initial_sequence: Initial sequence
        **kwargs: Additional parameters for create_encoder

    Returns:
        DNAResult with encoding results

    Example:
        result = quick_encode("Hello World")
        if result.success:
            print(f"DNA: {result.sequence}")
    """
    encoder = create_encoder(**kwargs)
    return encoder.encode_text(initial_sequence, text)


def quick_decode(dna_sequence: str, initial_sequence: str = "ATGCATGCATGCATGC") -> str:
    """
    Quick function to decode DNA sequence with default settings.

    Args:
        dna_sequence: DNA sequence to decode
        initial_sequence: Initial sequence

    Returns:
        Decoded text

    Example:
        text = quick_decode("ATGCATGCATGCATGCTATACGCG...")
        print(f"Text: {text}")
    """
    encoder = DNAEncoder()
    return encoder.decode_sequence(dna_sequence, initial_sequence)


# Logging configuration for the entire package
import logging

def setup_logging(level: str = "INFO", format_string: str = None):
    """
    Configure logging for the dna_encoder package.

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        format_string: Custom log format (None = default)
    """
    if format_string is None:
        format_string = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format=format_string
    )

    # Set level for all package loggers
    package_logger = logging.getLogger(__name__)
    package_logger.setLevel(getattr(logging, level.upper()))


# Automatic logging configuration on import
setup_logging()
