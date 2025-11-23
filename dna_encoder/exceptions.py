"""
Exception definitions used by the DNA encoder.

This module contains a hierarchy of exceptions specific to the DNA encoder,
enabling precise handling of different types of errors.
"""

from typing import Optional, Dict, Any


class DNAEncodingError(Exception):
    """
    Base class for all DNA encoding errors.

    Attributes:
        message: Error message
        details: Additional error details in dictionary form
    """
    
    def __init__(self, message: str, details: Optional[Dict[str, Any]] = None):
        """
        Initialize DNAEncodingError exception.

        Args:
            message: Error description message
            details: Optional dictionary with additional error information
        """
        super().__init__(message)
        self.message = message
        self.details = details or {}

    def __str__(self) -> str:
        """Returns string representation of error with additional details."""
        if self.details:
            details_str = ", ".join(f"{k}={v}" for k, v in self.details.items())
            return f"{self.message} (details: {details_str})"
        return self.message


class ConfigurationError(DNAEncodingError):
    """
    Exception raised for configuration errors.

    Used when configuration parameters are invalid or inconsistent.
    """
    pass


class ValidationError(DNAEncodingError):
    """
    Exception raised for DNA sequence validation errors.

    Used when a sequence does not meet required quality criteria.
    """

    def __init__(self, message: str, sequence: Optional[str] = None,
                 failed_checks: Optional[Dict[str, Any]] = None):
        """
        Initialize ValidationError exception.

        Args:
            message: Validation error description message
            sequence: Sequence that failed validation
            failed_checks: Dictionary with failed check results
        """
        details = {}
        if sequence:
            details['sequence'] = sequence
        if failed_checks:
            details.update(failed_checks)
        
        super().__init__(message, details)
        self.sequence = sequence
        self.failed_checks = failed_checks or {}


class BacktrackingError(DNAEncodingError):
    """
    Exception raised when backtracking algorithm cannot find a solution.

    Used when all backtracking possibilities have been exhausted
    without finding a valid sequence.
    """

    def __init__(self, message: str, attempts: Optional[int] = None,
                 last_position: Optional[int] = None):
        """
        Initialize BacktrackingError exception.

        Args:
            message: Error description message
            attempts: Number of attempts made
            last_position: Last position where algorithm stopped
        """
        details = {}
        if attempts is not None:
            details['attempts'] = attempts
        if last_position is not None:
            details['last_position'] = last_position
        
        super().__init__(message, details)
        self.attempts = attempts
        self.last_position = last_position


class EncodingError(DNAEncodingError):
    """
    Exception raised for errors in data to DNA encoding process.

    Used when input data (text/bits) cannot be encoded into DNA sequence.
    """

    def __init__(self, message: str, input_data: Optional[str] = None,
                 bit_position: Optional[int] = None):
        """
        Initialize EncodingError exception.

        Args:
            message: Encoding error description message
            input_data: Input data that could not be encoded
            bit_position: Bit position where error occurred
        """
        details = {}
        if input_data:
            details['input_data'] = input_data
        if bit_position is not None:
            details['bit_position'] = bit_position
        
        super().__init__(message, details)
        self.input_data = input_data
        self.bit_position = bit_position


class DecodingError(DNAEncodingError):
    """
    Exception raised for errors in DNA to data decoding process.

    Used when DNA sequence cannot be properly decoded.
    """

    def __init__(self, message: str, dna_sequence: Optional[str] = None,
                 invalid_pair: Optional[str] = None):
        """
        Initialize DecodingError exception.

        Args:
            message: Decoding error description message
            dna_sequence: DNA sequence that cannot be decoded
            invalid_pair: Nucleotide pair that cannot be decoded
        """
        details = {}
        if dna_sequence:
            details['dna_sequence'] = dna_sequence
        if invalid_pair:
            details['invalid_pair'] = invalid_pair
        
        super().__init__(message, details)
        self.dna_sequence = dna_sequence
        self.invalid_pair = invalid_pair


class Primer3Error(DNAEncodingError):
    """
    Exception raised for errors related to primer3 library.

    Used when primer3 functions return errors or unexpected results.
    """

    def __init__(self, message: str, primer3_function: Optional[str] = None,
                 sequence: Optional[str] = None):
        """
        Initialize Primer3Error exception.

        Args:
            message: Primer3 error description message
            primer3_function: Name of primer3 function that caused the error
            sequence: Sequence being processed when error occurred
        """
        details = {}
        if primer3_function:
            details['primer3_function'] = primer3_function
        if sequence:
            details['sequence'] = sequence[:50] + ('...' if len(sequence) > 50 else '')
        
        super().__init__(message, details)
        self.primer3_function = primer3_function
        self.sequence = sequence


class InputError(DNAEncodingError):
    """
    Exception raised for input data errors.

    Used when input data has invalid format or content.
    """

    def __init__(self, message: str, input_value: Optional[Any] = None,
                 expected_type: Optional[str] = None):
        """
        Initialize InputError exception.

        Args:
            message: Error description message
            input_value: Invalid input value
            expected_type: Expected data type
        """
        details = {}
        if input_value is not None:
            details['input_value'] = str(input_value)
        if expected_type:
            details['expected_type'] = expected_type
        
        super().__init__(message, details)
        self.input_value = input_value
        self.expected_type = expected_type