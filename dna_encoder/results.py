"""
DNA encoding result structures.

This module contains data classes for representing results of DNA encoding operations,
including success/failure status, generated sequences, and detailed statistics.
"""

from dataclasses import dataclass
from typing import Optional, Dict, Any


@dataclass
class DNAResult:
    """
    Result of DNA encoding operation.

    This class encapsulates all information about a DNA encoding attempt,
    including whether it succeeded, the generated sequence, error details,
    and comprehensive statistics about the encoding process.

    Attributes:
        success: Whether the operation succeeded
        sequence: Generated DNA sequence (None if failed)
        error_message: Error message (None if successful)
        encoding_stats: Statistics from the encoding process
        validation_details: Details about sequence validation
    """
    success: bool
    sequence: Optional[str] = None
    error_message: Optional[str] = None
    encoding_stats: Optional[Dict[str, Any]] = None
    validation_details: Optional[Dict[str, Any]] = None

    def to_dict(self, include_details: bool = True) -> Dict[str, Any]:
        """
        Convert result to dictionary representation.

        Args:
            include_details: Whether to include detailed statistics

        Returns:
            Dictionary representation of the result
        """
        result = {
            'success': self.success,
            'sequence': self.sequence,
            'error_message': self.error_message,
        }

        if include_details:
            result['encoding_stats'] = self.encoding_stats
            result['validation_details'] = self.validation_details

        return result

    def get_sequence_length(self) -> int:
        """
        Get length of generated sequence.

        Returns:
            Length of sequence, 0 if no sequence generated
        """
        return len(self.sequence) if self.sequence else 0

    def get_encoding_attempts(self) -> int:
        """
        Get number of encoding attempts made.

        Returns:
            Number of attempts, 0 if no stats available
        """
        if self.encoding_stats and 'attempts' in self.encoding_stats:
            return self.encoding_stats['attempts']
        return 0

    def get_backtrack_count(self) -> int:
        """
        Get number of backtracking steps performed.

        Returns:
            Number of backtracks, 0 if no stats available
        """
        if self.encoding_stats and 'backtrack_count' in self.encoding_stats:
            return self.encoding_stats['backtrack_count']
        return 0

    def has_validation_details(self) -> bool:
        """
        Check if validation details are available.

        Returns:
            True if validation details are present
        """
        return self.validation_details is not None and len(self.validation_details) > 0

    def __str__(self) -> str:
        """String representation of the result."""
        if self.success:
            length = self.get_sequence_length()
            attempts = self.get_encoding_attempts()
            backtracks = self.get_backtrack_count()
            return f"DNAResult(success=True, length={length}, attempts={attempts}, backtracks={backtracks})"
        else:
            return f"DNAResult(success=False, error='{self.error_message}')"

    def __repr__(self) -> str:
        """Detailed representation of the result."""
        return self.__str__()


# Export the main class
__all__ = ['DNAResult']