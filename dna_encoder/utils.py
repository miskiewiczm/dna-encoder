"""
Helper functions for DNA encoder.

This module contains utility functions used by various components
of the DNA encoder, such as data conversions, pseudorandom generators,
mathematical functions, etc.
"""

import random
import hashlib
import math
from typing import Optional, List, Any, Iterator, Generator, Dict, Tuple
import logging
from collections import Counter

from .config import DNAEncoderConfig
from .exceptions import InputError
from dna_commons import normalize_dna_sequence

logger = logging.getLogger(__name__)



def string_to_bits(text: str, encoding: str = 'utf-8') -> str:
    """
    Converts string to binary representation.
    
    Each character is converted to 8-bit binary representation.
    
    Args:
        text: Tekst do konwersji
        encoding: Kodowanie tekstu (domyślnie UTF-8)
        
    Returns:
        String reprezentacji binarnej (np. "0100100001100101")
        
    Raises:
        InputError: Gdy tekst nie może być zakodowany
    """
    if not isinstance(text, str):
        raise InputError(f"Oczekiwano string, otrzymano {type(text)}", 
                        input_value=text, expected_type="str")
    
    try:
        # Konwertuj każdy znak na 8-bitową reprezentację
        bits = ''.join(format(byte, '08b') for byte in text.encode(encoding))
        logger.debug(f"Skonwertowano tekst (długość: {len(text)}) na {len(bits)} bitów")
        return bits
        
    except UnicodeEncodeError as e:
        raise InputError(f"Nie można zakodować tekstu w {encoding}: {str(e)}", 
                        input_value=text, expected_type=f"string encodable in {encoding}")


def bits_to_string(bits: str, encoding: str = 'utf-8') -> str:
    """
    Konwertuje reprezentację binarną z powrotem na string.
    
    Args:
        bits: String bitów (długość musi być wielokrotnością 8)
        encoding: Kodowanie do użycia przy dekodowaniu
        
    Returns:
        Zdekodowany tekst
        
    Raises:
        InputError: Gdy bity mają nieprawidłowy format lub długość
    """
    if not isinstance(bits, str):
        raise InputError(f"Oczekiwano string bitów, otrzymano {type(bits)}", 
                        input_value=bits, expected_type="str")
    
    if len(bits) % 8 != 0:
        raise InputError(f"Długość bitów ({len(bits)}) musi być wielokrotnością 8", 
                        input_value=bits)
    
    if not all(bit in '01' for bit in bits):
        raise InputError("String bitów can only contain characters '0' i '1'", 
                        input_value=bits)
    
    try:
        # Konwertuj każde 8 bitów na bajt
        bytes_list = []
        for i in range(0, len(bits), 8):
            byte_bits = bits[i:i+8]
            byte_value = int(byte_bits, 2)
            bytes_list.append(byte_value)
        
        # Dekoduj bajty na tekst
        text = bytes(bytes_list).decode(encoding)
        logger.debug(f"Zdekodowano {len(bits)} bitów na tekst (długość: {len(text)})")
        return text
        
    except (ValueError, UnicodeDecodeError) as e:
        raise InputError(f"Nie można zdekodować bitów: {str(e)}", 
                        input_value=bits[:50] + "..." if len(bits) > 50 else bits)


def validate_dna_sequence(sequence: str) -> bool:
    """
    Checks if sequence contains only valid DNA nucleotides.
    
    Args:
        sequence: Sekwencja do sprawdzenia
        
    Returns:
        True jeśli sekwencja jest prawidłowa
    """
    if not isinstance(sequence, str):
        return False
    
    valid_nucleotides = set('ATGCatgc')
    return all(nucleotide in valid_nucleotides for nucleotide in sequence)




def calculate_sequence_statistics(sequence: str) -> dict:
    """
    Calculates basic DNA sequence statistics.
    
    Args:
        sequence: Sekwencja DNA
        
    Returns:
        Słownik ze statystykami sekwencji
    """
    sequence = normalize_dna_sequence(sequence)
    length = len(sequence)
    
    if length == 0:
        return {
            'length': 0,
            'gc_content': 0.0,
            'nucleotide_counts': {'A': 0, 'T': 0, 'G': 0, 'C': 0},
            'nucleotide_frequencies': {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}
        }
    
    # Zlicz nukleotydy
    nucleotide_counts = {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'G': sequence.count('G'),
        'C': sequence.count('C')
    }
    
    # Oblicz częstotliwości
    nucleotide_frequencies = {
        nuc: count / length for nuc, count in nucleotide_counts.items()
    }
    
    # Oblicz zawartość GC
    gc_content = (nucleotide_counts['G'] + nucleotide_counts['C']) / length
    
    return {
        'length': length,
        'gc_content': gc_content,
        'nucleotide_counts': nucleotide_counts,
        'nucleotide_frequencies': nucleotide_frequencies
    }


def chunked_sequence(sequence: str, chunk_size: int) -> Generator[str, None, None]:
    """
    Dzieli sekwencję na chunki o zadanej wielkości.
    
    Args:
        sequence: Sekwencja do podzielenia
        chunk_size: Rozmiar każdego chunka
        
    Yields:
        Kolejne fragmenty sekwencji
        
    Raises:
        InputError: Gdy chunk_size jest nieprawidłowy
    """
    if chunk_size <= 0:
        raise InputError("Rozmiar chunka musi być dodatni", 
                        input_value=chunk_size, expected_type="positive integer")
    
    for i in range(0, len(sequence), chunk_size):
        yield sequence[i:i + chunk_size]


def find_longest_homopolymer(sequence: str) -> tuple:
    """
    Znajduje najdłuższy homopolimer w sekwencji DNA.

    Args:
        sequence: Sekwencja DNA

    Returns:
        Tuple (nucleotide, length, position)
    """
    """
    Znajduje najdłuższy homopolimer w sekwencji.
    
    Args:
        sequence: Sekwencja DNA
        
    Returns:
        Tuple (nucleotide, length, position) najdłuższego homopolimeru
    """
    if not sequence:
        return ('', 0, -1)
    
    sequence = normalize_dna_sequence(sequence)
    
    max_length = 1
    max_nucleotide = sequence[0]
    max_position = 0
    
    current_length = 1
    current_nucleotide = sequence[0]
    current_position = 0
    
    for i in range(1, len(sequence)):
        if sequence[i] == current_nucleotide:
            current_length += 1
        else:
            if current_length > max_length:
                max_length = current_length
                max_nucleotide = current_nucleotide
                max_position = current_position
            
            current_nucleotide = sequence[i]
            current_length = 1
            current_position = i
    
    # Sprawdź ostatni homopolimer
    if current_length > max_length:
        max_length = current_length
        max_nucleotide = current_nucleotide
        max_position = current_position
    
    return (max_nucleotide, max_length, max_position)



class ProgressTracker:
    """
    Class for tracking progress of long-running operations.
    """
    
    def __init__(self, total: int, description: str = "Progress"):
        """
        Inicjalizuje tracker postępu.
        
        Args:
            total: Całkowita liczba kroków
            description: Opis operacji
        """
        self.total = total
        self.current = 0
        self.description = description
        self.last_logged_percent = -1
    
    def update(self, increment: int = 1) -> None:
        """
        Aktualizuje postęp o zadaną liczbę kroków.
        
        Args:
            increment: Liczba kroków do dodania
        """
        self.current += increment
        
        # Loguj co 10%
        percent = int((self.current / self.total) * 100) if self.total > 0 else 100
        if percent >= self.last_logged_percent + 10:
            logger.info(f"{self.description}: {percent}% ({self.current}/{self.total})")
            self.last_logged_percent = percent
    
    def is_complete(self) -> bool:
        """
        Checks if operation is complete.
        
        Returns:
            True jeśli current >= total
        """
        return self.current >= self.total
    
    def get_percent(self) -> float:
        """
        Returns current completion percentage.
        
        Returns:
            Procent ukończenia (0.0-100.0)
        """
        return (self.current / self.total) * 100 if self.total > 0 else 100.0


def format_sequence_for_display(sequence: str, line_length: int = 60) -> str:
    """
    Formatuje sekwencję DNA do wyświetlania z podziałem na linie.
    
    Args:
        sequence: Sekwencja DNA
        line_length: Liczba nukleotydów na linię
        
    Returns:
        Sformatowana sekwencja z podziałem na linie
    """
    if not sequence:
        return ""
    
    lines = []
    for i in range(0, len(sequence), line_length):
        line = sequence[i:i + line_length]
        # Dodaj numery pozycji co 10 nukleotydów
        formatted_line = ""
        for j, nucleotide in enumerate(line):
            if j > 0 and j % 10 == 0:
                formatted_line += " "
            formatted_line += nucleotide
        
        # Dodaj numer pozycji na początku linii
        position = i + 1
        lines.append(f"{position:>6}: {formatted_line}")
    
    return "\n".join(lines)
