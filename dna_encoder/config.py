"""
DNA encoder configuration.

This module contains all configurations and parameters used by the DNA encoder,
including validation profiles, primer3 parameters, sequence quality constraints,
and encoding schemes.

Validation profiles are loaded from JSON files:
- default_profiles.json (shipped with package)
- user_profiles.json (optional, in current directory)
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional
import logging

logger = logging.getLogger(__name__)


@dataclass
class DNAEncoderConfig:
    """
    DNA encoder configuration containing all parameters controlling
    the encoding and sequence validation process.

    Attributes:
        primer3_mv_conc: Mg2+ ion concentration in mM for primer3 calculations
        primer3_dv_conc: Divalent ion concentration in mM for primer3 calculations
        primer3_dntp_conc: dNTP concentration in mM for primer3 calculations
        primer3_dna_conc: DNA concentration in nM for primer3 calculations

        min_gc: Minimum GC content (0.0-1.0)
        max_gc: Maximum GC content (0.0-1.0)
        min_tm: Minimum melting temperature in °C
        max_tm: Maximum melting temperature in °C
        max_hairpin_tm: Maximum Tm for hairpin structures in °C
        max_homodimer_tm: Maximum Tm for homodimers in °C

        window_size: Sequence quality analysis window size
        max_homopolymer_length: Maximum length of identical nucleotide runs
        max_backtrack_attempts: Maximum number of backtracking attempts

        bit_to_nucleotides: Dictionary mapping bits to available nucleotide pairs
        enable_deterministic: Whether to enable deterministic generation mode
        default_seed: Default seed for pseudorandom generator

        log_level: Logging level
        enable_progress_logging: Whether to enable progress logging
    """
    
    # Primer3 parameters - used for thermodynamic calculations
    primer3_mv_conc: float = 50.0
    primer3_dv_conc: float = 4.0
    primer3_dntp_conc: float = 0.5
    primer3_dna_conc: float = 50.0
    primer3_temp_c: float = 37.0

    # DNA sequence quality constraints
    min_gc: float = 0.45
    max_gc: float = 0.55
    min_tm: float = 55.0
    max_tm: float = 65.0
    max_hairpin_tm: float = 30.0
    max_homodimer_tm: float = 30.0

    # Encoding algorithm parameters
    window_size: int = 20
    max_homopolymer_length: int = 4
    max_dinucleotide_repeats: int = 2
    max_3prime_gc: int = 3
    max_backtrack_attempts: int = 1000
    enable_backtrack_heuristics: bool = True
    heuristic_weights: Dict[str, float] = field(default_factory=lambda: {
        'gc_balance': 1.0,
        'diversity': 0.05,
        'pair_repeat': 0.3,
        'pair_complement': 0.2,
        'novelty': 0.01
    })

    # Bit encoding scheme to nucleotide pairs
    # Each bit ('0' or '1') is mapped to a list of available nucleotide pairs
    bit_to_nucleotides: Dict[str, List[str]] = None

    # Deterministic mode parameters
    enable_deterministic: bool = True
    default_seed: int = None

    # Validation profiles and rules (consistent with dna_generator)
    validation_profile: Optional[str] = None
    validation_rules: Dict[str, bool] = field(default_factory=lambda: {
        'gc_content': True,
        'melting_temperature': True,
        'hairpin_structures': True,
        'homodimer_structures': True,
        'homopolymer_runs': True,
        'dinucleotide_repeats': True,
        'three_prime_stability': True,
    })
    
    # Parametry logowania
    log_level: str = "INFO"
    enable_progress_logging: bool = True
    
    def __post_init__(self):
        """
        Inicjalizacja po utworzeniu instancji - ustawia domyślne wartości
        dla parametrów, które nie mogą być ustawione w definicji klasy.
        """
        if self.bit_to_nucleotides is None:
            # Domyślny schemat kodowania - każdy bit mapowany na 8 par nukleotydów
            self.bit_to_nucleotides = {
                '0': ['TA', 'TT', 'GC', 'CC', 'AC', 'AG', 'GT', 'CT'],
                '1': ['AT', 'AA', 'CG', 'GG', 'TC', 'TG', 'GA', 'CA']
            }
        
        if self.validation_profile:
            # Load profiles from JSON files
            from .profile_loader import get_profile_loader

            loader = get_profile_loader()
            profile = loader.get_profile(self.validation_profile)

            if profile is None:
                available = ', '.join(loader.list_profiles().keys())
                raise ValueError(
                    f"Nieznany profil walidacji '{self.validation_profile}'. Dostępne: {available}"
                )

            profile_rules = profile.get('rules')
            if profile_rules:
                # Kopia, aby użytkownik mógł lokalnie nadpisać po __post_init__
                self.validation_rules = profile_rules.copy()

            profile_params = profile.get('params', {})
            for attr, value in profile_params.items():
                if hasattr(self, attr):
                    setattr(self, attr, value)
                else:
                    logger.debug(
                        "Parametr '%s' z profilu '%s' został pominięty (brak w konfiguracji)",
                        attr,
                        self.validation_profile,
                    )
            
        # Walidacja parametrów
        self._validate_config()
        
        # Konfiguracja loggingu
        self._setup_logging()
    
    def _validate_config(self) -> None:
        """
        Waliduje poprawność parametrów konfiguracji.
        
        Raises:
            ValueError: Gdy parametry mają nieprawidłowe wartości
        """
        if not (0.0 <= self.min_gc <= 1.0):
            raise ValueError(f"min_gc musi być w zakresie 0.0-1.0, otrzymano: {self.min_gc}")
        
        if not (0.0 <= self.max_gc <= 1.0):
            raise ValueError(f"max_gc musi być w zakresie 0.0-1.0, otrzymano: {self.max_gc}")
        
        if self.min_gc >= self.max_gc:
            raise ValueError(f"min_gc ({self.min_gc}) musi być mniejsze od max_gc ({self.max_gc})")
        
        if self.min_tm >= self.max_tm:
            raise ValueError(f"min_tm ({self.min_tm}) musi być mniejsze od max_tm ({self.max_tm})")
        
        if self.window_size <= 0:
            raise ValueError(f"window_size musi być dodatnie, otrzymano: {self.window_size}")
        
        if self.max_homopolymer_length <= 0:
            raise ValueError(f"max_homopolymer_length musi być dodatnie, otrzymano: {self.max_homopolymer_length}")

        if self.max_dinucleotide_repeats < 0:
            raise ValueError(
                f"max_dinucleotide_repeats musi być >= 0, otrzymano: {self.max_dinucleotide_repeats}"
            )

        if not (0 <= self.max_3prime_gc <= 5):
            raise ValueError(f"max_3prime_gc musi być w zakresie 0-5, otrzymano: {self.max_3prime_gc}")

        # Walidacja flag reguł
        required_keys = {
            'gc_content',
            'melting_temperature',
            'hairpin_structures',
            'homodimer_structures',
            'homopolymer_runs',
            'dinucleotide_repeats',
            'three_prime_stability',
        }
        missing_keys = required_keys - set(self.validation_rules.keys())
        if missing_keys:
            raise ValueError(
                f"Missing validation flags: {', '.join(sorted(missing_keys))}"
            )
        
        # Validate encoding scheme
        if self.bit_to_nucleotides:
            required_bits = {'0', '1'}
            available_bits = set(self.bit_to_nucleotides.keys())
            if not required_bits.issubset(available_bits):
                missing = required_bits - available_bits
                raise ValueError(f"Encoding scheme must contain mappings for bits: {missing}")

            # Check if all nucleotide pairs are valid
            valid_nucleotides = {'A', 'T', 'G', 'C'}
            for bit, pairs in self.bit_to_nucleotides.items():
                for pair in pairs:
                    if len(pair) != 2 or not all(n in valid_nucleotides for n in pair):
                        raise ValueError(f"Invalid nucleotide pair '{pair}' for bit '{bit}'")

        logger.info("Configuration validated successfully")

    def _setup_logging(self) -> None:
        """Configure logging system according to settings."""
        numeric_level = getattr(logging, self.log_level.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f'Invalid logging level: {self.log_level}')
        
        logging.basicConfig(
            level=numeric_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
    
    def get_primer3_params(self) -> Dict[str, float]:
        """
        Returns primer3 parameters in dictionary format.

        Returns:
            Dict containing primer3 parameters ready to use in primer3 functions
        """
        return {
            'mv_conc': self.primer3_mv_conc,
            'dv_conc': self.primer3_dv_conc,
            'dntp_conc': self.primer3_dntp_conc,
            'dna_conc': self.primer3_dna_conc,
            'temp_c': self.primer3_temp_c
        }
    
    def get_quality_constraints(self) -> Dict[str, float]:
        """
        Returns sequence quality constraints in dictionary format.

        Returns:
            Dict containing all quality parameters
        """
        return {
            'min_gc': self.min_gc,
            'max_gc': self.max_gc,
            'min_tm': self.min_tm,
            'max_tm': self.max_tm,
            'max_hairpin_tm': self.max_hairpin_tm,
            'max_homodimer_tm': self.max_homodimer_tm,
            'max_homopolymer_length': self.max_homopolymer_length,
            'max_dinucleotide_repeats': self.max_dinucleotide_repeats,
            'max_3prime_gc': self.max_3prime_gc,
        }
    
    def copy_with_overrides(self, **kwargs) -> 'DNAEncoderConfig':
        """
        Creates a copy of configuration with overridden parameters.

        Args:
            **kwargs: Parameters to override

        Returns:
            New DNAEncoderConfig instance with overridden values
        """
        # Get all current values
        current_values = self.__dict__.copy()

        if 'validation_rules' in current_values and isinstance(current_values['validation_rules'], dict):
            current_values['validation_rules'] = current_values['validation_rules'].copy()
        
        # Nadpisz wybrane wartości
        current_values.update(kwargs)

        # Usuń bit_to_nucleotides jeśli nie został jawnie nadpisany
        # żeby __post_init__ mógł ustawić domyślną wartość
        if 'bit_to_nucleotides' not in kwargs and current_values.get('bit_to_nucleotides') is None:
            current_values.pop('bit_to_nucleotides', None)
        
        # Utwórz nową instancję
        return DNAEncoderConfig(**current_values)


# Domyślna instancja konfiguracji
DEFAULT_CONFIG = DNAEncoderConfig()
