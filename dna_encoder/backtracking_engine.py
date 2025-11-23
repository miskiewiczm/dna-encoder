"""
Backtracking engine for DNA encoding.

This module contains the core backtracking algorithm logic used for encoding
data into DNA sequences with quality control and heuristic optimization.
"""

import logging
from typing import Optional, List, Dict, Any, Tuple

from .config import DNAEncoderConfig
from dna_commons import DNAValidator, DeterministicRandom, generate_seed_from_string, SequenceAnalyzer
from .exceptions import BacktrackingError

logger = logging.getLogger(__name__)

# DNA complement table for heuristics
COMPLEMENT_TABLE = str.maketrans('ATGC', 'TACG')

# Bit to nucleotide mapping - must match encoder.py
BIT_TO_NUCLEOTIDES = {
    '0': ['TA', 'TT', 'GC', 'CC', 'AC', 'AG', 'GT', 'CT'],
    '1': ['AT', 'AA', 'CG', 'GG', 'TC', 'TG', 'GA', 'CA']
}


class BacktrackingEngine:
    """
    Core backtracking engine for DNA sequence generation.

    This class implements the backtracking algorithm used to encode bit sequences
    into DNA while maintaining biochemical quality constraints.
    """

    def __init__(self, config: DNAEncoderConfig, validator: DNAValidator,
                 analyzer: SequenceAnalyzer):
        """
        Initialize the backtracking engine.

        Args:
            config: Configuration for generation parameters
            validator: DNA sequence validator
            analyzer: Sequence analysis utility
        """
        self.config = config
        self.validator = validator
        self.analyzer = analyzer

    def encode_bits_single_attempt(
        self,
        initial_sequence: str,
        bit_sequence: str
    ) -> Tuple[Optional[str], Dict[str, Any]]:
        """
        Single attempt to encode bits using backtracking algorithm.

        Args:
            initial_sequence: Starting DNA sequence
            bit_sequence: Bit string to encode

        Returns:
            Tuple of (generated sequence or None, detailed statistics)
        """
        # Setup encoding environment
        seed, rng, base_seed = self._setup_encoding_environment(initial_sequence, bit_sequence)
        stats = self._initialize_stats(initial_sequence, bit_sequence, seed)

        # Execute main encoding algorithm
        final_sequence, backtrack_count = self._execute_encoding_algorithm(
            initial_sequence, bit_sequence, base_seed, rng, stats
        )

        # Finalize and return results
        return self._finalize_encoding_result(final_sequence, backtrack_count, stats)

    def _setup_encoding_environment(self, initial_sequence: str, bit_sequence: str) -> Tuple[Optional[int], DeterministicRandom, Optional[int]]:
        """Setup encoding environment with seed and RNG."""
        seed = None
        if self.config.enable_deterministic:
            seed = self.config.default_seed
            if seed is None:
                seed = generate_seed_from_string(initial_sequence, additional_data=bit_sequence)
            logger.debug(f"Deterministic mode: seed = {seed}")

        rng = DeterministicRandom(
            seed=seed,
            deterministic=self.config.enable_deterministic
        )
        base_seed = seed if self.config.enable_deterministic else None
        return seed, rng, base_seed

    def _execute_encoding_algorithm(self, initial_sequence: str, bit_sequence: str,
                                   base_seed: Optional[int], rng: DeterministicRandom,
                                   stats: Dict[str, Any]) -> Tuple[str, int]:
        """Execute main backtracking algorithm."""
        nucleotide_pairs = [BIT_TO_NUCLEOTIDES[bit].copy() for bit in bit_sequence]
        sequence_parts = [initial_sequence]
        available_options = []

        if nucleotide_pairs:
            available_options.append(nucleotide_pairs[0].copy())

        position = 0
        backtrack_count = 0
        max_backtracks = getattr(self.config, 'max_backtrack_attempts', 1000)

        while position < len(bit_sequence):
            stats['total_attempts'] += 1

            if not available_options or not available_options[-1]:
                backtrack_result = self._handle_backtrack(
                    sequence_parts, available_options, position,
                    backtrack_count, max_backtracks, stats
                )
                if backtrack_result is None:
                    return ''.join(sequence_parts), backtrack_count
                position, backtrack_count = backtrack_result
                continue

            chosen_pair = self._choose_and_validate_pair(
                sequence_parts, available_options, position, base_seed, rng, stats
            )

            if chosen_pair:
                sequence_parts.append(chosen_pair)
                position += 1
                self._update_acceptance_stats(stats, ''.join(sequence_parts), chosen_pair, position)

                if position < len(bit_sequence):
                    available_options.append(nucleotide_pairs[position].copy())

        return ''.join(sequence_parts), backtrack_count

    def _choose_and_validate_pair(self, sequence_parts: List[str], available_options: List[List[str]],
                                 position: int, base_seed: Optional[int], rng: DeterministicRandom,
                                 stats: Dict[str, Any]) -> Optional[str]:
        """Choose and validate nucleotide pair."""
        current_options = available_options[-1]
        current_sequence = ''.join(sequence_parts)

        chosen_pair = self._choose_pair_with_heuristics(
            current_sequence=current_sequence,
            options=current_options,
            position=position,
            base_seed=base_seed,
            rng=rng
        )
        current_options.remove(chosen_pair)

        candidate_sequence = ''.join(sequence_parts + [chosen_pair])
        if self._validate_candidate(candidate_sequence, stats):
            return chosen_pair
        return None

    def _finalize_encoding_result(self, final_sequence: str, backtrack_count: int,
                                 stats: Dict[str, Any]) -> Tuple[str, Dict[str, Any]]:
        """Finalize encoding result with statistics."""
        self._finalize_stats(stats, final_sequence, backtrack_count)
        logger.debug(
            "Generated sequence (length: %s, backtracks: %s, attempts: %s)",
            len(final_sequence), backtrack_count, stats['total_attempts']
        )
        return final_sequence, stats

    def _initialize_stats(self, initial_sequence: str, bit_sequence: str,
                         seed: Optional[int]) -> Dict[str, Any]:
        """Initialize comprehensive statistics tracking."""
        return {
            'seed': seed,
            'heuristics_enabled': getattr(self.config, 'enable_backtrack_heuristics', True),
            'initial_sequence_length': len(initial_sequence),
            'target_bits': len(bit_sequence),
            'backtrack_count': 0,
            'total_attempts': 0,
            'max_depth_reached': len(initial_sequence),
            'validation_failures': {
                'gc_content': 0,
                'melting_temp': 0,
                'homopolymers': 0,
                'dinucleotide_repeats': 0,
                'three_prime_stability': 0,
                'hairpin': 0,
                'homodimer': 0,
            },
            'window_rollup': {
                'gc_min': None,
                'gc_max': None,
                'tm_min': None,
                'tm_max': None,
                'hairpin_tm_max': None,
                'homodimer_tm_max': None,
            },
            'attempt_history': [],
        }

    def _handle_backtrack(self, sequence_parts: List[str], available_options: List[List[str]],
                         position: int, backtrack_count: int, max_backtracks: int,
                         stats: Dict[str, Any]) -> Optional[Tuple[int, int]]:
        """
        Handle backtracking step.

        Returns:
            Tuple of (new_position, new_backtrack_count) or None if impossible
        """
        # Check if we can backtrack
        if len(sequence_parts) <= 1:  # Don't backtrack beyond initial sequence
            logger.debug("Cannot generate sequence - impossible encoding.")
            return None

        # Perform backtrack
        sequence_parts.pop()
        available_options.pop()
        position -= 1
        backtrack_count += 1
        stats['backtrack_count'] = backtrack_count

        # Check backtrack limit
        if backtrack_count > max_backtracks:
            logger.debug(f"Exceeded maximum backtracks: {max_backtracks}")
            return None

        return position, backtrack_count

    def _validate_candidate(self, candidate_sequence: str, stats: Dict[str, Any]) -> bool:
        """Validate candidate sequence using the central validator."""
        window_size = getattr(self.config, 'window_size', 16)

        if len(candidate_sequence) >= window_size:
            window_seq = candidate_sequence[-window_size:]
        else:
            window_seq = candidate_sequence

        return self.validator.validate_window(window_seq)

    def _update_acceptance_stats(self, stats: Dict[str, Any], candidate_sequence: str,
                               chosen_pair: str, position: int) -> None:
        """Update statistics when a candidate is accepted."""
        stats['max_depth_reached'] = max(stats['max_depth_reached'], len(candidate_sequence))
        stats['attempt_history'].append({
            'position': position,
            'pair': chosen_pair
        })

    def _finalize_stats(self, stats: Dict[str, Any], final_sequence: str,
                       backtrack_count: int) -> None:
        """Finalize statistics after successful generation."""
        stats['backtrack_count'] = backtrack_count
        stats['max_depth_reached'] = max(stats['max_depth_reached'], len(final_sequence))
        stats['final_sequence_length'] = len(final_sequence)

        if final_sequence:
            stats['sequence_analysis'] = self.analyzer.analyze_sequence(final_sequence)

    def _choose_pair_with_heuristics(
        self,
        current_sequence: str,
        options: List[str],
        position: int,
        base_seed: Optional[int],
        rng: DeterministicRandom
    ) -> str:
        """Choose nucleotide pair using heuristics adapted for 2-nt pairs."""
        if not options:
            raise BacktrackingError("No available nucleotide pairs to choose from")

        heuristics_enabled = getattr(self.config, 'enable_backtrack_heuristics', True)
        if not heuristics_enabled or len(options) <= 1:
            return self._select_with_deterministic_tiebreak(options, base_seed, position, rng)

        heuristic_context = self._prepare_heuristic_context(current_sequence)
        scored_options = self._score_all_options(options, current_sequence, heuristic_context)
        scored_options.sort(key=lambda item: item[0])
        return self._select_from_scored_options(scored_options, base_seed, position, rng)

    def _prepare_heuristic_context(self, current_sequence: str) -> Dict[str, Any]:
        """Prepare context data for heuristic evaluation."""
        current_sequence = current_sequence.upper()
        window_size = max(2, getattr(self.config, 'window_size', 2))

        # Only calculate target_gc if GC content validation is enabled
        if self.validator.rules.gc_content:
            target_gc = (self.config.min_gc + self.config.max_gc) / 2.0
        else:
            target_gc = None  # No GC targeting when validation disabled

        recent_len = max(100, window_size)
        recent_context = current_sequence[-recent_len:]

        weights = self.config.heuristic_weights or {}
        last_pair = current_sequence[-2:] if len(current_sequence) >= 2 else ""
        last_pair_complement = last_pair.translate(COMPLEMENT_TABLE) if last_pair else ""

        return {
            'window_size': window_size,
            'target_gc': target_gc,
            'recent_context': recent_context,
            'weights': weights,
            'last_pair': last_pair,
            'last_pair_complement': last_pair_complement,
            'hard_penalty_base': 1_000_000.0
        }

    def _score_all_options(self, options: List[str], current_sequence: str, context: Dict[str, Any]) -> List[Tuple[float, str]]:
        """Score all nucleotide pair options using heuristics."""
        scored: List[Tuple[float, str]] = []
        weights = context['weights']

        for pair in options:
            candidate_sequence = (current_sequence + pair).upper()
            window = candidate_sequence[-context['window_size']:]

            # Check hard constraints first
            hard_penalty = self._calculate_hard_penalties(window)
            if hard_penalty > 0:
                scored.append((context['hard_penalty_base'] + hard_penalty, pair))
                continue

            # Calculate soft heuristic scores
            score = self._calculate_heuristic_score(
                pair, window, context['target_gc'], context['recent_context'], candidate_sequence,
                context['last_pair'], context['last_pair_complement'],
                weights.get('gc_balance', 1.0), weights.get('diversity', 0.05),
                weights.get('pair_repeat', 0.3), weights.get('pair_complement', 0.2),
                weights.get('novelty', 0.01)
            )
            scored.append((score, pair))

        return scored

    def _calculate_hard_penalties(self, window: str) -> float:
        """Calculate hard constraint penalties (only for enabled rules)."""
        hard_penalty = 0.0

        # Only apply penalties for enabled validation rules
        if self.validator.rules.homopolymer_runs:
            homopoly_valid, _ = self.validator._check_homopolymer_runs(window)
            if not homopoly_valid:
                hard_penalty += 300_000.0

        if self.validator.rules.three_prime_stability:
            stability_valid = self.validator._check_3_prime_stability(window)
            if not stability_valid:
                hard_penalty += 200_000.0

        if self.validator.rules.dinucleotide_repeats:
            repeats_valid, _ = self.validator._check_dinucleotide_repeats(window)
            if not repeats_valid:
                hard_penalty += 150_000.0

        return hard_penalty

    def _calculate_heuristic_score(
        self, pair: str, window: str, target_gc: float, recent_context: str,
        candidate_sequence: str, last_pair: str, last_pair_complement: str,
        gc_weight: float, diversity_weight: float, pair_repeat_weight: float,
        pair_complement_weight: float, novelty_weight: float
    ) -> float:
        """Calculate comprehensive heuristic score for a nucleotide pair."""
        # GC content score (only if validation enabled)
        gc_score = 0.0
        if target_gc is not None:
            gc_fraction = (window.count('G') + window.count('C')) / len(window) if window else 0.0
            gc_score = gc_weight * abs(gc_fraction - target_gc)

        # Diversity score
        diversity_score = 0.0
        if window:
            diversity_score = (window.count(pair[0]) + window.count(pair[1])) / max(1, len(window))
        diversity_contribution = diversity_weight * diversity_score

        # Pair repeat score
        pair_repeat_score = 0.0
        if len(window) >= 2:
            occurrences = window.count(pair)
            pair_repeat_score = occurrences / max(1, len(window) // 2)
        immediate_repeat = 1.0 if last_pair and pair == last_pair else 0.0
        repeat_contribution = pair_repeat_weight * (pair_repeat_score + immediate_repeat)

        # Complement contribution
        complement_contribution = 0.0
        if last_pair_complement and pair == last_pair_complement:
            complement_contribution = pair_complement_weight

        # Novelty score (anti-periodicity)
        novelty_hits = 0
        for k in (4, 6, 8):
            if len(candidate_sequence) >= k:
                kmer = candidate_sequence[-k:]
                novelty_hits += recent_context.count(kmer)
        novelty_contribution = novelty_weight * novelty_hits

        return gc_score + diversity_contribution + repeat_contribution + complement_contribution + novelty_contribution

    def _select_from_scored_options(
        self, scored: List[Tuple[float, str]], base_seed: Optional[int],
        position: int, rng: DeterministicRandom
    ) -> str:
        """Select final option from scored candidates."""
        best_score = scored[0][0]

        # Non-deterministic mode: select from top candidates within margin
        if not self.config.enable_deterministic:
            margin = 0.02
            eligible = [pair for score, pair in scored if score <= best_score + margin]
            return self._select_with_deterministic_tiebreak(
                eligible, base_seed, position, rng, deterministic_override=False
            )

        # Deterministic mode: select from exactly tied candidates
        tolerance = 1e-6
        best_candidates = [pair for score, pair in scored if abs(score - best_score) <= tolerance]
        if not best_candidates:
            best_candidates = [scored[0][1]]

        return self._select_with_deterministic_tiebreak(best_candidates, base_seed, position, rng)

    def _select_with_deterministic_tiebreak(
        self,
        candidates: List[str],
        base_seed: Optional[int],
        position: int,
        rng: DeterministicRandom,
        deterministic_override: Optional[bool] = None
    ) -> str:
        """Select element from candidates using deterministic tie-breaking."""
        if not candidates:
            raise BacktrackingError("No candidates to choose from")

        deterministic_mode = (
            self.config.enable_deterministic
            if deterministic_override is None
            else deterministic_override
        )

        if len(candidates) == 1:
            return candidates[0]

        if deterministic_mode and base_seed is not None:
            # Deterministic selection based on position and seed
            position_seed = (base_seed + position) % (2**31)
            local_rng = DeterministicRandom(seed=position_seed, deterministic=True)
            return local_rng.choice(candidates)
        else:
            # Use provided RNG
            return rng.choice(candidates)


# Export the main class
__all__ = ['BacktrackingEngine', 'BIT_TO_NUCLEOTIDES']