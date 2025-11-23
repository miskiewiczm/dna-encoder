#!/usr/bin/env python3
"""
Entry point for running the package as a module.

Usage:
    python -m dna_encoder
    python -m dna_encoder --test
"""

import sys
import argparse
import logging
from pathlib import Path
import unittest

from .config import DNAEncoderConfig
from .profile_loader import get_profile_loader
from .encoder import DNAEncoder
from .decoder import DNADecoder
from .utils import format_sequence_for_display


def configure_logging(quiet: bool, verbose: bool) -> None:
    level = logging.INFO
    if quiet:
        level = logging.ERROR
    elif verbose:
        level = logging.DEBUG
    logging.basicConfig(level=level, format='%(levelname)s:%(name)s:%(message)s', force=True)
    logging.getLogger('dna_encoder').setLevel(level)


def run_tests() -> bool:
    """Run unit tests via unittest discovery with package context."""
    pkg_dir = Path(__file__).resolve().parent
    top_level = pkg_dir.parent
    # Ensure package import context (for relative imports in tests)
    if str(top_level) not in sys.path:
        sys.path.insert(0, str(top_level))
    # Discover tests named test*.py in package directory, define top_level_dir
    suite = unittest.defaultTestLoader.discover(
        start_dir=str(pkg_dir), pattern="test*.py", top_level_dir=str(top_level)
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    success = result.wasSuccessful()
    if success:
        print("\n✓ All tests passed successfully!")
    else:
        print("\n✗ Some tests failed.")
    return success


def show_info():
    """Show package information."""
    try:
        from . import __version__, __description__, __author__
        print(f"DNA Encoder Refactored v{__version__}")
        print(f"Author: {__author__}")
        print(f"Description: {__description__}")
        print()
        print("Available commands:")
        print("  python -m dna_encoder --test    # Run tests")
        print("  python -m dna_encoder --help    # Show help")
        print()
        print("Code usage example:")
        print("  from dna_encoder import quick_encode")
        print("  result = quick_encode('Hello World!')")
        print("  if result.success:")
        print("      print(f'DNA: {result.sequence}')")
        print()
        print("See tests in test_encoder.py for more usage examples")
    except ImportError as e:
        print(f"Import error: {e}")
        return False
    return True


def build_encoder_from_args(args: argparse.Namespace) -> DNAEncoder:
    config_kwargs = {
        'min_gc': args.min_gc,
        'max_gc': args.max_gc,
        'min_tm': args.min_tm,
        'max_tm': args.max_tm,
        'max_hairpin_tm': args.max_hairpin,
        'max_homodimer_tm': args.max_homodimer,
        'window_size': args.window_size,
        'enable_deterministic': not getattr(args, 'random', False),
        'enable_backtrack_heuristics': not getattr(args, 'no_heuristics', False)
    }

    if args.profile:
        config_kwargs['validation_profile'] = args.profile

    if args.seed is not None:
        config_kwargs['default_seed'] = args.seed

    config = DNAEncoderConfig(**config_kwargs)

    using_custom_profile = args.profile in (None, 'user')
    if using_custom_profile:
        overrides = {
            'gc_content': not args.no_gc_check,
            'melting_temperature': not args.no_tm_check,
            'hairpin_structures': not args.no_hairpin_check,
            'homodimer_structures': not args.no_homodimer_check,
            'homopolymer_runs': not args.no_homopolymer_check,
            'dinucleotide_repeats': not args.no_dinucleotide_check,
            'three_prime_stability': not args.no_3prime_check,
        }
        config.validation_rules = {**config.validation_rules, **overrides}

    if args.profile == 'user':
        if not args.profile_file:
            raise ValueError("--profile user requires --profile-file <path to JSON> option")
        with open(args.profile_file, 'r', encoding='utf-8') as fh:
            data = json.load(fh)
        rules = data.get('rules')
        if isinstance(rules, dict):
            config.validation_rules = {**config.validation_rules, **rules}
        params = data.get('params')
        if isinstance(params, dict):
            for key, value in params.items():
                if hasattr(config, key):
                    setattr(config, key, value)

    return DNAEncoder(config)


def dna_result_to_dict(result, raw: bool = False) -> dict:
    data = {
        'success': result.success,
        'sequence': result.sequence,
        'error_message': result.error_message,
        'encoding_stats': result.encoding_stats,
        'validation_details': result.validation_details,
    }
    return data


def format_encode_output(result, args) -> str:
    if not result.success:
        if args.format == 'json':
            return json.dumps(dna_result_to_dict(result, raw=args.json_raw), indent=2, ensure_ascii=False)
        message = [f"✗ Encoding failed: {result.error_message or 'unknown error'}"]
        if result.validation_details and not args.sequences_only:
            message.append(f"Details: {result.validation_details}")
        return "\n".join(message)

    if args.format == 'json':
        payload = dna_result_to_dict(result, raw=args.json_raw)
        if args.sequences_only:
            payload = payload.get('sequence', '')
        return json.dumps(payload, indent=2, ensure_ascii=False)

    if args.format == 'fasta':
        if not result.sequence:
            return ''
        header = ">encoded|success=1"
        body = result.sequence if args.sequences_only else result.sequence
        return f"{header}\n{body}"

    # format text
    if args.sequences_only:
        return result.sequence or ''

    lines = ["✓ Encoding completed successfully"]
    lines.append(format_sequence_for_display(result.sequence))
    stats = result.encoding_stats or {}
    lines.append("Statistics:")
    lines.append(f"  Attempts: {stats.get('attempts_used')}")
    lines.append(f"  Backtracks: {stats.get('backtrack_count')}")
    lines.append(f"  Total choices: {stats.get('total_attempts')}")
    lines.append(f"  Max depth: {stats.get('max_depth_reached')}")
    analysis = stats.get('sequence_analysis') if isinstance(stats, dict) else None
    if analysis:
        complexity = analysis.get('complexity_k2')
        if complexity is not None:
            lines.append(f"  Complexity k2: {complexity:.3f}")
        longest = analysis.get('longest_homopolymer')
        if longest and isinstance(longest, dict):
            lines.append(f"  Longest homopolymer: {longest.get('nucleotide', '?')} ({longest.get('length', 0)} nt)")
    return "\n".join(lines)


def emit_output(text: str, args) -> None:
    if text is None:
        text = ''
    if args.output:
        out_text = text if text.endswith('\n') else text + '\n'
        Path(args.output).write_text(out_text, encoding='utf-8')
        if not args.quiet:
            print(f"Output saved to {args.output}")
    else:
        print(text)


def export_windows_to_csv(sequence: str, window_size: int, validator, csv_path: str):
    """
    Export sliding window analysis to CSV file.

    Args:
        sequence: Full DNA sequence to analyze
        window_size: Size of sliding window (typically 20)
        validator: DNAValidator instance for quality checks
        csv_path: Path to output CSV file
    """
    import csv

    rows = []

    # Analyze each window (sliding with step=1)
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]

        # Get comprehensive metrics for this window
        metrics = validator.validate_sequence(window)

        row = {
            'window_start': i,
            'window_end': i + window_size,
            'sequence': window,
            'gc_content': f"{metrics.gc_content:.4f}",
            'melting_temperature': f"{metrics.melting_temperature:.2f}",
            'hairpin_tm': f"{metrics.hairpin_tm:.2f}",
            'homodimer_tm': f"{metrics.homodimer_tm:.2f}",
            'has_homopolymers': metrics.has_homopolymers,
            'has_dinucleotide_repeats': metrics.has_dinucleotide_repeats,
            'three_prime_gc_count': metrics.three_prime_gc_count,
            'is_valid': metrics.is_valid,
            'longest_homopolymer': str(metrics.longest_homopolymer) if metrics.longest_homopolymer else '',
            'max_dinucleotide_repeat': str(metrics.max_dinucleotide_repeat) if metrics.max_dinucleotide_repeat else '',
        }
        rows.append(row)

    # Write to CSV
    fieldnames = [
        'window_start', 'window_end', 'sequence',
        'gc_content', 'melting_temperature', 'hairpin_tm', 'homodimer_tm',
        'has_homopolymers', 'has_dinucleotide_repeats', 'three_prime_gc_count',
        'is_valid', 'longest_homopolymer', 'max_dinucleotide_repeat'
    ]

    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Analiza sliding window zapisana do {csv_path} ({len(rows)} okien)")


def main():
    """Główna funkcja modułu."""
    parser = argparse.ArgumentParser(
        description="DNA Encoder - Pakiet do kodowania danych w sekwencjach DNA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Przykłady użycia:
  python -m dna_encoder                                    # Pokaż informacje o pakiecie
  python -m dna_encoder --test                             # Uruchom testy jednostkowe
  python -m dna_encoder --encode "test" --initial ATGCATGC      # Zakoduj tekst
  python -m dna_encoder --encode-bits "01001100" --initial ATGC # Zakoduj bity
  python -m dna_encoder --decode "ATGC..." --initial ATGC       # Zdekoduj sekwencję
        """
    )

    parser.add_argument('--test', action='store_true', help='Uruchom testy jednostkowe')
    parser.add_argument('--version', action='store_true', help='Pokaż wersję pakietu')
    parser.add_argument('--encode', metavar='TEXT', help='Zakoduj podany tekst do DNA')
    parser.add_argument('--encode-bits', metavar='BITS', help='Zakoduj ciąg bitów do DNA (np. "01001100")')
    parser.add_argument('--decode', metavar='DNA', help='Zdekoduj podaną sekwencję DNA')
    parser.add_argument('--initial', '-i', metavar='SEQ', required=False,
                        help='Sekwencja początkowa DNA (wymagana dla --encode i --encode-bits)')
    parser.add_argument('--max-attempts', type=int, default=5, help='Limit prób kodowania (domyślnie 5)')
    parser.add_argument('--profile', choices=sorted(get_profile_loader().list_profiles().keys()),
                        help='Wybierz profil walidacji (np. strict, relaxed, sequence_only)')
    parser.add_argument('--profile-file', metavar='JSON', help='Ścieżka do profilu użytkownika (dla --profile user)')
    parser.add_argument('--min-gc', type=float, default=0.45, help='Minimalna zawartość GC (0.0-1.0)')
    parser.add_argument('--max-gc', type=float, default=0.55, help='Maksymalna zawartość GC (0.0-1.0)')
    parser.add_argument('--min-tm', type=float, default=55.0, help='Minimalna temperatura topnienia (°C)')
    parser.add_argument('--max-tm', type=float, default=65.0, help='Maksymalna temperatura topnienia (°C)')
    parser.add_argument('--max-hairpin', type=float, default=30.0, help='Limit Tm hairpinów (°C)')
    parser.add_argument('--max-homodimer', type=float, default=30.0, help='Limit Tm homodimerów (°C)')
    parser.add_argument('--window-size', type=int, default=20, help='Rozmiar okna analizy walidatora')
    parser.add_argument('--no-gc-check', action='store_true', help='Wyłącz kontrolę zawartości GC')
    parser.add_argument('--no-tm-check', action='store_true', help='Wyłącz kontrolę temperatury topnienia')
    parser.add_argument('--no-hairpin-check', action='store_true', help='Wyłącz kontrolę hairpinów')
    parser.add_argument('--no-homodimer-check', action='store_true', help='Wyłącz kontrolę homodimerów')
    parser.add_argument('--no-homopolymer-check', action='store_true', help='Wyłącz kontrolę homopolimerów')
    parser.add_argument('--no-dinucleotide-check', action='store_true', help='Wyłącz kontrolę powtórzeń dinukleotydów')
    parser.add_argument('--no-3prime-check', action='store_true', help='Wyłącz kontrolę końca 3\' DNA')
    parser.add_argument('--seed', type=int, help='Wymuś konkretny seed deterministyczny')
    parser.add_argument('--random', action='store_true', help='Wyłącz tryb deterministyczny (użyj losowego RNG)')
    parser.add_argument('--no-heuristics', action='store_true', help='Wyłącz heurystyki wyboru par nukleotydów')
    parser.add_argument('--as-bits', action='store_true', help='Podczas dekodowania zwróć ciąg bitów')
    parser.add_argument('--output', '-o', metavar='PLIK', help='Zapisz wynik do pliku')
    parser.add_argument('--format', choices=['text', 'json', 'fasta'], default='text', help='Format wyjścia (domyślnie text)')
    parser.add_argument('--sequences-only', action='store_true', help='Wyświetl tylko sekwencje (bez statystyk)')
    parser.add_argument('--json-raw', action='store_true', help='W formacie JSON zwróć surowe metryki')
    parser.add_argument('--csv-file', metavar='PLIK', help='Eksportuj analizę sliding window do CSV (okna 20nt, krok=1)')
    parser.add_argument('--quiet', action='store_true', help='Ogranicz logowanie do błędów')
    parser.add_argument('--verbose', action='store_true', help='Włącz szczegółowe logowanie (DEBUG)')

    args = parser.parse_args()

    configure_logging(args.quiet, args.verbose)

    if args.version:
        try:
            from . import __version__
            print(f"dna_encoder v{__version__}")
        except ImportError:
            print("Nie można określić wersji pakietu.")
        return

    performed_action = False

    if args.test:
        performed_action = True
        print("Uruchamianie testów...")
        success = run_tests()
        if not success:
            sys.exit(1)

    if args.encode:
        performed_action = True
        if not args.initial:
            print("✗ Error: --initial is required for encoding")
            print("Example: python -m dna_encoder --encode 'test' --initial ATGCATGC")
            sys.exit(1)
        encoder = build_encoder_from_args(args)
        initial_sequence = args.initial
        result = encoder.encode_text(initial_sequence, args.encode, max_attempts=args.max_attempts)
        header = f"Encoding text: '{args.encode}'\nInitial sequence: {initial_sequence}"
        if not args.output and not args.sequences_only and not args.quiet:
            print(header)
        output_text = format_encode_output(result, args)
        if not result.success:
            emit_output(output_text, args)
            sys.exit(1)
        emit_output(output_text, args)

        # Export CSV if requested
        if args.csv_file and result.sequence:
            window_size = encoder.config.window_size
            export_windows_to_csv(
                result.sequence,
                window_size=window_size,
                validator=encoder.validator,
                csv_path=args.csv_file
            )

    if args.encode_bits:
        performed_action = True
        if not args.initial:
            print("✗ Error: --initial is required for encoding bits")
            print("Example: python -m dna_encoder --encode-bits '01001100' --initial ATGCATGC")
            sys.exit(1)

        # Validate bit string
        if not all(c in '01' for c in args.encode_bits):
            print("✗ Error: --encode-bits requires a string containing only '0' and '1'")
            print(f"Received: '{args.encode_bits}'")
            sys.exit(1)

        encoder = build_encoder_from_args(args)
        initial_sequence = args.initial
        result = encoder.encode_bits(initial_sequence, args.encode_bits, max_attempts=args.max_attempts)
        header = f"Encoding bits: '{args.encode_bits}' ({len(args.encode_bits)} bits)\nInitial sequence: {initial_sequence}"
        if not args.output and not args.sequences_only and not args.quiet:
            print(header)
        output_text = format_encode_output(result, args)
        if not result.success:
            emit_output(output_text, args)
            sys.exit(1)
        emit_output(output_text, args)

        # Export CSV if requested
        if args.csv_file and result.sequence:
            window_size = encoder.config.window_size
            export_windows_to_csv(
                result.sequence,
                window_size=window_size,
                validator=encoder.validator,
                csv_path=args.csv_file
            )

    if args.decode:
        performed_action = True
        if not args.initial:
            print("✗ Error: --initial is required for decoding")
            print("Example: python -m dna_encoder --decode 'ATGC...' --initial ATGCATGC")
            sys.exit(1)
        decoder = DNADecoder()
        try:
            decoded = decoder.decode_sequence(args.decode, args.initial)
            decode_output = decoded
            if args.as_bits:
                decode_output = decoded
            if args.output:
                emit_output(decode_output, args)
            else:
                print(f"Decoding sequence (len={len(args.decode)}):")
                print("Bits:" if args.as_bits else "Text:", decode_output)
        except Exception as exc:
            print("✗ Decoding failed with error:", exc)
            sys.exit(1)

    if not performed_action:
        show_info()


if __name__ == "__main__":
    main()
