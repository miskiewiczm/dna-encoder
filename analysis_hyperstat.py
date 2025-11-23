"""
Hyper-statistical analysis of DNA encoder performance.

Tests all primers from primers.txt with random bit sequences to generate
comprehensive success rate statistics across all validation profiles.

Analysis:
- 500+ primers from primers.txt
- 100 random 200-bit sequences per primer
- All 4 validation profiles
- Random mode (non-deterministic)

Output: CSV file with primer × profile success rates
"""

import random
import time
import csv
import logging
import warnings
from pathlib import Path
from typing import Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
from os import cpu_count
from dna_encoder import DNAEncoder, DNAEncoderConfig

# Disable all logging from libraries
logging.basicConfig(level=logging.CRITICAL)
logging.getLogger('dna_encoder').setLevel(logging.CRITICAL)
logging.getLogger('dna_commons').setLevel(logging.CRITICAL)
warnings.filterwarnings('ignore')


def generate_random_bits(length: int) -> str:
    """Generate random bit string of specified length."""
    return ''.join(random.choice('01') for _ in range(length))


def process_primer(args: Tuple[str, int, int, int]) -> Tuple[str, Dict[str, int]]:
    """
    Process a single primer with random sequences.

    Args:
        args: (primer, num_sequences, bit_length, max_attempts)

    Returns:
        (primer, {profile: success_count})
    """
    primer, num_sequences, bit_length, max_attempts = args

    profiles = ['sequence_only', 'relaxed', 'pcr_friendly', 'strict']

    # Initialize encoders for each profile
    encoders = {
        profile: DNAEncoder(DNAEncoderConfig(
            validation_profile=profile,
            enable_deterministic=False
        ))
        for profile in profiles
    }

    # Generate random sequences
    random_sequences = [generate_random_bits(bit_length) for _ in range(num_sequences)]

    # Test each profile
    primer_results = {}
    for profile in profiles:
        encoder = encoders[profile]
        successes = 0

        for seq in random_sequences:
            result = encoder.encode_bits(primer, seq, max_attempts=max_attempts)
            if result.success:
                successes += 1

        primer_results[profile] = successes

    return (primer, primer_results)


def run_hyperstat_analysis(
    primers_file: Path,
    num_sequences: int = 100,
    bit_length: int = 200,
    max_attempts: int = 10,
    output_file: str = 'hyperstat_results.csv',
    num_workers: int = None
) -> Dict[str, Dict[str, int]]:
    """
    Run comprehensive statistical analysis with multiprocessing.

    Args:
        primers_file: Path to primers.txt
        num_sequences: Number of random sequences to test per primer
        bit_length: Length of each random bit sequence
        max_attempts: Maximum encoding attempts per sequence
        output_file: Output CSV filename
        num_workers: Number of parallel workers (None = cpu_count)

    Returns:
        Dictionary: {primer: {profile: success_count}}
    """
    # Load primers
    print(f"Loading primers from {primers_file}...")
    with open(primers_file) as f:
        primers = [line.strip() for line in f if line.strip()]

    if num_workers is None:
        num_workers = cpu_count()

    print(f"Loaded {len(primers)} primers")
    print(f"Testing {num_sequences} random {bit_length}-bit sequences per primer")
    print(f"Profiles: sequence_only, relaxed, pcr_friendly, strict")
    print(f"Total tests: {len(primers)} × {num_sequences} × 4 = {len(primers) * num_sequences * 4:,}")
    print(f"Using {num_workers} parallel workers")
    print(f"\nStarting analysis...\n")

    # Prepare arguments for workers
    worker_args = [(primer, num_sequences, bit_length, max_attempts) for primer in primers]

    # Results storage
    results = {}
    total_primers = len(primers)
    start_time = time.time()

    # Process primers in parallel using ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit all tasks
        future_to_primer = {
            executor.submit(process_primer, args): args[0]
            for args in worker_args
        }

        # Process results as they complete
        for idx, future in enumerate(as_completed(future_to_primer), 1):
            primer, primer_results = future.result()
            results[primer] = primer_results

            # Progress reporting
            if idx % 10 == 0 or idx == total_primers:
                elapsed = time.time() - start_time
                avg_time_per_primer = elapsed / idx
                remaining_primers = total_primers - idx
                eta = remaining_primers * avg_time_per_primer

                print(f"Progress: {idx}/{total_primers} primers "
                      f"({idx/total_primers*100:.1f}%) | "
                      f"Elapsed: {elapsed:.1f}s | "
                      f"ETA: {eta:.1f}s")

                # Show sample results for last processed primer
                print(f"  Primer {primer[:20]}...:")
                for profile in ['sequence_only', 'relaxed', 'pcr_friendly', 'strict']:
                    rate = primer_results[profile] / num_sequences * 100
                    print(f"    {profile:15}: {primer_results[profile]:3}/{num_sequences} ({rate:5.1f}%)")
                print()

    total_time = time.time() - start_time
    print(f"\nAnalysis complete!")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
    print(f"Average time per primer: {total_time/len(primers):.2f}s")
    print(f"Speedup: ~{3400/total_time:.1f}x faster than sequential")

    # Save results to CSV
    save_results_csv(results, output_file, num_sequences)
    print(f"\nResults saved to: {output_file}")

    # Print summary statistics
    print_summary_statistics(results, num_sequences)

    return results


def save_results_csv(
    results: Dict[str, Dict[str, int]],
    output_file: str,
    num_sequences: int
):
    """Save results to CSV file."""
    profiles = ['sequence_only', 'relaxed', 'pcr_friendly', 'strict']

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)

        # Header
        header = ['primer'] + profiles + [f'{p}_rate' for p in profiles]
        writer.writerow(header)

        # Data rows
        for primer, profile_results in results.items():
            row = [primer]
            # Success counts
            row.extend([profile_results[p] for p in profiles])
            # Success rates
            row.extend([f"{profile_results[p]/num_sequences*100:.1f}" for p in profiles])
            writer.writerow(row)


def print_summary_statistics(results: Dict[str, Dict[str, int]], num_sequences: int):
    """Print summary statistics."""
    profiles = ['sequence_only', 'relaxed', 'pcr_friendly', 'strict']

    print(f"\n{'='*80}")
    print(f"SUMMARY STATISTICS")
    print(f"{'='*80}\n")

    # Calculate aggregate statistics per profile
    print(f"{'Profile':<15} {'Avg Success Rate':<20} {'Min':<10} {'Max':<10} {'Primers with 100%':<20}")
    print(f"{'-'*80}")

    for profile in profiles:
        # Collect all success rates for this profile
        rates = [results[primer][profile] / num_sequences * 100 for primer in results]

        avg_rate = sum(rates) / len(rates)
        min_rate = min(rates)
        max_rate = max(rates)
        perfect_primers = sum(1 for r in rates if r == 100.0)

        print(f"{profile:<15} {avg_rate:>6.1f}% "
              f"({sum(results[p][profile] for p in results)}/{len(results)*num_sequences:>6}) "
              f"{min_rate:>6.1f}% "
              f"{max_rate:>6.1f}% "
              f"{perfect_primers:>6} ({perfect_primers/len(results)*100:5.1f}%)")

    print(f"\n{'='*80}")

    # Find best and worst performing primers per profile
    print(f"\nBEST PERFORMING PRIMERS (per profile):\n")

    for profile in profiles:
        # Sort primers by success rate for this profile
        sorted_primers = sorted(
            results.items(),
            key=lambda x: x[1][profile],
            reverse=True
        )

        print(f"{profile}:")
        for i, (primer, profile_results) in enumerate(sorted_primers[:3], 1):
            rate = profile_results[profile] / num_sequences * 100
            print(f"  {i}. {primer[:30]:30} {profile_results[profile]:3}/{num_sequences} ({rate:5.1f}%)")
        print()

    print(f"\nWORST PERFORMING PRIMERS (per profile):\n")

    for profile in profiles:
        # Sort primers by success rate for this profile
        sorted_primers = sorted(
            results.items(),
            key=lambda x: x[1][profile]
        )

        print(f"{profile}:")
        for i, (primer, profile_results) in enumerate(sorted_primers[:3], 1):
            rate = profile_results[profile] / num_sequences * 100
            print(f"  {i}. {primer[:30]:30} {profile_results[profile]:3}/{num_sequences} ({rate:5.1f}%)")
        print()


def main():
    """Run hyper-statistical analysis."""
    # Configuration
    primers_file = Path(__file__).parent / 'primers.txt'

    if not primers_file.exists():
        print(f"Error: primers.txt not found at {primers_file}")
        return

    # Run analysis
    results = run_hyperstat_analysis(
        primers_file=primers_file,
        num_sequences=100,      # 100 random sequences per primer
        bit_length=200,         # 200 bits per sequence
        max_attempts=10,        # Max 10 attempts per encoding
        output_file='hyperstat_results.csv'
    )

    print("\nAnalysis complete! Check hyperstat_results.csv for detailed results.")


if __name__ == '__main__':
    main()
