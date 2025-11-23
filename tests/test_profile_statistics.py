"""
Comprehensive profile statistics tests.

Tests all validation profiles with 100 random primers and various data sizes
to generate real-world statistics about success rates and backtracking performance.
"""

import unittest
import random
import time
import statistics
from pathlib import Path
from typing import List, Dict, Any
from dna_encoder import DNAEncoder, DNAEncoderConfig


class TestProfileStatistics(unittest.TestCase):
    """Statistical analysis of validation profiles with random primers."""

    @classmethod
    def setUpClass(cls):
        """Load primers and prepare test data."""
        primers_file = Path(__file__).parent.parent / 'primers.txt'
        if not primers_file.exists():
            raise FileNotFoundError(f"primers.txt not found at {primers_file}")

        with open(primers_file) as f:
            all_primers = [line.strip() for line in f if line.strip()]

        # Select 100 random primers for testing
        random.seed(42)  # For reproducibility
        cls.test_primers = random.sample(all_primers, min(100, len(all_primers)))

        # Test messages of different sizes
        cls.test_messages = {
            'short': ['Hi', 'OK', 'DNA', 'Test'],
            'medium': ['Hello World', 'DNA Storage', 'Encode This'],
            'long': ['This is a longer message for testing', 'DNA encoding test data']
        }

        print(f"\n{'='*80}")
        print(f"PROFILE STATISTICS TEST")
        print(f"Testing with {len(cls.test_primers)} random primers")
        print(f"{'='*80}\n")

    def _run_profile_test(self, profile_name: str, max_attempts: int = 10) -> Dict[str, Any]:
        """
        Run comprehensive test for a single profile.

        Args:
            profile_name: Name of validation profile
            max_attempts: Maximum encoding attempts per test

        Returns:
            Dictionary with statistics
        """
        config = DNAEncoderConfig(validation_profile=profile_name)
        encoder = DNAEncoder(config)

        results = {
            'profile': profile_name,
            'total_tests': 0,
            'successes': 0,
            'failures': 0,
            'backtracks': [],
            'execution_times': [],
            'failed_primers': [],
            'by_message_size': {
                'short': {'success': 0, 'total': 0},
                'medium': {'success': 0, 'total': 0},
                'long': {'success': 0, 'total': 0}
            }
        }

        # Test with random subset of primers and all message sizes
        test_primers = random.sample(self.test_primers, 20)  # 20 primers per profile

        for primer in test_primers:
            for size_category, messages in self.test_messages.items():
                for message in messages:
                    results['total_tests'] += 1
                    results['by_message_size'][size_category]['total'] += 1

                    start_time = time.time()
                    result = encoder.encode_text(primer, message, max_attempts=max_attempts)
                    elapsed = time.time() - start_time

                    results['execution_times'].append(elapsed)

                    if result.success:
                        results['successes'] += 1
                        results['by_message_size'][size_category]['success'] += 1

                        # Collect backtracking stats
                        if result.encoding_stats:
                            backtrack_count = result.encoding_stats.get('backtrack_count', 0)
                            results['backtracks'].append(backtrack_count)
                    else:
                        results['failures'] += 1
                        results['failed_primers'].append((primer, message))

        return results

    def test_all_profiles_statistics(self):
        """Generate comprehensive statistics for all profiles."""
        profiles = ['sequence_only', 'relaxed', 'pcr_friendly', 'strict']
        all_results = {}

        print(f"\n{'='*80}")
        print(f"RUNNING TESTS FOR ALL PROFILES")
        print(f"{'='*80}\n")

        for profile in profiles:
            print(f"Testing profile: {profile}...", end=' ', flush=True)
            start = time.time()
            results = self._run_profile_test(profile, max_attempts=10)
            elapsed = time.time() - start
            print(f"Done ({elapsed:.1f}s)")
            all_results[profile] = results

        # Print results
        self._print_statistics(all_results)

        # Store for test assertions
        self.all_results = all_results

        # Basic assertions - sequence_only should have highest success rate
        seq_only_rate = all_results['sequence_only']['successes'] / all_results['sequence_only']['total_tests']
        self.assertGreater(seq_only_rate, 0.8,
            f"sequence_only should have >80% success rate, got {seq_only_rate:.1%}")

    def _print_statistics(self, all_results: Dict[str, Dict[str, Any]]):
        """Print formatted statistics table."""
        print(f"\n{'='*80}")
        print(f"RESULTS SUMMARY")
        print(f"{'='*80}\n")

        # Header
        print(f"{'Profile':<15} {'Success Rate':<15} {'Backtracks (avg/med/max)':<30} {'Time (avg)':<15}")
        print(f"{'-'*80}")

        for profile, results in all_results.items():
            total = results['total_tests']
            successes = results['successes']
            success_rate = successes / total if total > 0 else 0

            # Backtrack statistics
            if results['backtracks']:
                bt_avg = statistics.mean(results['backtracks'])
                bt_median = statistics.median(results['backtracks'])
                bt_max = max(results['backtracks'])
                backtrack_str = f"{bt_avg:.1f} / {bt_median:.0f} / {bt_max}"
            else:
                backtrack_str = "N/A"

            # Time statistics
            if results['execution_times']:
                time_avg = statistics.mean(results['execution_times']) * 1000  # ms
                time_str = f"{time_avg:.1f}ms"
            else:
                time_str = "N/A"

            print(f"{profile:<15} {successes}/{total} ({success_rate:.1%})  {backtrack_str:<30} {time_str:<15}")

        # Detailed breakdown by message size
        print(f"\n{'='*80}")
        print(f"SUCCESS RATE BY MESSAGE SIZE")
        print(f"{'='*80}\n")

        for size_category in ['short', 'medium', 'long']:
            print(f"\n{size_category.upper()} messages:")
            print(f"{'Profile':<15} {'Success Rate':<20}")
            print(f"{'-'*40}")

            for profile, results in all_results.items():
                size_data = results['by_message_size'][size_category]
                total = size_data['total']
                success = size_data['success']
                rate = success / total if total > 0 else 0
                print(f"{profile:<15} {success}/{total} ({rate:.1%})")

        # Failed primers analysis
        print(f"\n{'='*80}")
        print(f"FAILURE ANALYSIS")
        print(f"{'='*80}\n")

        for profile, results in all_results.items():
            if results['failed_primers']:
                print(f"\n{profile}: {len(results['failed_primers'])} failures")
                # Show first 3 failures
                for i, (primer, message) in enumerate(results['failed_primers'][:3]):
                    print(f"  {i+1}. Primer: {primer[:20]}... Message: '{message}'")
                if len(results['failed_primers']) > 3:
                    print(f"  ... and {len(results['failed_primers']) - 3} more")
            else:
                print(f"\n{profile}: No failures ✓")

        print(f"\n{'='*80}\n")

    def test_sequence_only_high_success_rate(self):
        """Test that sequence_only profile has high success rate."""
        config = DNAEncoderConfig(validation_profile='sequence_only')
        encoder = DNAEncoder(config)

        successes = 0
        total = 30  # Test 30 random cases

        for _ in range(total):
            primer = random.choice(self.test_primers)
            message = random.choice(self.test_messages['short'])
            result = encoder.encode_text(primer, message, max_attempts=5)
            if result.success:
                successes += 1

        success_rate = successes / total
        self.assertGreater(success_rate, 0.85,
            f"sequence_only should have >85% success rate, got {success_rate:.1%}")

    def test_backtracking_reduces_with_attempts(self):
        """Test that more attempts don't necessarily mean more backtracks per attempt."""
        primer = random.choice(self.test_primers)
        message = "Test"

        config = DNAEncoderConfig(validation_profile='pcr_friendly')
        encoder = DNAEncoder(config)

        # Try with more attempts
        result = encoder.encode_text(primer, message, max_attempts=20)

        # If successful, backtracking should be reasonable
        if result.success and result.encoding_stats:
            backtracks = result.encoding_stats.get('backtrack_count', 0)
            # Reasonable backtrack count
            self.assertLess(backtracks, 1000,
                f"Backtrack count ({backtracks}) seems excessive")


class TestProfileComparison(unittest.TestCase):
    """Compare profiles on identical test cases."""

    @classmethod
    def setUpClass(cls):
        """Load primers."""
        primers_file = Path(__file__).parent.parent / 'primers.txt'
        with open(primers_file) as f:
            all_primers = [line.strip() for line in f if line.strip()]

        random.seed(42)
        cls.test_primers = random.sample(all_primers, 10)

    def test_profile_ordering(self):
        """Test that profiles follow expected success rate ordering."""
        primer = self.test_primers[0]
        message = "Test"

        profiles = ['sequence_only', 'relaxed', 'pcr_friendly', 'strict']
        success_rates = {}

        for profile in profiles:
            config = DNAEncoderConfig(validation_profile=profile)
            encoder = DNAEncoder(config)

            successes = 0
            attempts = 10

            for _ in range(attempts):
                result = encoder.encode_text(primer, message, max_attempts=10)
                if result.success:
                    successes += 1

            success_rates[profile] = successes / attempts

        print(f"\n{'='*60}")
        print(f"PROFILE COMPARISON (same test case, 10 runs)")
        print(f"{'='*60}")
        for profile in profiles:
            rate = success_rates[profile]
            print(f"  {profile:<15}: {rate:.1%}")
        print(f"{'='*60}\n")

        # sequence_only should be most permissive
        self.assertGreaterEqual(success_rates['sequence_only'],
                                success_rates['strict'],
                                "sequence_only should have ≥ success rate than strict")


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)
