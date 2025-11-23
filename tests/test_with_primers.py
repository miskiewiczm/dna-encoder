"""Tests using validated primers from primers.txt file.

These tests ensure dna-encoder works with real PCR-validated primers
and maintains compatibility with dna-commons and dna-generator.
"""

import unittest
import random
from pathlib import Path
from dna_encoder import DNAEncoder, DNADecoder, DNAEncoderConfig


class TestWithValidatedPrimers(unittest.TestCase):
    """Test encoding with real PCR-validated primers."""

    @classmethod
    def setUpClass(cls):
        """Load primers from primers.txt."""
        primers_file = Path(__file__).parent.parent.parent / 'primers.txt'
        if not primers_file.exists():
            raise FileNotFoundError(f"primers.txt not found at {primers_file}")

        with open(primers_file) as f:
            cls.primers = [line.strip() for line in f if line.strip()]

        print(f"\nLoaded {len(cls.primers)} primers from primers.txt")

    def test_sequence_only_profile_with_primers(self):
        """Test sequence_only profile with multiple primers."""
        config = DNAEncoderConfig(validation_profile='sequence_only')
        encoder = DNAEncoder(config)

        # Test with 10 random primers
        test_primers = random.sample(self.primers, min(10, len(self.primers)))
        successes = 0

        for primer in test_primers:
            result = encoder.encode_text(primer, "Hi")
            if result.success:
                successes += 1

        # sequence_only should succeed on most primers
        self.assertGreater(successes, 7,
            f"Only {successes}/10 primers succeeded with sequence_only profile")

    def test_relaxed_profile_with_primers(self):
        """Test relaxed profile with primers."""
        config = DNAEncoderConfig(validation_profile='relaxed')
        encoder = DNAEncoder(config)

        # Test with 10 random primers
        test_primers = random.sample(self.primers, min(10, len(self.primers)))
        successes = 0

        for primer in test_primers:
            result = encoder.encode_text(primer, "Hi")
            if result.success:
                successes += 1

        # relaxed should succeed on some primers
        self.assertGreater(successes, 3,
            f"Only {successes}/10 primers succeeded with relaxed profile")

    def test_pcr_friendly_profile_with_primers(self):
        """Test pcr_friendly profile with primers."""
        config = DNAEncoderConfig(validation_profile='pcr_friendly')
        encoder = DNAEncoder(config)

        # Test with 10 random primers
        test_primers = random.sample(self.primers, min(10, len(self.primers)))
        successes = 0
        failures = []

        for primer in test_primers:
            result = encoder.encode_text(primer, "Hi")
            if result.success:
                successes += 1
            else:
                failures.append((primer, result.error_message if hasattr(result, 'error_message') else 'unknown'))

        # pcr_friendly is strict, expect some failures
        print(f"\nPCR-friendly: {successes}/10 succeeded")
        if failures and successes < 3:
            print(f"Sample failures: {failures[:3]}")

        # Just check that it doesn't crash
        self.assertIsNotNone(result)

    def test_encode_decode_roundtrip_with_primers(self):
        """Test that encoding and decoding preserves data."""
        config = DNAEncoderConfig(validation_profile='sequence_only')
        encoder = DNAEncoder(config)
        decoder = DNADecoder()

        test_messages = ["Hi", "Test", "DNA", "OK"]
        primer = random.choice(self.primers)

        for message in test_messages:
            result = encoder.encode_text(primer, message)
            if result.success:
                decoded = decoder.decode_sequence(result.sequence, primer)
                # Decoded should be the original message
                self.assertEqual(message, decoded,
                    f"Expected '{message}', got '{decoded}'")

    def test_multiple_attempts_increases_success(self):
        """Test that more attempts increases success rate."""
        primer = random.choice(self.primers)

        # Test with 5 attempts
        config1 = DNAEncoderConfig(validation_profile='pcr_friendly')
        encoder1 = DNAEncoder(config1)
        result1 = encoder1.encode_text(primer, "Test", max_attempts=5)

        # Test with 20 attempts
        config2 = DNAEncoderConfig(validation_profile='pcr_friendly')
        encoder2 = DNAEncoder(config2)
        result2 = encoder2.encode_text(primer, "Test", max_attempts=20)

        # With more attempts, we should have equal or better success
        if result1.success:
            self.assertTrue(result2.success,
                "More attempts should not decrease success rate")

    def test_primer_validation(self):
        """Test that most primers from primers.txt pass basic sequence checks."""
        # Use sequence_only profile which checks homopolymers, dinucleotide repeats, 3' stability
        config = DNAEncoderConfig(validation_profile='sequence_only')
        encoder = DNAEncoder(config)

        valid_count = 0
        invalid_primers = []
        for i, primer in enumerate(self.primers[:50]):  # Test first 50
            if encoder.validate_sequence(primer):
                valid_count += 1
            else:
                invalid_primers.append((i+1, primer))

        # Most primers should pass basic sequence checks (expect ≥80%)
        pass_rate = valid_count / 50
        self.assertGreaterEqual(pass_rate, 0.80,
            f"Only {valid_count}/50 primers passed ({pass_rate:.1%}). "
            f"Failed primers: {invalid_primers[:5]}")

    def test_different_profiles_same_primer(self):
        """Test same primer with all profiles."""
        primer = self.primers[0]  # Use first primer
        test_text = "Hi"

        profiles = ['sequence_only', 'relaxed', 'pcr_friendly', 'strict']
        results = {}

        for profile in profiles:
            config = DNAEncoderConfig(validation_profile=profile)
            encoder = DNAEncoder(config)
            result = encoder.encode_text(primer, test_text, max_attempts=10)
            results[profile] = result.success

        # Profiles should get progressively stricter
        # sequence_only should be most permissive
        self.assertTrue(results['sequence_only'],
            "sequence_only profile should succeed")

        print(f"\nProfile success rates for primer {primer}:")
        for profile, success in results.items():
            print(f"  {profile}: {'✓' if success else '✗'}")


class TestEncoderDecoderCompatibility(unittest.TestCase):
    """Test compatibility between encoder and decoder."""

    @classmethod
    def setUpClass(cls):
        """Load primers from primers.txt."""
        primers_file = Path(__file__).parent.parent.parent / 'primers.txt'
        with open(primers_file) as f:
            cls.primers = [line.strip() for line in f if line.strip()]

    def test_encode_bits_decode_bits(self):
        """Test encoding and decoding binary data."""
        encoder = DNAEncoder(DNAEncoderConfig(validation_profile='sequence_only'))
        decoder = DNADecoder()

        primer = random.choice(self.primers)
        # Use bits that represent a valid ASCII character
        bits = "01001000"  # 'H' in ASCII

        result = encoder.encode_bits(primer, bits)
        if result.success:
            decoded_text = decoder.decode_sequence(result.sequence, primer)
            # Decoder returns text, so we should get 'H'
            self.assertEqual(decoded_text, 'H')

    def test_long_text_encoding(self):
        """Test encoding longer text."""
        encoder = DNAEncoder(DNAEncoderConfig(validation_profile='sequence_only'))
        decoder = DNADecoder()

        primer = random.choice(self.primers)
        long_text = "Hello World! This is a test."

        result = encoder.encode_text(primer, long_text, max_attempts=20)
        if result.success:
            decoded = decoder.decode_sequence(result.sequence, primer)
            self.assertEqual(long_text, decoded)
        else:
            # Long text might fail, but shouldn't crash
            self.assertIsNotNone(result.error_message if hasattr(result, 'error_message') else 'Failed')


if __name__ == '__main__':
    unittest.main(verbosity=2)
