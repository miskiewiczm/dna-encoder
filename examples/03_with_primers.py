"""
Using validated PCR primers.

This example shows how to use real PCR-validated primers from primers.txt
for encoding data into DNA sequences.
"""

from pathlib import Path
from dna_encoder import DNAEncoder, DNADecoder, DNAEncoderConfig
import random

# Load primers from file
primers_file = Path(__file__).parent.parent / 'primers.txt'

if primers_file.exists():
    with open(primers_file) as f:
        primers = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(primers)} validated primers\n")
else:
    print("primers.txt not found, using default primer")
    primers = ["CATCTATCCCTTCGAACGAC"]

# Use a random primer
primer = random.choice(primers)
print(f"Selected primer: {primer}")

# Configure encoder with PCR-friendly profile
config = DNAEncoderConfig(validation_profile='pcr_friendly')
encoder = DNAEncoder(config)
decoder = DNADecoder()

# Encode message
message = "DNA Storage"
print(f"Encoding: '{message}'")

result = encoder.encode_text(primer, message, max_attempts=20)

if result.success:
    print(f"\n✓ Encoding successful!")
    print(f"  Final sequence: {result.sequence}")
    print(f"  Length: {len(result.sequence)} nt")

    # Show encoding stats if available
    if result.encoding_stats:
        print(f"  Backtracks: {result.encoding_stats.get('backtrack_count', 0)}")
        print(f"  Total choices: {result.encoding_stats.get('total_choices', 0)}")

    # Decode to verify
    decoded = decoder.decode_sequence(result.sequence, primer)
    print(f"\n  Decoded: '{decoded}'")
    print(f"  Verified: {decoded == message}")

else:
    print(f"\n✗ Encoding failed")
    print(f"  Try using 'sequence_only' or 'relaxed' profile for better success rate")
