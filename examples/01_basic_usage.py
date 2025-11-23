"""
Basic usage example for dna-encoder.

This example shows how to encode text into DNA and decode it back.
"""

from dna_encoder import DNAEncoder, DNADecoder

# Initialize encoder and decoder
encoder = DNAEncoder()
decoder = DNADecoder()

# Text to encode
message = "Hello, DNA!"

# Starting primer sequence (20 nt)
primer = "CATCTATCCCTTCGAACGAC"

# Encode the message
print(f"Encoding: '{message}'")
result = encoder.encode_text(primer, message)

if result.success:
    print(f"✓ Success!")
    print(f"  Sequence: {result.sequence}")
    print(f"  Length: {len(result.sequence)} nt")
    if result.encoding_stats:
        print(f"  Backtracks: {result.encoding_stats.get('backtrack_count', 0)}")

    # Decode back to verify
    decoded = decoder.decode_sequence(result.sequence, primer)
    print(f"  Decoded: '{decoded}'")
    print(f"  Match: {decoded == message}")
else:
    print(f"✗ Encoding failed: {result.error_message}")
