"""
Encoding binary data.

This example shows how to encode arbitrary binary data into DNA sequences.
"""

from dna_encoder import DNAEncoder, DNADecoder, DNAEncoderConfig

# Initialize with sequence_only profile for higher success rate
config = DNAEncoderConfig(validation_profile='sequence_only')
encoder = DNAEncoder(config)
decoder = DNADecoder()

primer = "CATCTATCCCTTCGAACGAC"

# Example 1: Encode a byte directly as bits
print("Example 1: Single byte")
byte_value = 65  # ASCII 'A'
bits = format(byte_value, '08b')  # Convert to 8-bit string

print(f"  Byte: {byte_value} ('{chr(byte_value)}')")
print(f"  Bits: {bits}")

result = encoder.encode_bits(primer, bits)
if result.success:
    print(f"  ✓ Encoded: {result.sequence}")
    decoded = decoder.decode_sequence(result.sequence, primer)
    print(f"  ✓ Decoded: '{decoded}' (matches: {decoded == chr(byte_value)})")
else:
    print(f"  ✗ Failed")

print()

# Example 2: Encode multiple bytes
print("Example 2: Multiple bytes")
message = "DNA"
bits = ''.join(format(ord(c), '08b') for c in message)

print(f"  Message: '{message}'")
print(f"  Bits: {bits} ({len(bits)} bits)")

result = encoder.encode_bits(primer, bits)
if result.success:
    print(f"  ✓ Encoded: {result.sequence[:50]}... ({len(result.sequence)} nt)")
    decoded = decoder.decode_sequence(result.sequence, primer)
    print(f"  ✓ Decoded: '{decoded}' (matches: {decoded == message})")
else:
    print(f"  ✗ Failed")

print()

# Example 3: Using encode_text (higher level)
print("Example 3: Using encode_text (recommended)")
message = "Hello"
print(f"  Message: '{message}'")

result = encoder.encode_text(primer, message)
if result.success:
    print(f"  ✓ Encoded: {result.sequence[:50]}...")
    decoded = decoder.decode_sequence(result.sequence, primer)
    print(f"  ✓ Decoded: '{decoded}' (matches: {decoded == message})")
