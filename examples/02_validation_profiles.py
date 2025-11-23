"""
Using different validation profiles.

This example demonstrates the different quality profiles available:
- sequence_only: Basic sequence checks only
- relaxed: Relaxed thermodynamic constraints
- pcr_friendly: Moderate constraints for PCR
- strict: Strict quality requirements
"""

from dna_encoder import DNAEncoder, DNAEncoderConfig

# Test primer and message
primer = "CATCTATCCCTTCGAACGAC"
message = "Test"

# Try different profiles
profiles = ['sequence_only', 'relaxed', 'pcr_friendly', 'strict']

print(f"Encoding '{message}' with different profiles:\n")

for profile_name in profiles:
    config = DNAEncoderConfig(validation_profile=profile_name)
    encoder = DNAEncoder(config)

    result = encoder.encode_text(primer, message, max_attempts=10)

    status = "✓" if result.success else "✗"
    print(f"{status} {profile_name:15s}: ", end="")

    if result.success:
        print(f"Success (length: {len(result.sequence)} nt)")
    else:
        attempts = result.encoding_stats.get('attempts', 10) if result.encoding_stats else 10
        print(f"Failed after {attempts} attempts")

print("\nNote: Stricter profiles may fail more often due to biochemical constraints.")
