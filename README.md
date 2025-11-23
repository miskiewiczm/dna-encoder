# DNA Encoder

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Encode data into DNA sequences with biochemical quality control and backtracking.

## Features

- **Text and binary encoding** into DNA sequences
- **Multiple validation profiles** (sequence_only, relaxed, pcr_friendly, strict)
- **Backtracking algorithm** for high-quality sequence generation
- **PCR-validated primers** support (500+ primers included)
- **Biochemical quality control** (GC content, Tm, hairpins, homodimers)
- **Decode back to original data** with verification

## Installation

```bash
pip install dna-encoder
```

Or install in development mode:

```bash
git clone https://github.com/miskiewiczm/dna-encoder.git
cd dna-encoder
pip install -e .
```

## Quick Start

```python
from dna_encoder import DNAEncoder, DNADecoder

# Initialize
encoder = DNAEncoder()
decoder = DNADecoder()

# Encode message
primer = "CATCTATCCCTTCGAACGAC"
result = encoder.encode_text(primer, "Hello, DNA!")

if result.success:
    print(f"Encoded: {result.sequence}")

    # Decode back
    decoded = decoder.decode_sequence(result.sequence, primer)
    print(f"Decoded: {decoded}")  # "Hello, DNA!"
```

## Validation Profiles

The encoder supports different quality profiles:

| Profile | Description | Use Case |
|---------|-------------|----------|
| `sequence_only` | Basic sequence checks only | Maximum success rate |
| `relaxed` | Relaxed thermodynamic constraints | General purpose |
| `pcr_friendly` | Moderate PCR constraints | PCR amplification |
| `strict` | Strict quality requirements | High-quality applications |

```python
from dna_encoder import DNAEncoderConfig, DNAEncoder

config = DNAEncoderConfig(validation_profile='pcr_friendly')
encoder = DNAEncoder(config)
```

## Using Validated Primers

The package includes 500+ PCR-validated primers in `primers.txt`:

```python
from pathlib import Path

# Load primers
with open('primers.txt') as f:
    primers = [line.strip() for line in f]

# Use any primer
primer = primers[0]
result = encoder.encode_text(primer, "Data")
```

## CLI Usage

```bash
# Encode text
dna-encoder --encode "Hello" --initial ATGCATGC --profile pcr_friendly

# Decode sequence
dna-encoder --decode ATGCATGCAAGGTTCC --initial ATGCATGC

# Show help
dna-encoder --help
```

## Examples

See the `examples/` directory for more usage examples:

- `01_basic_usage.py` - Basic encoding and decoding
- `02_validation_profiles.py` - Using different profiles
- `03_with_primers.py` - Using validated primers
- `04_binary_data.py` - Encoding binary data

## Testing

```bash
# Run tests
pytest tests/

# With coverage
pytest tests/ --cov=dna_encoder --cov-report=term-missing
```

Current test coverage: **50%** (24 tests passing)

## Dependencies

- `dna-commons>=0.1.0` - Core DNA operations library

## Project Structure

```
dna-encoder/
├── dna_encoder/          # Main package
│   ├── encoder.py        # Encoding logic
│   ├── decoder.py        # Decoding logic
│   ├── backtracking_engine.py  # Backtracking algorithm
│   ├── config.py         # Configuration
│   └── default_profiles.json   # Validation profiles
├── tests/                # Test suite
│   └── test_with_primers.py    # Primer-based tests
├── examples/             # Usage examples
└── primers.txt           # 500+ validated primers
```

## License

MIT License

## Contributing

Contributions welcome! Please open an issue or submit a pull request.

## Links

- GitHub: https://github.com/miskiewiczm/dna-encoder
- Issues: https://github.com/miskiewiczm/dna-encoder/issues

## Citation

If you use this software in your research, please cite:

```bibtex
@software{dna_encoder,
  author = {Miskiewicz, Marek},
  title = {DNA Encoder: Encode Data into DNA Sequences with Biochemical Quality Control},
  year = {2025},
  url = {https://github.com/miskiewiczm/dna-encoder},
  version = {0.1.0}
}
```
