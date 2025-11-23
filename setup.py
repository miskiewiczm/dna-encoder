"""Setup configuration for dna-encoder package."""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

setup(
    name="dna-encoder",
    version="0.1.0",
    author="Marek Miskiewicz",
    author_email="miskiewiczm@users.noreply.github.com",
    description="Encode data into DNA sequences with biochemical quality control and backtracking",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/miskiewiczm/dna-encoder",
    packages=find_packages(exclude=["tests", "tests.*", "examples"]),
    package_data={
        'dna_encoder': ['*.json'],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    python_requires=">=3.8",
    install_requires=[
        "dna-commons>=0.1.0",  # Our base library
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "dna-encoder=dna_encoder.__main__:main",
        ],
    },
    include_package_data=True,
)
