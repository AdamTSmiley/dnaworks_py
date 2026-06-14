"""
DNA sequence handling: validation, reverse complement, and input parsing.
"""

import sys
from pathlib import Path

VALID_BASES = set("ACGT")

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def validate(seq: str) -> str:
    """
    Validate and clean a DNA sequence. Strips whitespace, uppercases,
    and errors with non-ACGT characters.

    Parameters
    seq : str
        Raw DNA sequence.

    Returns
    str
        Cleaned, uppercase DNA sequence.

    Raises
    SystemExit
        If the sequence contains non-ACGT characters.
    """
    seq = seq.upper().replace(" ", "").replace("\n", "").replace("\r", "")

    bad = set(seq) - VALID_BASES
    if bad:
        bad_str = ", ".join(sorted(bad))
        # Find first bad position for a helpful error message
        for i, base in enumerate(seq):
            if base not in VALID_BASES:
                sys.exit(
                    f"ERROR: Invalid character(s) in sequence: {bad_str}\n"
                    f"  First occurrence: position {i + 1} ('{base}')\n"
                    f"  Only A, C, G, T are allowed."
                )

    if len(seq) == 0:
        sys.exit("ERROR: Empty sequence.")

    return seq


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def read_fasta(path: str | Path) -> str:
    """
    Read a single sequence from a FASTA file. If the file contains
    multiple records, they are concatenated (with a warning).

    Parameters
    path : str or Path
        Path to the FASTA file.

    Returns
    str
        Validated DNA sequence.
    """
    path = Path(path)
    if not path.exists():
        sys.exit(f"ERROR: File not found: {path}")

    records = []
    current = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    records.append("".join(current))
                    current = []
            elif line:
                current.append(line)

    if current:
        records.append("".join(current))

    if not records:
        sys.exit(f"ERROR: No sequence found in {path}")

    if len(records) > 1:
        print(
            f"WARNING: {path} contains {len(records)} records. "
            f"Concatenating into a single sequence.",
            file=sys.stderr,
        )

    return validate("".join(records))


def read_input(source: str) -> str:
    """
    Read a DNA sequence from either a file path or a raw string.

    If source looks like a file path (ends in .fa, .fasta, .fna, .txt,
    or the file exists), reads it as FASTA. Otherwise treats it as a
    raw sequence string.

    Parameters
    source : str
        File path or raw DNA sequence.

    Returns
    str
        Validated DNA sequence.
    """
    fasta_extensions = {".fa", ".fasta", ".fna", ".txt"}
    p = Path(source)

    if p.suffix.lower() in fasta_extensions or p.exists():
        return read_fasta(p)

    return validate(source)
