#!/usr/bin/env python3
"""Generate DNAWORKS.inp files from FASTA files for Fortran benchmarking."""

import os
import sys
from pathlib import Path


def fasta_to_seq(path):
    """Read a single sequence from a FASTA file."""
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def write_dnaworks_inp(name, seq, outdir, tm=62, tol=2, oligo_len=60):
    """Write a DNAWORKS.inp file for DNA-only mode."""
    gene_dir = Path(outdir) / name
    gene_dir.mkdir(parents=True, exist_ok=True)

    with open(gene_dir / "DNAWORKS.inp", "w") as f:
        f.write(f"title  {name} oligo design (DNA-only benchmark)\n\n")
        f.write(f"melting {tm}\n")
        f.write(f"tolerance {tol}\n")
        f.write(f"length {oligo_len}\n")
        f.write("concentration oligo 2.5e-7  sodium 5e-2  magnesium 2e-3\n\n")
        f.write("nucleotide\n")

        # Write sequence in 60-char lines, indented by one space
        for i in range(0, len(seq), 60):
            f.write(f" {seq[i:i+60]}\n")

        f.write("//\n")

    return gene_dir


def main():
    fasta_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("example_fastas")
    outdir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("fortran_inputs")

    fastas = sorted(fasta_dir.glob("*.fasta"))
    if not fastas:
        print(f"No .fasta files found in {fasta_dir}")
        sys.exit(1)

    print(f"Generating DNAWORKS.inp files for {len(fastas)} genes:")
    for fa in fastas:
        name = fa.stem
        seq = fasta_to_seq(fa)
        gene_dir = write_dnaworks_inp(name, seq, outdir)
        gc = 100 * sum(1 for b in seq if b in "GC") / len(seq)
        print(f"  {name:12s}  {len(seq):5d} bp  {gc:5.1f}% GC  -> {gene_dir}/DNAWORKS.inp")

    print(f"\nDone. To run Fortran on all genes:")
    print(f"  cd <gene_dir> && /path/to/dnaworks")


if __name__ == "__main__":
    main()
