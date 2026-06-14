#!/usr/bin/env python3
"""
dnaworks_py: Design overlapping oligonucleotides for PCR-based gene assembly.

Usage:
  python cli.py ATGAAA...TAA                    # raw sequence
  python cli.py my_gene.fa                      # FASTA file
  python cli.py my_gene.fa --tm 60 --length 50  # custom params
"""

import argparse
import sys
import time

from sequence import read_input
from scoring import ScoringConfig
from overlaps import design_oligos


def main():
    parser = argparse.ArgumentParser(
        prog="dnaworks_py",
        description="Design overlapping oligos for PCR-based gene assembly.",
    )

    parser.add_argument(
        "sequence",
        help="DNA sequence (raw string) or path to FASTA file",
    )
    parser.add_argument(
        "--tm", type=float, default=62.0,
        help="Target overlap Tm in °C (default: 62.0)",
    )
    parser.add_argument(
        "--tm-tolerance", type=float, default=2.0,
        help="Tm tolerance in °C (default: 2.0)",
    )
    parser.add_argument(
        "--length", type=int, default=60,
        help="Target oligo length in nt (default: 60)",
    )
    parser.add_argument(
        "--method", choices=["breslauer86", "santalucia97"],
        default="breslauer86",
        help="Tm calculation method (default: breslauer86)",
    )
    parser.add_argument(
        "--na", type=float, default=50.0,
        help="Sodium concentration in mM (default: 50)",
    )
    parser.add_argument(
        "--mg", type=float, default=0.0,
        help="Magnesium concentration in mM (default: 0)",
    )
    parser.add_argument(
        "--dnac", type=float, default=250.0,
        help="Total primer concentration in nM (default: 250)",
    )
    parser.add_argument(
        "--nogaps", action="store_true",
        help="No gaps between overlaps (minimize oligo length)",
    )
    parser.add_argument(
        "--shifts", type=int, default=1000,
        help="Number of starting offsets to try (default: 1000)",
    )
    parser.add_argument(
        "--json", action="store_true",
        help="Output in JSON format",
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Suppress progress messages",
    )

    args = parser.parse_args()

    # Read and validate sequence
    seq = read_input(args.sequence)

    if not args.quiet:
        print(f"Sequence: {len(seq)} bp", file=sys.stderr)
        gc = 100 * sum(1 for b in seq if b in "GC") / len(seq)
        print(f"GC content: {gc:.1f}%", file=sys.stderr)
        print(f"Target Tm: {args.tm}°C (±{args.tm_tolerance}°C)", file=sys.stderr)
        print(f"Target oligo length: {args.length} nt", file=sys.stderr)
        print(f"Tm method: {args.method}", file=sys.stderr)
        print(f"Designing...", file=sys.stderr)

    config = ScoringConfig(
        target_tm=args.tm,
        tm_tolerance=args.tm_tolerance,
        target_length=args.length,
        tm_method=args.method,
        na=args.na,
        mg=args.mg,
        dnac=args.dnac,
    )

    start = time.time()

    result = design_oligos(
        seq, config,
        oligo_length=args.length,
        nogaps=args.nogaps,
        max_shifts=args.shifts,
    )

    elapsed = time.time() - start

    if not args.quiet:
        print(f"Done in {elapsed:.1f}s\n", file=sys.stderr)

    if args.json:
        _print_json(result, seq, args)
    else:
        _print_text(result, seq, args, elapsed)


def _print_text(result, seq, args, elapsed):
    """Print human-readable output."""
    n = len(seq)
    gc = 100 * sum(1 for b in seq if b in "GC") / n

    print("=" * 70)
    print("dnaworks_py — Oligonucleotide Design Report")
    print("=" * 70)
    print()
    print(f"Sequence length:    {n} bp")
    print(f"GC content:         {gc:.1f}%")
    print(f"Tm method:          {args.method}")
    print(f"Target Tm:          {args.tm}°C (±{args.tm_tolerance}°C)")
    print(f"Target oligo len:   {args.length} nt")
    print(f"[Na+]:              {args.na} mM")
    print(f"[Mg2+]:             {args.mg} mM")
    print(f"Primer conc:        {args.dnac} nM")
    print(f"Design time:        {elapsed:.1f}s")
    print()

    # Scores
    s = result.scores
    print(f"Scores:")
    print(f"  Tm deviation:     {s.tm:.2f}")
    print(f"  Mispriming:       {s.misprime:.2f}")
    print(f"  Oligo length:     {s.length:.2f}")
    print(f"  Repeats:          {s.repeat:.2f}")
    print(f"  GC stretches:     {s.gc:.2f}")
    print(f"  AT stretches:     {s.at:.2f}")
    print(f"  TOTAL:            {s.total:.2f}")
    print()

    # Overlaps
    print(f"Overlaps ({result.num_overlaps}):")
    print(f"  {'#':>3s}  {'Start':>5s}  {'End':>5s}  {'Len':>3s}  {'Tm':>6s}")
    print(f"  {'---':>3s}  {'-----':>5s}  {'-----':>5s}  {'---':>3s}  {'------':>6s}")
    for i, (ov, t) in enumerate(zip(result.overlaps, result.overlap_tms)):
        print(f"  {i+1:>3d}  {ov[0]+1:>5d}  {ov[1]+1:>5d}  {ov[1]-ov[0]+1:>3d}  {t:>6.1f}")
    print()

    # Oligos
    print(f"Oligos ({result.num_oligos}):")
    print(f"  {'#':>3s}  {'Strand':>9s}  {'Len':>3s}  {'Pos':>11s}  {'Sequence'}")
    print(f"  {'---':>3s}  {'---------':>9s}  {'---':>3s}  {'-----------':>11s}  {'--------'}")
    for o in result.oligos:
        print(
            f"  {o['index']:>3d}  {o['strand']:>9s}  {o['length']:>3d}  "
            f"{o['start']+1:>5d}-{o['end']+1:<5d}  {o['sequence']}"
        )
    print()

    # Summary for ordering
    print("Oligo Order List:")
    for o in result.oligos:
        print(f"  oligo_{o['index']:02d}  {o['sequence']}")


def _print_json(result, seq, args):
    """Print JSON output for programmatic use."""
    import json

    output = {
        "sequence_length": len(seq),
        "gc_content": round(100 * sum(1 for b in seq if b in "GC") / len(seq), 1),
        "parameters": {
            "target_tm": args.tm,
            "tm_tolerance": args.tm_tolerance,
            "oligo_length": args.length,
            "tm_method": args.method,
            "na_mm": args.na,
            "mg_mm": args.mg,
            "dnac_nm": args.dnac,
        },
        "scores": {
            "tm": round(result.scores.tm, 3),
            "misprime": round(result.scores.misprime, 3),
            "length": round(result.scores.length, 3),
            "repeat": round(result.scores.repeat, 3),
            "gc": round(result.scores.gc, 3),
            "at": round(result.scores.at, 3),
            "total": round(result.scores.total, 3),
        },
        "num_overlaps": result.num_overlaps,
        "num_oligos": result.num_oligos,
        "overlaps": [
            {
                "index": i + 1,
                "start": ov[0] + 1,
                "end": ov[1] + 1,
                "length": ov[1] - ov[0] + 1,
                "tm": round(t, 1),
            }
            for i, (ov, t) in enumerate(zip(result.overlaps, result.overlap_tms))
        ],
        "oligos": [
            {
                "index": o["index"],
                "strand": o["strand"],
                "length": o["length"],
                "start": o["start"] + 1,
                "end": o["end"] + 1,
                "sequence": o["sequence"],
            }
            for o in result.oligos
        ],
    }

    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
