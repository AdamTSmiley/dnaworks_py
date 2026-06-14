#!/usr/bin/env python3
"""
Parse Fortran DNAWorks LOGFILE.txt output into structured JSON.

Robustly extracts the oligo listing section by anchoring to the
"oligonucleotides need to be synthesized" marker, then reconstructs
overlap positions by mapping oligos back to the gene sequence.
"""

import json
import re
import sys
from pathlib import Path


def parse_logfile(logfile_path: str, fasta_path: str) -> dict:
    with open(logfile_path) as f:
        lines = f.readlines()

    gene_seq = _read_fasta(fasta_path)
    oligos = _parse_oligos(lines)
    scores = _parse_scores(lines)
    summary = _parse_summary(lines)
    overlaps = _reconstruct_overlaps(oligos, gene_seq)

    return {
        "gene_length": len(gene_seq),
        "num_oligos": len(oligos),
        "num_overlaps": len(overlaps),
        "oligos": oligos,
        "overlaps": overlaps,
        "scores": scores,
        "summary": summary,
    }


def _read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def _parse_oligos(lines: list[str]) -> list[dict]:
    """
    Parse the oligo listing by finding the section marker first.

    The section looks like:
        N oligonucleotides need to be synthesized
     ================================================================
      1 ATGGTG...CTGGAT 60
      2 CGATC...TGCAC 60
      ...
    """
    oligos = []

    # Find the marker line
    marker_idx = None
    for i, line in enumerate(lines):
        if "oligonucleotides need to be synthesized" in line:
            marker_idx = i
            break

    if marker_idx is None:
        return oligos

    # Skip the bar line (===...) after the marker
    start_idx = marker_idx + 1
    while start_idx < len(lines):
        if lines[start_idx].strip().startswith("=") or lines[start_idx].strip().startswith("-"):
            start_idx += 1
            break
        start_idx += 1

    # Parse oligo lines: format is "  N SEQUENCE LEN"
    oligo_re = re.compile(r"^\s*(\d+)\s+([ACGTacgt]+)\s+(\d+)\s*$")

    for i in range(start_idx, len(lines)):
        line = lines[i]

        # Stop at blank line or next section
        stripped = line.strip()
        if not stripped:
            break
        if stripped.startswith("=") or stripped.startswith("-"):
            break

        match = oligo_re.match(line)
        if match:
            idx = int(match.group(1))
            seq = match.group(2).upper()
            length = int(match.group(3))

            # Sanity check: length should match sequence and be reasonable
            if length == len(seq) and length < 500:
                oligos.append({
                    "index": idx,
                    "sequence": seq,
                    "length": length,
                })

    return oligos


def _parse_scores(lines: list[str]) -> dict:
    scores = {}
    text = "".join(lines)
    patterns = {
        "codon": r"total codon usage score\s*\.+\s+([\d.]+)",
        "length": r"total length score\s*\.+\s+([\d.]+)",
        "tm": r"total melting temperature score\s*\.+\s+([\d.]+)",
        "repeat": r"total repeat score\s*\.+\s+([\d.]+)",
        "pattern": r"total pattern score\s*\.+\s+([\d.]+)",
        "misprime": r"total mispriming score\s*\.+\s+([\d.]+)",
        "at": r"total AT content score\s*\.+\s+([\d.]+)",
        "gc": r"total GC content score\s*\.+\s+([\d.]+)",
        "overall": r"OVERALL score\s*\.+\s+([\d.]+)",
    }
    for key, pat in patterns.items():
        match = re.search(pat, text, re.IGNORECASE)
        if match:
            scores[key] = float(match.group(1))
    return scores


def _parse_summary(lines: list[str]) -> dict:
    text = "".join(lines)
    summary = {}
    pattern = re.compile(
        r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+\|\s+"
        r"([\d.]+)\s+([\d.]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)",
        re.MULTILINE,
    )
    match = pattern.search(text)
    if match:
        summary = {
            "trial": int(match.group(1)),
            "target_tm": int(match.group(2)),
            "target_length": int(match.group(3)),
            "overall_score": float(match.group(4)),
            "tm_range": float(match.group(5)),
            "shortest_oligo": int(match.group(6)),
            "longest_oligo": int(match.group(7)),
            "num_oligos": int(match.group(8)),
            "num_repeats": int(match.group(9)),
            "num_misprimes": int(match.group(10)),
        }
    return summary


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def _reconstruct_overlaps(oligos: list[dict], gene_seq: str) -> list[dict]:
    """Reconstruct overlap positions by mapping oligos to the gene sequence."""
    if not oligos:
        return []

    # Map each oligo to its position on the gene
    positions = []
    for oligo in oligos:
        seq = oligo["sequence"]

        # Try sense
        pos = gene_seq.find(seq)
        if pos >= 0:
            positions.append({"start": pos, "end": pos + len(seq) - 1, "strand": "sense"})
            continue

        # Try antisense (RC)
        rc = _reverse_complement(seq)
        pos = gene_seq.find(rc)
        if pos >= 0:
            positions.append({"start": pos, "end": pos + len(rc) - 1, "strand": "antisense"})
            continue

        # Couldn't map — skip
        positions.append(None)

    # Find overlaps between adjacent oligos
    overlaps = []
    for i in range(len(positions) - 1):
        if positions[i] is None or positions[i + 1] is None:
            continue

        p1 = positions[i]
        p2 = positions[i + 1]

        ov_start = max(p1["start"], p2["start"])
        ov_end = min(p1["end"], p2["end"])

        if ov_start <= ov_end:
            ov_seq = gene_seq[ov_start : ov_end + 1]
            overlaps.append({
                "index": i + 1,
                "start": ov_start + 1,  # 1-indexed
                "end": ov_end + 1,
                "length": ov_end - ov_start + 1,
                "sequence": ov_seq,
            })

    return overlaps


def main():
    fortran_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("fortran_outputs")
    fasta_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("example_fastas")
    outdir = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("fortran_outputs")

    logfiles = sorted(fortran_dir.glob("*.txt"))
    if not logfiles:
        print(f"No .txt files found in {fortran_dir}")
        sys.exit(1)

    print(f"Parsing {len(logfiles)} Fortran output files:")
    for lf in logfiles:
        name = lf.stem
        fasta = fasta_dir / f"{name}.fasta"
        if not fasta.exists():
            print(f"  {name}: FASTA not found, skipping")
            continue

        result = parse_logfile(str(lf), str(fasta))
        out_json = outdir / f"{name}.json"
        with open(out_json, "w") as f:
            json.dump(result, f, indent=2)

        # Validate
        expected_overlaps = result["num_oligos"] - 1
        actual_overlaps = result["num_overlaps"]
        status = "OK" if actual_overlaps == expected_overlaps else f"WARN: expected {expected_overlaps} overlaps"

        print(f"  {name}: {result['num_oligos']} oligos, "
              f"{result['num_overlaps']} overlaps  [{status}]")

    print("Done.")


if __name__ == "__main__":
    main()
