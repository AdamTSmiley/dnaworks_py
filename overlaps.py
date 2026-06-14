"""
Overlap placement for oligonucleotide assembly.

Places overlaps across a DNA sequence such that each overlap region has
a melting temperature close to the target. Tries many starting offsets
and keeps the arrangement with the best overall score.

Algorithm:
  1. For each starting shift (0 to ~1000):
     a. Place first overlap: binary search forward (forolap)
     b. Advance by oligo_length to estimate next oligo boundary
     c. Place next overlap: binary search backward (revolap) from estimated end
     d. Repeat until end of sequence
     e. Require odd number of overlaps (needed for assembly PCR)
     f. Score the placement, keep if best
  2. If all shifts give even overlaps, increase oligo length and retry
"""

import random as _random
from dataclasses import dataclass
from tm import tm as calc_tm
from scoring import (
    SequenceAnalysis,
    ScoringConfig,
    Scores,
    analyze_sequence,
    score_placement,
)


@dataclass
class DesignResult:
    """Result of overlap placement optimization."""
    overlaps: list[tuple[int, int]]
    oligos: list[dict]
    scores: Scores
    overlap_tms: list[float]
    num_overlaps: int
    num_oligos: int


def design_oligos(
    seq: str,
    config: ScoringConfig,
    oligo_length: int = 60,
    nogaps: bool = False,
    random_length: bool = False,
    max_shifts: int = 1000,
    seed: int | None = None,
) -> DesignResult:
    """
    Design overlapping oligonucleotides for assembly PCR.

    Parameters
    ----------
    seq : str
        Validated DNA sequence (uppercase ACGT only).
    config : ScoringConfig
        Scoring parameters (target Tm, weights, etc.).
    oligo_length : int
        Target oligo length in nt. Default 60.
    nogaps : bool
        If True, no gaps between overlaps (oligos as short as possible).
    random_length : bool
        If True, randomize oligo lengths for SA optimization.
    max_shifts : int
        Number of starting offsets to try. Default 1000.
    seed : int or None
        Random seed for reproducibility.

    Returns
    -------
    DesignResult
    """
    if seed is not None:
        _random.seed(seed)

    n = len(seq)
    tm_kw = dict(method=config.tm_method, na=config.na, mg=config.mg, dnac=config.dnac)
    target_tm = config.target_tm
    tm_tol = config.tm_tolerance

    # Pre-compute sequence-dependent scores
    analysis = analyze_sequence(seq, config)

    best_overlaps = None
    best_score = float("inf")
    best_scores_obj = None
    working_len = oligo_length

    for retry in range(10000):
        found_odd = False

        for shift in range(max_shifts):
            overlaps = _place_overlaps(
                seq, n, shift, working_len, nogaps, random_length,
                target_tm, tm_tol, tm_kw,
            )

            if not overlaps or len(overlaps) % 2 == 0:
                continue

            found_odd = True
            scores = score_placement(analysis, overlaps, config)

            if scores.total < best_score:
                best_score = scores.total
                best_overlaps = overlaps
                best_scores_obj = scores

        if found_odd and best_overlaps is not None:
            break

        if (retry + 1) % 200 == 0:
            working_len += 1

    if best_overlaps is None:
        raise RuntimeError(
            "Failed to find a valid overlap arrangement. "
            "Try increasing oligo length or adjusting Tm target."
        )

    oligos = _extract_oligos(seq, best_overlaps)

    overlap_tms = [
        calc_tm(seq[s : e + 1], **tm_kw) for s, e in best_overlaps
    ]

    return DesignResult(
        overlaps=best_overlaps,
        oligos=oligos,
        scores=best_scores_obj,
        overlap_tms=overlap_tms,
        num_overlaps=len(best_overlaps),
        num_oligos=len(oligos),
    )


def _place_overlaps(
    seq: str, n: int, shift: int, oligo_len: int,
    nogaps: bool, random_length: bool,
    target_tm: float, tm_tol: float, tm_kw: dict,
) -> list[tuple[int, int]]:
    """
    Place overlaps for a single starting shift.

    Follows the Fortran Generate_Overlaps / Make_Olap logic:
      1. ForOlap for first overlap (search forward for end)
      2. Advance by oligo_length, estimate next overlap region
      3. RevOlap for subsequent overlaps (search backward for start)
      4. Adjust if overlap collides with previous
    """
    overlaps = []

    def get_olength():
        if nogaps:
            return 0
        if random_length:
            return _random.randint(20, oligo_len)
        return oligo_len

    # --- First overlap ---
    first = shift
    if first >= n - 8:
        return []

    last = _forolap(seq, first, n, target_tm, tm_tol, tm_kw)
    if last is None or last >= n:
        return []

    overlaps.append((first, last))

    # Advance: estimate where the first oligo ends
    olength = get_olength()
    # first still points to start of overlap 1; advance by oligo_length
    est_last = first + olength - 1  # estimated end of oligo 1
    est_first = est_last - 7        # next overlap starts ~7 nt before oligo end

    # Don't collide with previous overlap
    if est_first <= overlaps[-1][1]:
        est_first = overlaps[-1][1] + 1
        est_last = est_first + 7

    # Bail if first overlap is already past where an oligo would fit
    if shift > 0 and overlaps[0][1] >= oligo_len:
        return []

    # --- Subsequent overlaps ---
    for _ in range(999):
        olength = get_olength()

        if est_last >= n:
            break

        # RevOlap: search backward from est_last to find overlap start
        new_first = _revolap(seq, est_last, n, target_tm, tm_tol, tm_kw)
        if new_first is None:
            break

        new_last = est_last

        # If RevOlap result collides with previous overlap, use ForOlap instead
        if new_first <= overlaps[-1][1]:
            new_first = overlaps[-1][1] + 1
            new_last = _forolap(seq, new_first, n, target_tm, tm_tol, tm_kw)
            if new_last is None or new_last >= n:
                break

        overlaps.append((new_first, new_last))

        # Advance for next iteration
        # The Fortran updates first=new_first (via pass-by-ref from Make_Olap)
        # then: est_last = first + olength - 1; est_first = est_last - 7
        est_last = new_first + olength - 1
        est_first = est_last - 7

        if est_first <= overlaps[-1][1]:
            est_first = overlaps[-1][1] + 1
            est_last = est_first + 7

    return overlaps


def _forolap(
    seq: str, first: int, n: int,
    target_tm: float, tm_tol: float, tm_kw: dict,
) -> int | None:
    """Binary search forward: find end position for an overlap starting at first."""
    step = 32
    last = first + step
    step //= 2

    while step >= 1:
        if last >= n:
            last -= step
        else:
            ov = seq[first : last + 1]
            if len(ov) <= 7:
                last += step
            else:
                diff = target_tm - calc_tm(ov, **tm_kw)
                if abs(diff) <= tm_tol:
                    break
                elif diff > 0:
                    last += step
                else:
                    last -= step
        step //= 2

    # Clamp
    last = max(first + 7, min(last, n - 1))
    if last >= n:
        return None

    # Final step: compare last vs last+1 (or last-1)
    d1 = abs(target_tm - calc_tm(seq[first : last + 1], **tm_kw))
    if last + 1 < n:
        d2 = abs(target_tm - calc_tm(seq[first : last + 2], **tm_kw))
        if d2 < d1:
            last += 1
    if n - last <= 2:
        last = n - 1
    return last


def _revolap(
    seq: str, last: int, n: int,
    target_tm: float, tm_tol: float, tm_kw: dict,
) -> int | None:
    """Binary search backward: find start position for an overlap ending at last."""
    step = 32
    first = last - step
    step //= 2

    while step >= 1:
        if first < 0:
            first += step
        else:
            ov = seq[first : last + 1]
            if len(ov) <= 7:
                first -= step
            else:
                diff = target_tm - calc_tm(ov, **tm_kw)
                if abs(diff) <= tm_tol:
                    break
                elif diff > 0:
                    first -= step
                else:
                    first += step
        step //= 2

    # Clamp
    first = max(0, min(first, last - 7))
    if first < 0:
        return None

    # Final step: compare first vs first-1
    d1 = abs(target_tm - calc_tm(seq[first : last + 1], **tm_kw))
    if first > 0:
        d2 = abs(target_tm - calc_tm(seq[first - 1 : last + 1], **tm_kw))
        if d2 < d1:
            first -= 1
    if first <= 1:
        first = 0
    return first


def _extract_oligos(seq: str, overlaps: list[tuple[int, int]]) -> list[dict]:
    """
    Extract oligo sequences from overlap positions.

    Oligos alternate sense/antisense:
      Oligo 0 (sense):     seq[0 .. overlap[0].end]
      Oligo 1 (antisense): RC(seq[overlap[0].start .. overlap[1].end])
      Oligo 2 (sense):     seq[overlap[1].start .. overlap[2].end]
      ...
      Last oligo (antisense): RC(seq[overlap[-1].start .. seqlen-1])
    """
    from sequence import reverse_complement

    n_overlaps = len(overlaps)
    oligos = []

    for i in range(n_overlaps + 1):
        if i == 0:
            start, end = 0, overlaps[0][1]
        elif i == n_overlaps:
            start, end = overlaps[-1][0], len(seq) - 1
        else:
            start, end = overlaps[i - 1][0], overlaps[i][1]

        oligo_seq = seq[start : end + 1]
        is_sense = (i % 2 == 0)

        if not is_sense:
            oligo_seq = reverse_complement(oligo_seq)

        oligos.append({
            "index": i + 1,
            "sequence": oligo_seq,
            "length": len(oligo_seq),
            "start": start,
            "end": end,
            "strand": "sense" if is_sense else "antisense",
        })

    return oligos
