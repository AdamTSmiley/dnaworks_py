"""
Scoring functions for oligonucleotide overlap placement.

Scores are split into two categories:
  - Sequence-dependent (computed once): repeats, GC content, AT content
  - Placement-dependent (computed per trial): Tm deviation, misprimes, length

The misprime detection follows DNAWorks logic:
  1. Pre-scan for all position pairs with similar 3' tips and overall homology
  2. For each trial placement, check if any of these pairs coincide with
     an oligo 3' end

Oligo 3' ends exist only at every other overlap junction:
  - In standard assembly, odd overlaps (0-indexed: 0, 2, 4, ...) are where
    a sense oligo's 3' end and an antisense oligo's 3' end both terminate.
  - Even overlaps (0-indexed: 1, 3, 5, ...) are 5' junctions — no priming.
"""

from dataclasses import dataclass, field
from sequence import reverse_complement
from tm import tm


# ============================================================================
# Data structures
# ============================================================================

@dataclass
class PotentialMisprime:
    """A pair of positions that could misprime if an oligo 3' end lands there."""
    pos1: int       # start of window 1 (0-indexed)
    pos2: int       # start of window 2 (0-indexed)
    mtype: int      # 1=direct-sense, 2=inverse-sense, 3=inverse-antisense, 4=direct-antisense


@dataclass
class SequenceAnalysis:
    """Pre-computed, placement-independent analysis of a DNA sequence."""
    seq: str
    rc: str
    length: int

    # Potential misprimes (placement-independent catalog)
    potential_misprimes: list[PotentialMisprime] = field(default_factory=list)

    # Per-position scores (placement-independent)
    repeat_scores: list[float] = field(default_factory=list)
    gc_scores: list[float] = field(default_factory=list)
    at_scores: list[float] = field(default_factory=list)

    # Totals (pre-normalized by 20/len)
    total_repeat: float = 0.0
    total_gc: float = 0.0
    total_at: float = 0.0


@dataclass
class Scores:
    """Complete scores for a specific overlap placement."""
    tm: float = 0.0
    misprime: float = 0.0
    length: float = 0.0
    repeat: float = 0.0
    gc: float = 0.0
    at: float = 0.0
    total: float = 0.0


@dataclass
class ScoringConfig:
    """Parameters controlling the scoring functions."""
    # Tm parameters
    target_tm: float = 62.0
    tm_tolerance: float = 1.0
    tm_method: str = "breslauer86"
    na: float = 50.0
    mg: float = 0.0
    dnac: float = 250.0

    # Oligo length
    target_length: int = 60

    # Repeat detection
    repeat_len: int = 8

    # Misprime detection
    misprime_len: int = 18
    misprime_tip: int = 6
    max_mismatches: int = 8

    # Scoring weights
    wt_tm: float = 1.0
    wt_misprime: float = 1.0
    wt_length: float = 1.0
    wt_repeat: float = 1.0
    wt_gc: float = 1.0
    wt_at: float = 1.0


# ============================================================================
# Sequence analysis (computed once)
# ============================================================================

def analyze_sequence(seq: str, config: ScoringConfig) -> SequenceAnalysis:
    """
    Pre-compute all placement-independent scores for a DNA sequence.

    This is O(n^2) for repeats and misprimes but only runs once.
    """
    rc = reverse_complement(seq)
    n = len(seq)
    analysis = SequenceAnalysis(seq=seq, rc=rc, length=n)

    analysis.potential_misprimes = _find_potential_misprimes(
        seq, rc, config.misprime_len, config.misprime_tip, config.max_mismatches
    )
    analysis.repeat_scores = _find_repeats(seq, rc, n, config.repeat_len)
    analysis.gc_scores = _find_gc_stretches(seq, n)
    analysis.at_scores = _find_at_stretches(seq, n)

    norm = 20.0 / n
    analysis.total_repeat = sum(analysis.repeat_scores) * norm
    analysis.total_gc = sum(analysis.gc_scores) * norm
    analysis.total_at = sum(analysis.at_scores) * norm

    return analysis


def _find_potential_misprimes(
    seq: str, rc: str, mp_len: int, mp_tip: int, max_mm: int
) -> list[PotentialMisprime]:
    """
    Scan for all position pairs that could cause mispriming.

    A potential misprime requires:
      1. The 3' tips (last mp_tip nt for sense, first mp_tip for antisense)
         are identical
      2. The full mp_len window has <= max_mm mismatches

    Four geometries are checked:
      Type 1 (direct-sense):      sense 3' tip at i matches sense at j
      Type 2 (inverse-sense):     antisense 3' tip at i matches sense RC at j
      Type 3 (inverse-antisense): sense 3' tip at i matches antisense RC at j
      Type 4 (direct-antisense):  antisense 3' tip at i matches antisense at j
    """
    n = len(seq)
    results = []
    tip_offset = mp_len - mp_tip  # offset from window start to sense 3' tip

    for i in range(n - mp_len + 1):
        for j in range(i, n - mp_len + 1):
            win_i = seq[i : i + mp_len]
            win_j = seq[j : j + mp_len]

            # Direct comparison (same strand)
            # Exclude i==j: a sequence trivially matches itself, not a real misprime.
            # (The Fortran's HMatchNum returns FALSE when pos1==pos2 for dir=1.)
            if i != j and _homologous(win_i, win_j, max_mm):
                # Type 1: sense 3' tips match
                if win_i[tip_offset:] == win_j[tip_offset:]:
                    results.append(PotentialMisprime(i, j, 1))
                # Type 4: antisense 3' tips match (first mp_tip nt)
                if win_i[:mp_tip] == win_j[:mp_tip]:
                    results.append(PotentialMisprime(i, j, 4))

            # Inverse comparison (complementary strands)
            # i==j IS allowed here: detects palindromic sequences that could
            # self-prime by forming hairpins.
            rc_win_j = reverse_complement(win_j)
            if _homologous(win_i, rc_win_j, max_mm):
                # Type 2: antisense 3' tip at i matches RC tip at j
                if win_i[:mp_tip] == rc_win_j[:mp_tip]:
                    results.append(PotentialMisprime(i, j, 2))
                # Type 3: sense 3' tip at i matches RC tip at j
                if win_i[tip_offset:] == rc_win_j[tip_offset:]:
                    results.append(PotentialMisprime(i, j, 3))

    return results


def _homologous(a: str, b: str, max_mm: int) -> bool:
    """Check if two equal-length sequences have <= max_mm mismatches."""
    mm = 0
    for ca, cb in zip(a, b):
        if ca != cb:
            mm += 1
            if mm > max_mm:
                return False
    return True


def _find_repeats(seq: str, rc: str, n: int, rep_len: int) -> list[float]:
    """
    Find all direct and inverted repeats of length >= rep_len.
    Repeats are expanded to their full extent once a seed match is found.
    Returns per-position repeat scores.
    """
    scores = [0.0] * n
    found = set()  # track (pos1, pos2, direction) to avoid duplicates

    for i in range(n - rep_len + 1):
        seed_i = seq[i : i + rep_len]

        for j in range(i, n - rep_len + 1):
            # Direct repeat: exclude i==j (a sequence is trivially identical to itself)
            if i != j and seq[j : j + rep_len] == seed_i:
                key = (i, j, 1)
                if key not in found:
                    found.add(key)
                    start1, start2, length = _expand_repeat(seq, i, j, rep_len, 1)
                    for k in range(length):
                        if start1 + k < n:
                            scores[start1 + k] += 1
                        if start2 + k < n:
                            scores[start2 + k] += 1

            # Inverted repeat: i==j IS included (palindromic sequences like
            # restriction sites can form hairpins during assembly PCR)
            rc_seed_j = reverse_complement(seq[j : j + rep_len])
            if rc_seed_j == seed_i:
                key = (i, j, -1)
                if key not in found:
                    found.add(key)
                    start1, start2, length = _expand_repeat(seq, i, j, rep_len, -1)
                    for k in range(length):
                        if start1 + k < n:
                            scores[start1 + k] += 1
                        if start2 + k < n:
                            scores[start2 + k] += 1

    return scores


def _expand_repeat(
    seq: str, pos1: int, pos2: int, rep_len: int, direction: int
) -> tuple[int, int, int]:
    """
    Expand a repeat match to its full extent in both directions.

    Parameters
    ----------
    seq : str
        Full DNA sequence.
    pos1, pos2 : int
        Starting positions of the seed match.
    rep_len : int
        Length of the initial seed.
    direction : int
        1 for direct repeat, -1 for inverted repeat.

    Returns
    -------
    tuple of (start1, start2, length)
        Expanded repeat coordinates.
    """
    n = len(seq)
    s1, s2 = pos1, pos2
    length = rep_len

    if direction == 1:
        # Expand left
        while s1 > 0 and s2 > 0:
            if seq[s1 - 1] == seq[s2 - 1]:
                s1 -= 1
                s2 -= 1
                length += 1
            else:
                break
        # Expand right
        e1, e2 = pos1 + rep_len, pos2 + rep_len
        while e1 < n and e2 < n:
            if seq[e1] == seq[e2]:
                e1 += 1
                e2 += 1
                length += 1
            else:
                break
    else:
        # Inverted repeat: expand outward
        # pos1 aligns with RC of pos2, so expanding left from pos1
        # means expanding right from pos2+rep_len and vice versa
        while s1 > 0 and (pos2 + rep_len + (pos1 - s1)) < n:
            left = s1 - 1
            right = pos2 + rep_len + (pos1 - s1)
            if seq[left] == _complement_base(seq[right]):
                s1 -= 1
                length += 1
            else:
                break
        while (pos1 + rep_len + (pos2 - s2)) < n and s2 > 0:
            right = pos1 + rep_len + (pos2 - s2)
            left = s2 - 1
            if seq[right] == _complement_base(seq[left]):
                s2 -= 1
                length += 1
            else:
                break

    return s1, s2, length


_COMP_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _complement_base(b: str) -> str:
    return _COMP_MAP[b]


def _find_gc_stretches(seq: str, n: int) -> list[float]:
    """Find 8-nt windows of pure GC content. Returns per-position scores."""
    scores = [0.0] * n
    gc = {"G", "C"}
    for i in range(n - 7):
        if all(b in gc for b in seq[i : i + 8]):
            for j in range(i, i + 8):
                scores[j] += 1
    return scores


def _find_at_stretches(seq: str, n: int) -> list[float]:
    """Find 8-nt windows of pure AT content. Returns per-position scores."""
    scores = [0.0] * n
    at = {"A", "T"}
    for i in range(n - 7):
        if all(b in at for b in seq[i : i + 8]):
            for j in range(i, i + 8):
                scores[j] += 1
    return scores


# ============================================================================
# Placement-dependent scoring (computed per trial)
# ============================================================================

def score_placement(
    analysis: SequenceAnalysis,
    overlaps: list[tuple[int, int]],
    config: ScoringConfig,
) -> Scores:
    """
    Score a specific overlap placement.

    Parameters
    ----------
    analysis : SequenceAnalysis
        Pre-computed sequence analysis.
    overlaps : list of (start, end) tuples
        0-indexed positions of each overlap region [start, end] inclusive.
    config : ScoringConfig
        Scoring parameters and weights.

    Returns
    -------
    Scores
        Per-component and total scores.
    """
    n = analysis.length
    norm = 20.0 / n

    scores = Scores()

    # Placement-dependent scores
    scores.tm = _score_tm(analysis.seq, overlaps, config) * norm
    scores.misprime = _score_misprimes(analysis, overlaps, config) * norm
    scores.length = _score_length(overlaps, config.target_length, n) * norm

    # Placement-independent scores (pre-computed)
    scores.repeat = analysis.total_repeat
    scores.gc = analysis.total_gc
    scores.at = analysis.total_at

    # Weighted total
    scores.total = (
        config.wt_tm * scores.tm
        + config.wt_misprime * scores.misprime
        + config.wt_length * scores.length
        + config.wt_repeat * scores.repeat
        + config.wt_gc * scores.gc
        + config.wt_at * scores.at
    )

    return scores


def _score_tm(seq: str, overlaps: list[tuple[int, int]], config: ScoringConfig) -> float:
    """
    Score Tm deviation across all overlaps.

    Within the tolerance band: score = 0
    Outside: score = ((deviation - tolerance)^2) / 10
    """
    total = 0.0
    for start, end in overlaps:
        overlap_seq = seq[start : end + 1]
        overlap_tm = tm(overlap_seq, method=config.tm_method,
                        na=config.na, mg=config.mg, dnac=config.dnac)
        diff = abs(config.target_tm - overlap_tm)
        if diff > config.tm_tolerance:
            penalty = max(1.0, diff - config.tm_tolerance)
            total += (penalty ** 2) / 10.0
    return total


def _score_misprimes(
    analysis: SequenceAnalysis,
    overlaps: list[tuple[int, int]],
    config: ScoringConfig,
) -> float:
    """
    Check which potential misprimes are 'actual' given the overlap placement.

    A potential misprime becomes actual when one of its positions aligns
    with an oligo 3' end. Oligo 3' ends exist only at even-indexed overlaps
    (0-indexed: 0, 2, 4, ...).

    At each such overlap:
      - Sense oligo 3' end is at the right edge (overlap end position)
      - Antisense oligo 3' end is at the left edge (overlap start position)
    """
    mp_len = config.misprime_len
    total = 0.0

    # Collect oligo 3' end positions as misprime window starts
    # sense_ends: start of mp_len window whose RIGHT end is the sense 3' end
    # anti_ends: start of mp_len window whose LEFT end is the antisense 3' end
    sense_ends = set()
    anti_ends = set()

    for idx, (start, end) in enumerate(overlaps):
        if idx % 2 == 0:  # only even-indexed overlaps have 3' ends
            sense_window_start = end - mp_len + 1
            if sense_window_start >= 0:
                sense_ends.add(sense_window_start)
            anti_ends.add(start)

    # Check each potential misprime
    for mp in analysis.potential_misprimes:
        scored = False
        if mp.mtype == 1:  # direct-sense: sense 3' primes at sense location
            if mp.pos1 in sense_ends or mp.pos2 in sense_ends:
                scored = True
        elif mp.mtype == 4:  # direct-antisense: antisense 3' primes at antisense
            if mp.pos1 in anti_ends or mp.pos2 in anti_ends:
                scored = True
        elif mp.mtype == 2:  # inverse-sense
            if mp.pos1 in anti_ends or mp.pos2 in sense_ends:
                scored = True
        elif mp.mtype == 3:  # inverse-antisense
            if mp.pos1 in sense_ends or mp.pos2 in anti_ends:
                scored = True

        if scored:
            # Each actual misprime contributes mp_len * 2 to per-position scores
            # (both windows get scored)
            total += 2 * mp_len

    return total


def _score_length(
    overlaps: list[tuple[int, int]],
    target_length: int,
    seq_len: int,
) -> float:
    """
    Score oligo length deviations.

    Each oligo spans from one overlap start to the next overlap end.
    Penalty = (overrun + 2)^2 applied to every position in the oligo.
    First and last oligos are exempt.
    """
    total = 0.0
    n_overlaps = len(overlaps)

    # Internal oligos span from OlapsPos[i-1][0] to OlapsPos[i][1]
    for i in range(1, n_overlaps):
        oligo_len = overlaps[i][1] - overlaps[i - 1][0] + 1
        overrun = oligo_len - target_length
        if overrun > 0:
            penalty = (overrun + 2) ** 2
            total += penalty * oligo_len

    return total
