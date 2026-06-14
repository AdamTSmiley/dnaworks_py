"""
Simulated annealing optimizer for overlap placement.

When oligo lengths are randomized, tries many different overlap arrangements
and uses the Metropolis criterion to escape local minima. This finds better
placements than the deterministic fixed-length mode for tricky sequences
(high GC, many repeats, misprime-prone).

SA parameters match the Fortran DNAWorks defaults:
  - Initial temperature: 0.5
  - Cooling factor: 0.96 per drop
  - Up to 1000 temperature drops × 500 rounds per drop
  - Early exit on convergence
"""

import math
import random as _random
import sys
import time
from dataclasses import dataclass

from tm import tm as calc_tm
from scoring import (
    SequenceAnalysis,
    ScoringConfig,
    Scores,
    analyze_sequence,
    score_placement,
)
from overlaps import _place_overlaps, _extract_oligos, DesignResult


@dataclass
class SAConfig:
    """Simulated annealing parameters."""
    initial_temp: float = 0.5
    cooling_factor: float = 0.96
    max_temp_drops: int = 1000
    rounds_per_drop: int = 500
    max_success_per_drop: int = 50
    max_stale_rounds: int = 1500
    min_temp: float = 0.0001
    min_score: float = 0.001


def design_oligos_sa(
    seq: str,
    config: ScoringConfig,
    oligo_length: int = 60,
    nogaps: bool = False,
    sa_config: SAConfig | None = None,
    max_shifts: int = 100,
    seed: int | None = None,
    verbose: bool = False,
) -> DesignResult:
    """
    Design oligos using simulated annealing with random oligo lengths.

    Each SA round generates a new overlap arrangement with randomized
    oligo lengths (between 20 and oligo_length), evaluates it, and
    accepts or rejects by the Metropolis criterion.

    Parameters
    ----------
    seq : str
        Validated DNA sequence.
    config : ScoringConfig
        Scoring parameters.
    oligo_length : int
        Maximum oligo length. Actual lengths randomized between 20 and this.
    nogaps : bool
        No gaps between overlaps.
    sa_config : SAConfig or None
        SA parameters. Uses defaults if None.
    max_shifts : int
        Shifts per Generate_Overlaps call. Lower than deterministic mode
        (100 vs 1000) since we're running many SA rounds. Default 100.
    seed : int or None
        Random seed for reproducibility.
    verbose : bool
        Print progress to stderr.

    Returns
    -------
    DesignResult
    """
    if seed is not None:
        _random.seed(seed)
    if sa_config is None:
        sa_config = SAConfig()

    n = len(seq)
    tm_kw = dict(method=config.tm_method, na=config.na, mg=config.mg, dnac=config.dnac)

    # Pre-compute sequence analysis (done once)
    if verbose:
        print("Pre-computing sequence analysis...", file=sys.stderr)
    analysis = analyze_sequence(seq, config)

    # Get initial solution using deterministic mode
    if verbose:
        print("Finding initial solution...", file=sys.stderr)
    current_overlaps = _find_initial_placement(
        seq, n, oligo_length, nogaps, config, tm_kw, max_shifts, analysis
    )
    if current_overlaps is None:
        raise RuntimeError("Failed to find initial overlap arrangement.")

    current_scores = score_placement(analysis, current_overlaps, config)
    best_overlaps = current_overlaps
    best_scores = current_scores
    best_score_val = best_scores.total

    # SA loop
    temp = sa_config.initial_temp
    stale_count = 0
    total_rounds = 0

    start_time = time.time()

    for temp_drop in range(sa_config.max_temp_drops):
        n_success = 0

        for round_num in range(sa_config.rounds_per_drop):
            total_rounds += 1

            # Generate new placement with random lengths
            new_overlaps = _random_placement(
                seq, n, oligo_length, nogaps, config, tm_kw, max_shifts, analysis
            )

            if new_overlaps is None or len(new_overlaps) % 2 == 0:
                continue

            new_scores = score_placement(analysis, new_overlaps, config)
            gain = new_scores.total - current_scores.total

            # Metropolis acceptance
            if gain < 0 or _random.random() < math.exp(-gain / temp):
                current_overlaps = new_overlaps
                current_scores = new_scores
                n_success += 1

            # Track best ever
            if current_scores.total < best_score_val:
                best_overlaps = current_overlaps
                best_scores = current_scores
                best_score_val = best_scores.total
                stale_count = 0
            else:
                stale_count += 1

            # Progress reporting
            if verbose and total_rounds % 100 == 0:
                elapsed = time.time() - start_time
                print(
                    f"  {total_rounds:5d} rounds, best={best_score_val:.3f}, "
                    f"temp={temp:.4f}, elapsed={elapsed:.1f}s",
                    file=sys.stderr,
                )

            # Early exit: too many rounds without improvement
            if stale_count > sa_config.max_stale_rounds:
                if verbose:
                    print("Converged (no improvement).", file=sys.stderr)
                break

            # Early exit: enough successes this temperature
            if n_success >= sa_config.max_success_per_drop:
                break

        # Cool down
        temp *= sa_config.cooling_factor

        # Check exit conditions
        if stale_count > sa_config.max_stale_rounds:
            break
        if n_success == 0:
            if verbose:
                print("Converged (no acceptances at current temp).", file=sys.stderr)
            break
        if temp < sa_config.min_temp:
            if verbose:
                print("Converged (temperature floor).", file=sys.stderr)
            break
        if best_score_val < sa_config.min_score:
            if verbose:
                print("Converged (near-zero score).", file=sys.stderr)
            break

    if verbose:
        elapsed = time.time() - start_time
        print(
            f"SA complete: {total_rounds} rounds in {elapsed:.1f}s, "
            f"best score={best_score_val:.3f}",
            file=sys.stderr,
        )

    # Build result
    oligos = _extract_oligos(seq, best_overlaps)
    overlap_tms = [
        calc_tm(seq[s : e + 1], **tm_kw) for s, e in best_overlaps
    ]

    return DesignResult(
        overlaps=best_overlaps,
        oligos=oligos,
        scores=best_scores,
        overlap_tms=overlap_tms,
        num_overlaps=len(best_overlaps),
        num_oligos=len(oligos),
    )


def _find_initial_placement(seq, n, oligo_len, nogaps, config, tm_kw, max_shifts, analysis):
    """Find the best initial placement using deterministic mode."""
    best = None
    best_score = float("inf")

    for shift in range(max_shifts):
        overlaps = _place_overlaps(
            seq, n, shift, oligo_len, nogaps, False,
            config.target_tm, config.tm_tolerance, tm_kw,
        )
        if overlaps and len(overlaps) % 2 == 1:
            scores = score_placement(analysis, overlaps, config)
            if scores.total < best_score:
                best_score = scores.total
                best = overlaps

    return best


def _random_placement(seq, n, oligo_len, nogaps, config, tm_kw, max_shifts, analysis):
    """Generate a single placement with random oligo lengths from a random shift."""
    # Pick a random starting shift and generate one placement with random lengths.
    # The SA loop provides the optimization — we just need diverse candidates.
    for _ in range(max_shifts):
        shift = _random.randint(0, max(1, n // 10))
        overlaps = _place_overlaps(
            seq, n, shift, oligo_len, nogaps, True,
            config.target_tm, config.tm_tolerance, tm_kw,
        )
        if overlaps and len(overlaps) % 2 == 1:
            return overlaps
    return None
