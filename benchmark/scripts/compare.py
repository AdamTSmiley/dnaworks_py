#!/usr/bin/env python3
"""
Compare Fortran DNAWorks and Python dnaworks_py outputs.

Each tool is judged by its own Tm method:
  - Python overlaps measured with Breslauer86
  - Fortran overlaps measured with the exact Fortran TmCalc replication

Statistical tests: Mann-Whitney U (pooled overlaps) and Wilcoxon signed-rank
(paired per-gene means) for Tm deviation, overlap GC, overlap length, and ΔG.
"""

import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# ---- Tm calculators ----

def tm_breslauer(seq):
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.Seq import Seq
    if len(seq) <= 7: return 0.0
    return mt.Tm_NN(Seq(seq), nn_table=mt.DNA_NN1, saltcorr=5,
                    Na=50, Mg=0, dnac1=125, dnac2=125)

_NN_HT = {
    "AA":(-7.9,-22.2473),"AT":(-7.2,-20.38082),"AC":(-8.4,-22.44082),"AG":(-7.8,-21.02469),
    "TA":(-7.2,-21.34081),"TT":(-7.9,-22.2473),"TC":(-8.2,-22.24469),"TG":(-8.5,-22.73082),
    "CA":(-8.5,-22.73082),"CT":(-7.8,-21.02469),"CC":(-8.0,-19.8612),"CG":(-10.6,-27.17776),
    "GA":(-8.2,-22.24469),"GT":(-8.4,-22.44082),"GC":(-9.8,-24.37776),"GG":(-8.0,-19.8612),
}
_R = 1.9872
_OC = _R * math.log((2e-7/100)/2)

def tm_fortran_exact(seq):
    seq = seq.upper(); n = len(seq)
    if n <= 7: return 0.0
    dh, ds = 0.2, -5.68
    for i in range(n-1):
        pdh, pds = _NN_HT[seq[i:i+2]]; dh += pdh; ds += pds
    if seq[0] in "AT": dh += 2.2; ds += 6.935
    if seq[-1] in "AT": dh += 2.2; ds += 6.935
    comp = {"A":"T","T":"A","C":"G","G":"C"}
    if seq == "".join(comp[b] for b in reversed(seq)):
        ds += -1.4; ds += _R * math.log(2e-7/100)
    else:
        ds += _OC
    return (1000*dh/ds) - 273.15

def fold_dg(seq):
    import RNA; _, mfe = RNA.fold(seq); return mfe

# ---- Data loading ----

def load_results(fortran_dir, python_dir, fasta_dir):
    results = {}
    for pj in sorted(Path(python_dir).glob("*.json")):
        name = pj.stem
        fj = Path(fortran_dir) / f"{name}.json"
        fa = Path(fasta_dir) / f"{name}.fasta"
        if not fj.exists(): continue
        with open(pj) as f: py = json.load(f)
        with open(fj) as f: ft = json.load(f)
        gene_seq = ""
        with open(fa) as f:
            for line in f:
                if not line.startswith(">"): gene_seq += line.strip().upper()

        def proc(td):
            out = []
            for ov in td.get("overlaps", []):
                s = ov.get("sequence", gene_seq[ov["start"]-1:ov["end"]])
                if len(s) < 8 or len(s) > 200: continue
                gc = 100*sum(1 for b in s if b in "GC")/len(s)
                out.append({"seq":s,"length":len(s),"gc":gc,
                            "tm_bres":tm_breslauer(s),"tm_ft":tm_fortran_exact(s),
                            "dg":fold_dg(s)})
            return out

        results[name] = {
            "gene_length": len(gene_seq),
            "gc_pct": 100*sum(1 for b in gene_seq if b in "GC")/len(gene_seq),
            "python": {"num_oligos":py.get("num_oligos",0), "overlaps":proc(py),
                       "oligo_lengths":[o["length"] for o in py.get("oligos",[])]},
            "fortran": {"num_oligos":ft.get("num_oligos",0), "overlaps":proc(ft),
                        "oligo_lengths":[o["length"] for o in ft.get("oligos",[])]},
        }
    return results

def load_timing(d):
    t = {"fortran":{},"python":{}}
    for tool in t:
        p = Path(d)/f"{tool}_timing.csv"
        if not p.exists(): continue
        with open(p) as f:
            next(f)
            for line in f:
                pts = line.strip().split(",")
                if len(pts)==3: t[tool].setdefault(pts[0],[]).append(float(pts[2]))
    return t

# ---- Helpers ----
BLUE, ORANGE = "#2196F3", "#FF5722"
def gl(g,r): return f"{g}\n({r[g]['gene_length']}bp)"
def sg(r): return sorted(r.keys(), key=lambda g: r[g]["gene_length"])

def pval_str(p):
    if p < 0.001: return "p<0.001"
    elif p < 0.01: return f"p={p:.3f}"
    elif p < 0.05: return f"p={p:.2f}"
    else: return f"p={p:.2f} (ns)"

def dot_legend(ax):
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0],[0],marker="o",color="w",markerfacecolor=BLUE,markersize=8,label="Python (Breslauer)"),
        Line2D([0],[0],marker="o",color="w",markerfacecolor=ORANGE,markersize=8,label="Fortran (native Tm)"),
    ], fontsize=8)

# ---- Plots ----

def plot_tm_strip(results, outdir):
    """Single strip plot: Python overlaps by Breslauer, Fortran overlaps by Fortran-exact."""
    genes = sg(results); n = len(genes)
    fig, ax = plt.subplots(figsize=(max(10,n*1.2), 5))

    all_py, all_ft = [], []
    for i, gene in enumerate(genes):
        py_v = [o["tm_bres"] for o in results[gene]["python"]["overlaps"]]
        ft_v = [o["tm_ft"] for o in results[gene]["fortran"]["overlaps"]]
        all_py.extend(py_v); all_ft.extend(ft_v)
        if py_v:
            jx = np.full(len(py_v),i-0.15)+np.random.normal(0,0.03,len(py_v))
            ax.scatter(jx, py_v, color=BLUE, alpha=0.5, s=15, zorder=3)
        if ft_v:
            jx = np.full(len(ft_v),i+0.15)+np.random.normal(0,0.03,len(ft_v))
            ax.scatter(jx, ft_v, color=ORANGE, alpha=0.5, s=15, zorder=3)

    ax.axhline(y=62, color="gray", linestyle="--", linewidth=1)
    ax.axhspan(60, 64, alpha=0.08, color="gray")
    ax.set_xticks(range(n))
    ax.set_xticklabels([gl(g,results) for g in genes], fontsize=8)
    ax.set_ylabel("Overlap Tm (°C)")
    ax.set_title("Overlap Melting Temperatures (each tool measured by its own Tm method)")
    dot_legend(ax)

    # Stats annotation
    u, p = stats.mannwhitneyu(np.abs(np.array(all_py)-62), np.abs(np.array(all_ft)-62), alternative="two-sided")
    ax.text(0.02, 0.02, f"|Tm-62| Mann-Whitney: {pval_str(p)}", transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    plt.tight_layout()
    plt.savefig(outdir/"tm_strip.png", dpi=150); plt.close()
    print("  Saved tm_strip.png")


def plot_tm_deviation(results, outdir):
    """Bar plot of |Tm - 62| with stats."""
    genes = sg(results)
    fig, ax = plt.subplots(figsize=(max(10,len(genes)*1.2), 5))
    x = np.arange(len(genes)); w = 0.35

    py_m, ft_m, py_x, ft_x = [], [], [], []
    py_means_per_gene, ft_means_per_gene = [], []
    for gene in genes:
        py_d = [abs(o["tm_bres"]-62) for o in results[gene]["python"]["overlaps"]]
        ft_d = [abs(o["tm_ft"]-62) for o in results[gene]["fortran"]["overlaps"]]
        py_m.append(np.mean(py_d) if py_d else 0); ft_m.append(np.mean(ft_d) if ft_d else 0)
        py_x.append(np.max(py_d) if py_d else 0); ft_x.append(np.max(ft_d) if ft_d else 0)
        py_means_per_gene.append(np.mean(py_d) if py_d else 0)
        ft_means_per_gene.append(np.mean(ft_d) if ft_d else 0)

    ax.bar(x-w/2, py_m, w, label="Python (Breslauer)", color=BLUE, alpha=0.7)
    ax.bar(x+w/2, ft_m, w, label="Fortran (native Tm)", color=ORANGE, alpha=0.7)
    ax.scatter(x-w/2, py_x, color="#1565C0", marker="_", s=100, linewidths=2, zorder=4)
    ax.scatter(x+w/2, ft_x, color="#BF360C", marker="_", s=100, linewidths=2, zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels([gl(g,results) for g in genes], fontsize=8)
    ax.set_ylabel("|Tm - 62°C| (°C)")
    ax.set_title("Tm Deviation from Target (bars=mean, ticks=max)")
    ax.legend()
    ax.axhline(y=2.0, color="gray", linestyle=":", alpha=0.5)

    # Wilcoxon signed-rank on paired per-gene means
    if len(py_means_per_gene) >= 5:
        _, p = stats.wilcoxon(py_means_per_gene, ft_means_per_gene)
        ax.text(0.02, 0.95, f"Wilcoxon (paired means): {pval_str(p)}", transform=ax.transAxes,
                fontsize=8, va="top", bbox=dict(boxstyle="round",facecolor="wheat",alpha=0.5))

    plt.tight_layout()
    plt.savefig(outdir/"tm_deviation.png", dpi=150); plt.close()
    print("  Saved tm_deviation.png")


def plot_overlap_gc(results, outdir):
    genes = sg(results)
    fig, ax = plt.subplots(figsize=(max(10,len(genes)*1.2), 5))
    all_py, all_ft = [], []
    for i, gene in enumerate(genes):
        for vals, color, dx, collect in [
            ([o["gc"] for o in results[gene]["python"]["overlaps"]], BLUE, -0.15, all_py),
            ([o["gc"] for o in results[gene]["fortran"]["overlaps"]], ORANGE, 0.15, all_ft),
        ]:
            collect.extend(vals)
            if vals:
                jx = np.full(len(vals),i+dx)+np.random.normal(0,0.03,len(vals))
                ax.scatter(jx, vals, color=color, alpha=0.5, s=15, zorder=3)
    ax.axhspan(40, 60, alpha=0.05, color="green")
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels([gl(g,results) for g in genes], fontsize=8)
    ax.set_ylabel("Overlap GC Content (%)"); ax.set_title("GC Content of Overlap Regions")
    dot_legend(ax)

    u, p = stats.mannwhitneyu(all_py, all_ft, alternative="two-sided")
    ax.text(0.02, 0.02, f"Mann-Whitney: {pval_str(p)}", transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round",facecolor="wheat",alpha=0.5))
    plt.tight_layout()
    plt.savefig(outdir/"overlap_gc.png", dpi=150); plt.close()
    print("  Saved overlap_gc.png")


def plot_overlap_lengths(results, outdir):
    genes = sg(results)
    fig, ax = plt.subplots(figsize=(max(10,len(genes)*1.2), 5))
    all_py, all_ft = [], []
    for i, gene in enumerate(genes):
        for vals, color, dx, collect in [
            ([o["length"] for o in results[gene]["python"]["overlaps"]], BLUE, -0.15, all_py),
            ([o["length"] for o in results[gene]["fortran"]["overlaps"]], ORANGE, 0.15, all_ft),
        ]:
            collect.extend(vals)
            if vals:
                jx = np.full(len(vals),i+dx)+np.random.normal(0,0.03,len(vals))
                ax.scatter(jx, vals, color=color, alpha=0.5, s=15, zorder=3)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels([gl(g,results) for g in genes], fontsize=8)
    ax.set_ylabel("Overlap Length (nt)"); ax.set_title("Overlap Region Lengths")
    dot_legend(ax)

    u, p = stats.mannwhitneyu(all_py, all_ft, alternative="two-sided")
    ax.text(0.02, 0.02, f"Mann-Whitney: {pval_str(p)}", transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round",facecolor="wheat",alpha=0.5))
    plt.tight_layout()
    plt.savefig(outdir/"overlap_lengths.png", dpi=150); plt.close()
    print("  Saved overlap_lengths.png")


def plot_overlap_stability(results, outdir):
    genes = sg(results)
    fig, ax = plt.subplots(figsize=(max(10,len(genes)*1.2), 5))
    all_py, all_ft = [], []
    for i, gene in enumerate(genes):
        for vals, color, dx, collect in [
            ([o["dg"] for o in results[gene]["python"]["overlaps"]], BLUE, -0.15, all_py),
            ([o["dg"] for o in results[gene]["fortran"]["overlaps"]], ORANGE, 0.15, all_ft),
        ]:
            collect.extend(vals)
            if vals:
                jx = np.full(len(vals),i+dx)+np.random.normal(0,0.03,len(vals))
                ax.scatter(jx, vals, color=color, alpha=0.5, s=15, zorder=3)
    ax.axhline(y=0, color="gray", linestyle="--", linewidth=0.5)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels([gl(g,results) for g in genes], fontsize=8)
    ax.set_ylabel("ΔG (kcal/mol)")
    ax.set_title("Overlap Single-Strand Stability (ViennaRNA MFE; closer to 0 = better)")
    dot_legend(ax)

    u, p = stats.mannwhitneyu(all_py, all_ft, alternative="two-sided")
    ax.text(0.02, 0.02, f"Mann-Whitney: {pval_str(p)}", transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round",facecolor="wheat",alpha=0.5))
    plt.tight_layout()
    plt.savefig(outdir/"overlap_stability.png", dpi=150); plt.close()
    print("  Saved overlap_stability.png")


def plot_oligo_counts(results, outdir):
    """Oligo counts with count annotations on bars."""
    genes = sg(results)
    fig, ax = plt.subplots(figsize=(max(8,len(genes)*1.0), 4))
    x = np.arange(len(genes)); w = 0.35
    py_c = [results[g]["python"]["num_oligos"] for g in genes]
    ft_c = [results[g]["fortran"]["num_oligos"] for g in genes]
    bars1 = ax.bar(x-w/2, py_c, w, label="Python", color=BLUE, alpha=0.7)
    bars2 = ax.bar(x+w/2, ft_c, w, label="Fortran", color=ORANGE, alpha=0.7)
    # Annotate counts
    for bar, val in zip(bars1, py_c):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5, str(val),
                ha="center", va="bottom", fontsize=7, color="#1565C0")
    for bar, val in zip(bars2, ft_c):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5, str(val),
                ha="center", va="bottom", fontsize=7, color="#BF360C")
    ax.set_xticks(x)
    ax.set_xticklabels([gl(g,results) for g in genes], fontsize=8)
    ax.set_ylabel("Number of Oligos"); ax.set_title("Oligo Count Comparison"); ax.legend()
    plt.tight_layout()
    plt.savefig(outdir/"oligo_counts.png", dpi=150); plt.close()
    print("  Saved oligo_counts.png")


def plot_oligo_lengths(results, outdir):
    """Strip plot of internal oligo lengths (better than boxplot when IQR≈0)."""
    genes = sg(results)
    fig, ax = plt.subplots(figsize=(max(10,len(genes)*1.2), 5))

    all_py, all_ft = [], []
    for i, gene in enumerate(genes):
        py_l = results[gene]["python"]["oligo_lengths"]
        ft_l = results[gene]["fortran"]["oligo_lengths"]
        py_int = [l for l in (py_l[1:-1] if len(py_l)>2 else py_l) if l < 180]
        ft_int = [l for l in (ft_l[1:-1] if len(ft_l)>2 else ft_l) if l < 180]
        all_py.extend(py_int); all_ft.extend(ft_int)

        if py_int:
            jx = np.full(len(py_int),i-0.15)+np.random.normal(0,0.04,len(py_int))
            ax.scatter(jx, py_int, color=BLUE, alpha=0.4, s=12, zorder=3)
        if ft_int:
            jx = np.full(len(ft_int),i+0.15)+np.random.normal(0,0.04,len(ft_int))
            ax.scatter(jx, ft_int, color=ORANGE, alpha=0.4, s=12, zorder=3)

    ax.axhline(y=60, color="gray", linestyle="--", linewidth=1)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels([gl(g,results) for g in genes], fontsize=8)
    ax.set_ylabel("Internal Oligo Length (nt)"); ax.set_title("Internal Oligo Lengths (target=60 nt)")

    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0],[0],marker="o",color="w",markerfacecolor=BLUE,markersize=8,label="Python"),
        Line2D([0],[0],marker="o",color="w",markerfacecolor=ORANGE,markersize=8,label="Fortran"),
    ], fontsize=8)

    # Stats
    u, p = stats.mannwhitneyu(all_py, all_ft, alternative="two-sided")
    ax.text(0.02, 0.02, f"Mann-Whitney: {pval_str(p)}", transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round",facecolor="wheat",alpha=0.5))
    plt.tight_layout()
    plt.savefig(outdir/"oligo_lengths.png", dpi=150); plt.close()
    print("  Saved oligo_lengths.png")


def plot_runtime(timing, results, outdir):
    genes = sg(results)
    py_t, ft_t, labels = [], [], []
    for g in genes:
        if g in timing["python"] and g in timing["fortran"]:
            py_t.append(np.median(timing["python"][g]))
            ft_t.append(np.median(timing["fortran"][g]))
            labels.append(g)
    if not labels: print("  No timing data"); return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    x = np.arange(len(labels)); w = 0.35
    ax1.bar(x-w/2, py_t, w, label="Python", color=BLUE, alpha=0.7)
    ax1.bar(x+w/2, ft_t, w, label="Fortran", color=ORANGE, alpha=0.7)
    ax1.set_xticks(x)
    ax1.set_xticklabels([gl(g,results) for g in labels], fontsize=7)
    ax1.set_ylabel("Runtime (seconds)"); ax1.set_title("Runtime per Gene")
    ax1.legend(); ax1.set_yscale("log")

    lengths = [results[g]["gene_length"] for g in labels]
    ax2.scatter(lengths, py_t, color=BLUE, s=60, label="Python", zorder=3)
    ax2.scatter(lengths, ft_t, color=ORANGE, s=60, label="Fortran", zorder=3)
    ax2.set_xlabel("Sequence Length (bp)"); ax2.set_ylabel("Runtime (seconds)")
    ax2.set_title("Runtime vs Sequence Length"); ax2.legend(); ax2.set_yscale("log")
    plt.tight_layout()
    plt.savefig(outdir/"runtime.png", dpi=150); plt.close()
    print("  Saved runtime.png")


def print_summary_table(results, timing):
    genes = sg(results)
    print("\n" + "="*130)
    print(f"{'Gene':<12s} {'Len':>5s} {'GC%':>5s} | "
          f"{'Py #Ol':>6s} {'Ft #Ol':>6s} | "
          f"{'Py TmDv':>7s} {'Ft TmDv':>7s} | "
          f"{'Py TmSD':>7s} {'Ft TmSD':>7s} | "
          f"{'Py ΔG':>7s} {'Ft ΔG':>7s} | "
          f"{'Py OvGC':>7s} {'Ft OvGC':>7s} | "
          f"{'Py t(s)':>7s} {'Ft t(s)':>7s}")
    print("-"*130)
    for gene in genes:
        py = results[gene]["python"]; ft = results[gene]["fortran"]
        py_tm = [o["tm_bres"] for o in py["overlaps"]]
        ft_tm = [o["tm_ft"] for o in ft["overlaps"]]
        py_dg = [o["dg"] for o in py["overlaps"]]
        ft_dg = [o["dg"] for o in ft["overlaps"]]
        py_gc = [o["gc"] for o in py["overlaps"]]
        ft_gc = [o["gc"] for o in ft["overlaps"]]
        print(f"{gene:<12s} {results[gene]['gene_length']:>5d} {results[gene]['gc_pct']:>5.1f} | "
              f"{py['num_oligos']:>6d} {ft['num_oligos']:>6d} | "
              f"{np.mean(np.abs(np.array(py_tm)-62)):>7.2f} {np.mean(np.abs(np.array(ft_tm)-62)):>7.2f} | "
              f"{np.std(py_tm):>7.2f} {np.std(ft_tm):>7.2f} | "
              f"{np.mean(py_dg):>7.2f} {np.mean(ft_dg):>7.2f} | "
              f"{np.mean(py_gc):>7.1f} {np.mean(ft_gc):>7.1f} | "
              f"{np.median(timing['python'].get(gene,[0])):>7.3f} "
              f"{np.median(timing['fortran'].get(gene,[0])):>7.3f}")

    # Pooled statistics
    all_py_dev = [abs(o["tm_bres"]-62) for g in genes for o in results[g]["python"]["overlaps"]]
    all_ft_dev = [abs(o["tm_ft"]-62) for g in genes for o in results[g]["fortran"]["overlaps"]]
    all_py_dg = [o["dg"] for g in genes for o in results[g]["python"]["overlaps"]]
    all_ft_dg = [o["dg"] for g in genes for o in results[g]["fortran"]["overlaps"]]
    all_py_gc = [o["gc"] for g in genes for o in results[g]["python"]["overlaps"]]
    all_ft_gc = [o["gc"] for g in genes for o in results[g]["fortran"]["overlaps"]]
    all_py_len = [o["length"] for g in genes for o in results[g]["python"]["overlaps"]]
    all_ft_len = [o["length"] for g in genes for o in results[g]["fortran"]["overlaps"]]

    print("="*130)
    print(f"\nPooled statistics (Mann-Whitney U, two-sided):")
    print(f"  |Tm-62°C|  Py mean={np.mean(all_py_dev):.2f}  Ft mean={np.mean(all_ft_dev):.2f}  {pval_str(stats.mannwhitneyu(all_py_dev,all_ft_dev).pvalue)}")
    print(f"  ΔG         Py mean={np.mean(all_py_dg):.2f}  Ft mean={np.mean(all_ft_dg):.2f}  {pval_str(stats.mannwhitneyu(all_py_dg,all_ft_dg).pvalue)}")
    print(f"  GC%        Py mean={np.mean(all_py_gc):.1f}  Ft mean={np.mean(all_ft_gc):.1f}  {pval_str(stats.mannwhitneyu(all_py_gc,all_ft_gc).pvalue)}")
    print(f"  Ov length  Py mean={np.mean(all_py_len):.1f}  Ft mean={np.mean(all_ft_len):.1f}  {pval_str(stats.mannwhitneyu(all_py_len,all_ft_len).pvalue)}")

    # Paired per-gene Wilcoxon
    py_gene_dev = [np.mean([abs(o["tm_bres"]-62) for o in results[g]["python"]["overlaps"]]) for g in genes]
    ft_gene_dev = [np.mean([abs(o["tm_ft"]-62) for o in results[g]["fortran"]["overlaps"]]) for g in genes]
    py_gene_dg = [np.mean([o["dg"] for o in results[g]["python"]["overlaps"]]) for g in genes]
    ft_gene_dg = [np.mean([o["dg"] for o in results[g]["fortran"]["overlaps"]]) for g in genes]

    print(f"\nPaired statistics (Wilcoxon signed-rank on per-gene means, n={len(genes)}):")
    if len(genes) >= 5:
        print(f"  |Tm-62°C|  {pval_str(stats.wilcoxon(py_gene_dev,ft_gene_dev).pvalue)}")
        print(f"  ΔG         {pval_str(stats.wilcoxon(py_gene_dg,ft_gene_dg).pvalue)}")


def main():
    fortran_dir = Path(sys.argv[1]) if len(sys.argv)>1 else Path("fortran_outputs")
    python_dir = Path(sys.argv[2]) if len(sys.argv)>2 else Path("python_outputs")
    fasta_dir = Path(sys.argv[3]) if len(sys.argv)>3 else Path("example_fastas")
    timing_dir = Path(sys.argv[4]) if len(sys.argv)>4 else Path("results")
    outdir = Path(sys.argv[5]) if len(sys.argv)>5 else Path("results")
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading results and computing overlap properties...")
    results = load_results(fortran_dir, python_dir, fasta_dir)
    if not results: print("No results."); sys.exit(1)
    print(f"Loaded {len(results)} genes\n")
    timing = load_timing(timing_dir)

    print("Generating plots...")
    np.random.seed(42)
    plot_tm_strip(results, outdir)
    plot_tm_deviation(results, outdir)
    plot_overlap_gc(results, outdir)
    plot_overlap_lengths(results, outdir)
    plot_overlap_stability(results, outdir)
    plot_oligo_counts(results, outdir)
    plot_oligo_lengths(results, outdir)
    plot_runtime(timing, results, outdir)
    print_summary_table(results, timing)

if __name__ == "__main__":
    main()
