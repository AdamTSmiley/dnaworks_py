#!/bin/bash
#
# run_benchmark.sh — Run Fortran DNAWorks and Python dnaworks_py on all test genes
#
# Usage:
#   ./run_benchmark.sh /path/to/dnaworks /path/to/dnaworks_py/cli.py
#
# Before running:
#   1. cd into the benchmark directory (where example_fastas/ lives)
#   2. Run: python scripts/generate_inputs.py
#   3. Then run this script

set -e

DNAWORKS_BIN="${1:-$HOME/DNAWorks/dnaworks}"
PYTHON_CLI="${2:-$HOME/dnaworks_py/cli.py}"
TIMING_RUNS="${3:-10}"

if [ ! -f "$DNAWORKS_BIN" ]; then
    echo "ERROR: Fortran binary not found at $DNAWORKS_BIN"
    echo "Usage: $0 /path/to/dnaworks /path/to/cli.py [timing_runs]"
    exit 1
fi

if [ ! -f "$PYTHON_CLI" ]; then
    echo "ERROR: Python CLI not found at $PYTHON_CLI"
    exit 1
fi

echo "================================================================"
echo "DNAWorks Benchmark: Fortran vs Python"
echo "================================================================"
echo "Fortran binary: $DNAWORKS_BIN"
echo "Python CLI:     $PYTHON_CLI"
echo "Timing runs:    $TIMING_RUNS per gene"
echo ""

# Create output directories
mkdir -p fortran_outputs python_outputs results

GENES=$(ls fortran_inputs/)

# ---- Run Fortran ----
echo "Running Fortran DNAWorks..."
echo "gene,run,time_seconds" > results/fortran_timing.csv

for gene in $GENES; do
    echo -n "  $gene: "
    cd fortran_inputs/$gene

    for run in $(seq 1 $TIMING_RUNS); do
        # Remove old output
        rm -f LOGFILE.txt

        # Time the run
        start=$(python3 -c "import time; print(time.time())")
        "$DNAWORKS_BIN" > /dev/null 2>&1
        end=$(python3 -c "import time; print(time.time())")
        elapsed=$(python3 -c "print(round($end - $start, 4))")
        echo "$gene,$run,$elapsed" >> ../../results/fortran_timing.csv
    done

    # Copy the last LOGFILE.txt to outputs
    cp LOGFILE.txt ../../fortran_outputs/${gene}.txt 2>/dev/null || echo "(no LOGFILE.txt)"
    cd ../..
    echo "done"
done

echo ""

# ---- Run Python ----
echo "Running Python dnaworks_py..."
echo "gene,run,time_seconds" > results/python_timing.csv

for gene in $GENES; do
    echo -n "  $gene: "
    fasta="example_fastas/${gene}.fasta"

    for run in $(seq 1 $TIMING_RUNS); do
        start=$(python3 -c "import time; print(time.time())")
        python3 "$PYTHON_CLI" "$fasta" --tm 62 --tm-tolerance 2 --length 60 \
            --json --quiet > /dev/null 2>&1
        end=$(python3 -c "import time; print(time.time())")
        elapsed=$(python3 -c "print(round($end - $start, 4))")
        echo "$gene,$run,$elapsed" >> results/python_timing.csv
    done

    # Save one JSON output for comparison
    python3 "$PYTHON_CLI" "$fasta" --tm 62 --tm-tolerance 2 --length 60 \
        --json --quiet > python_outputs/${gene}.json 2>/dev/null
    echo "done"
done

echo ""
echo "================================================================"
echo "Benchmark complete!"
echo "  Fortran outputs: fortran_outputs/"
echo "  Python outputs:  python_outputs/"
echo "  Timing data:     results/"
echo ""
echo "Next: python scripts/parse_fortran.py"
echo "Then: python scripts/compare.py"
echo "================================================================"
