#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=compleasm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=compleasm.%j.out
#SBATCH --error=compleasm.%j.err
#SBATCH --account=lp_svbelleghem

set -euo pipefail

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

ASSEMBLY_DIR="/scratch/leuven/357/vsc35707/sw-assembly/T2T_assembly"
COMPLEASM_DIR="/data/leuven/357/vsc35707/compleasm_kit"
CONDA_DIR="/data/leuven/357/vsc35707/miniconda3"

PURGED_ASSEMBLY="${ASSEMBLY_DIR}/purged.fa"
UNPURGED_ASSEMBLY="${ASSEMBLY_DIR}/Pogonus_SW.asm.hic.p_ctg.fa"

RESULTS_DIR="${ASSEMBLY_DIR}/compleasm_results"
LIBRARY_DIR="/data/leuven/357/vsc35707/compleasm_lineages"

LINEAGE="coleoptera"
THREADS="${SLURM_CPUS_PER_TASK:-24}"

# ------------------------------------------------------------
# Activate environment
# ------------------------------------------------------------

source "${CONDA_DIR}/etc/profile.d/conda.sh"
conda activate compleasm_env

# Prevent user-installed Python packages from overriding Conda packages
unset PYTHONPATH
unset PYTHONHOME
export PYTHONNOUSERSITE=1

# ------------------------------------------------------------
# Create output directories
# ------------------------------------------------------------

mkdir -p "$RESULTS_DIR"
mkdir -p "$LIBRARY_DIR"

# ------------------------------------------------------------
# Job information
# ------------------------------------------------------------

echo "============================================================"
echo "Compleasm assembly assessment"
echo "============================================================"
echo "Node:       $(hostname)"
echo "Start time: $(date)"
echo "Threads:    $THREADS"
echo "Lineage:    $LINEAGE"
echo "Python:     $(which python)"
echo

# ------------------------------------------------------------
# Verify Python environment
# ------------------------------------------------------------

python - <<'PY'
import sys
import numpy
import pandas

print("Python executable:", sys.executable)
print("Python version:", sys.version)
print("NumPy version:", numpy.__version__)
print("NumPy path:", numpy.__file__)
print("pandas version:", pandas.__version__)
print("pandas path:", pandas.__file__)

major, minor = map(int, numpy.__version__.split(".")[:2])

if (major, minor) < (1, 26):
    raise RuntimeError(
        f"NumPy >=1.26 is required, but Python loaded {numpy.__version__}"
    )
PY

echo

# ------------------------------------------------------------
# Check required files
# ------------------------------------------------------------

for file in \
    "${COMPLEASM_DIR}/compleasm.py" \
    "${COMPLEASM_DIR}/miniprot" \
    "${COMPLEASM_DIR}/hmmsearch" \
    "$PURGED_ASSEMBLY" \
    "$UNPURGED_ASSEMBLY"
do
    if [[ ! -s "$file" ]]; then
        echo "ERROR: Required file is missing or empty:" >&2
        echo "$file" >&2
        exit 1
    fi
done

chmod +x \
    "${COMPLEASM_DIR}/compleasm.py" \
    "${COMPLEASM_DIR}/miniprot" \
    "${COMPLEASM_DIR}/hmmsearch"

echo "Input assemblies:"
ls -lh "$PURGED_ASSEMBLY" "$UNPURGED_ASSEMBLY"
echo

# ------------------------------------------------------------
# Remove stale temporary download lock
# ------------------------------------------------------------

rm -f "${LIBRARY_DIR}/placement_files.tmp"

# ------------------------------------------------------------
# Download/check Coleoptera lineage
# ------------------------------------------------------------

echo "============================================================"
echo "Downloading or checking Compleasm lineage: $LINEAGE"
echo "============================================================"

python "${COMPLEASM_DIR}/compleasm.py" download \
    "$LINEAGE" \
    -L "$LIBRARY_DIR"

echo

# ------------------------------------------------------------
# Run Compleasm on unpurged assembly
# ------------------------------------------------------------

echo "============================================================"
echo "Running Compleasm on UNPURGED assembly"
echo "============================================================"
echo "Assembly: $UNPURGED_ASSEMBLY"
echo "Start:    $(date)"

rm -rf "${RESULTS_DIR}/unpurged"

python "${COMPLEASM_DIR}/compleasm.py" run \
    -a "$UNPURGED_ASSEMBLY" \
    -o "${RESULTS_DIR}/unpurged" \
    -l "$LINEAGE" \
    -L "$LIBRARY_DIR" \
    -t "$THREADS"

echo "Unpurged assembly finished: $(date)"
echo

# ------------------------------------------------------------
# Run Compleasm on purged assembly
# ------------------------------------------------------------

echo "============================================================"
echo "Running Compleasm on PURGED assembly"
echo "============================================================"
echo "Assembly: $PURGED_ASSEMBLY"
echo "Start:    $(date)"

rm -rf "${RESULTS_DIR}/purged"

python "${COMPLEASM_DIR}/compleasm.py" run \
    -a "$PURGED_ASSEMBLY" \
    -o "${RESULTS_DIR}/purged" \
    -l "$LINEAGE" \
    -L "$LIBRARY_DIR" \
    -t "$THREADS"

echo "Purged assembly finished: $(date)"
echo

# ------------------------------------------------------------
# Print summaries
# ------------------------------------------------------------

echo "============================================================"
echo "COMPLEASM SUMMARY: UNPURGED ASSEMBLY"
echo "============================================================"

if [[ -s "${RESULTS_DIR}/unpurged/summary.txt" ]]; then
    cat "${RESULTS_DIR}/unpurged/summary.txt"
else
    echo "WARNING: No summary found for unpurged assembly."
    find "${RESULTS_DIR}/unpurged" -maxdepth 2 -type f -print
fi

echo
echo "============================================================"
echo "COMPLEASM SUMMARY: PURGED ASSEMBLY"
echo "============================================================"

if [[ -s "${RESULTS_DIR}/purged/summary.txt" ]]; then
    cat "${RESULTS_DIR}/purged/summary.txt"
else
    echo "WARNING: No summary found for purged assembly."
    find "${RESULTS_DIR}/purged" -maxdepth 2 -type f -print
fi

echo
echo "============================================================"
echo "Compleasm analyses completed successfully"
echo "End time: $(date)"
echo "Results:"
echo "  ${RESULTS_DIR}/unpurged"
echo "  ${RESULTS_DIR}/purged"
echo "============================================================"
