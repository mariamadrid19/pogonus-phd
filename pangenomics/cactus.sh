#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=pangenome
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH -o pangenome.%j.out
#SBATCH -A lp_svbelleghem

###############################################################################
# Paths
###############################################################################

WORKDIR="/lustre1/scratch/357/vsc35707/pangenome"

LW_FASTA="${WORKDIR}/Pchalceus_LW_final.fasta"
SW_FASTA="${WORKDIR}/Pchalceus_SW_final.fasta"
LW_FAI="${LW_FASTA}.fai"
CHRFILE="${WORKDIR}/LW.chromosomes.txt"

# Official Cactus 3.2.1 installation.
CACTUS_ROOT="/vsc-hard-mounts/leuven-data/357/vsc35707/cactus-bin-v3.2.1"
CACTUS_ENV="${CACTUS_ROOT}/venv-cactus-v3.2.1"
CACTUS="${CACTUS_ENV}/bin/cactus-pangenome"

# Use fresh names; do not reuse the failed development-version workflow.
SEQFILE="${WORKDIR}/Pchalceus_LWref-v2.seqfile"
JOBSTORE="${WORKDIR}/Pchalceus-LWref-jobstore-v2"
OUTDIR="${WORKDIR}/Pchalceus-LWref-pangenome-v2"
TOIL_WORKDIR="${WORKDIR}/toil-LWref-work-v2"
LOGFILE="${WORKDIR}/Pchalceus-LWref-cactus-v2.log"

###############################################################################
# Environment setup
###############################################################################

cd "${WORKDIR}"

# Activate the matched Cactus 3.2.1 environment. This also configures PATH,
# PYTHONPATH and LD_LIBRARY_PATH for the precompiled Cactus distribution.
[[ -s "${CACTUS_ENV}/bin/activate" ]] || {
    echo "ERROR: Cactus environment activation file not found:" >&2
    echo "${CACTUS_ENV}/bin/activate" >&2
    exit 1
}

source "${CACTUS_ENV}/bin/activate"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Working directory: ${WORKDIR}"
echo "Cactus executable: $(command -v cactus-pangenome)"
echo "Python executable: $(command -v python3)"

python3 -c '
import toil
print("Toil version:", toil.version.version)
'

###############################################################################
# Validate programs
###############################################################################

[[ -x "${CACTUS}" ]] || {
    echo "ERROR: cactus-pangenome is not executable:" >&2
    echo "${CACTUS}" >&2
    exit 1
}

command -v singularity >/dev/null 2>&1 || {
    echo "ERROR: singularity is unavailable in this batch environment." >&2
    echo "Load the appropriate Apptainer/Singularity module in this script." >&2
    exit 1
}

echo "Singularity executable: $(command -v singularity)"

###############################################################################
# Validate input files
###############################################################################

[[ -s "${LW_FASTA}" ]] || {
    echo "ERROR: LW FASTA not found or empty: ${LW_FASTA}" >&2
    exit 1
}

[[ -s "${SW_FASTA}" ]] || {
    echo "ERROR: SW FASTA not found or empty: ${SW_FASTA}" >&2
    exit 1
}

[[ -s "${LW_FAI}" ]] || {
    echo "ERROR: LW FASTA index not found or empty: ${LW_FAI}" >&2
    exit 1
}

[[ -s "${CHRFILE}" ]] || {
    echo "ERROR: LW chromosome list not found or empty: ${CHRFILE}" >&2
    exit 1
}

###############################################################################
# Validate the 11 LW reference chromosomes
###############################################################################

CHR_COUNT=$(awk 'NF {count++} END {print count+0}' "${CHRFILE}")

if [[ "${CHR_COUNT}" -ne 11 ]]; then
    echo "ERROR: ${CHRFILE} contains ${CHR_COUNT} nonempty entries; expected 11." >&2
    exit 1
fi

DUPLICATES=$(
    awk 'NF {print $1}' "${CHRFILE}" |
        sort |
        uniq -d
)

if [[ -n "${DUPLICATES}" ]]; then
    echo "ERROR: Duplicate chromosome identifiers found:" >&2
    echo "${DUPLICATES}" >&2
    exit 1
fi

# Check every chromosome identifier against the first column of the FASTA index.
if ! awk '
    NR == FNR {
        fasta_ids[$1] = 1
        next
    }

    NF && !($1 in fasta_ids) {
        print "ERROR: Chromosome ID not found in LW FASTA: " $1 > "/dev/stderr"
        missing = 1
    }

    END {
        exit missing
    }
' "${LW_FAI}" "${CHRFILE}"
then
    exit 1
fi

mapfile -t REF_CONTIGS < <(
    awk 'NF {print $1}' "${CHRFILE}"
)

###############################################################################
# Protect against mixing this run with previous attempts
###############################################################################

if [[ -e "${JOBSTORE}" ]]; then
    echo "ERROR: The new job store already exists:" >&2
    echo "${JOBSTORE}" >&2
    echo "Do not mix different Cactus runs in one job store." >&2
    exit 1
fi

if [[ -e "${OUTDIR}" ]]; then
    echo "ERROR: The new output directory already exists:" >&2
    echo "${OUTDIR}" >&2
    echo "Do not mix results from different Cactus runs." >&2
    exit 1
fi

# Toil requires the work directory to exist before it starts.
mkdir -p "${TOIL_WORKDIR}"

###############################################################################
# Create the Minigraph-Cactus sequence file
###############################################################################

# LW is first and is explicitly selected as the reference.
printf "LW\t%s\nSW\t%s\n" \
    "${LW_FASTA}" \
    "${SW_FASTA}" \
    > "${SEQFILE}"

echo
echo "Starting Minigraph-Cactus 3.2.1"
echo "Reference assembly: LW"
echo
echo "Sequence file:"
cat "${SEQFILE}"
echo
echo "LW reference chromosomes:"
printf '%s\n' "${REF_CONTIGS[@]}"
echo
echo "Output directory: ${OUTDIR}"
echo "Job store: ${JOBSTORE}"
echo "Toil work directory: ${TOIL_WORKDIR}"
echo "Cactus log: ${LOGFILE}"
echo

###############################################################################
# Run Minigraph-Cactus
###############################################################################

if "${CACTUS}" \
    "${JOBSTORE}" \
    "${SEQFILE}" \
    --outDir "${OUTDIR}" \
    --outName Pchalceus_LWref \
    --reference LW \
    --refContigs "${REF_CONTIGS[@]}" \
    --otherContig chrOther \
    --gfa \
    --gbz \
    --vcf \
    --odgi \
    --maxCores "${SLURM_CPUS_PER_TASK}" \
    --maxMemory 240G \
    --mapCores 12 \
    --consCores 12 \
    --indexCores 12 \
    --workDir "${TOIL_WORKDIR}" \
    --logFile "${LOGFILE}" \
    --binariesMode singularity
then
    echo
    echo "Minigraph-Cactus completed successfully."
    echo "Results: ${OUTDIR}"
    echo "Log: ${LOGFILE}"
else
    status=$?
    echo
    echo "ERROR: Minigraph-Cactus failed with exit status ${status}." >&2
    echo "Inspect the log: ${LOGFILE}" >&2
    exit "${status}"
fi
