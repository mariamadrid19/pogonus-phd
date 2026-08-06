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

WORKDIR="/scratch/leuven/357/vsc35707/pangenome"

LW_FASTA="${WORKDIR}/Pchalceus_LW_final.fasta"
SW_FASTA="${WORKDIR}/Pchalceus_SW_final.fasta"
LW_FAI="${LW_FASTA}.fai"

CACTUS="/data/leuven/357/vsc35707/miniconda3/bin/cactus-pangenome"

SEQFILE="${WORKDIR}/Pchalceus_LWref.seqfile"
CHRFILE="${WORKDIR}/LW.chromosomes.txt"

JOBSTORE="${WORKDIR}/Pchalceus-LWref-jobstore"
OUTDIR="${WORKDIR}/Pchalceus-LWref-pangenome"
TOIL_WORKDIR="${WORKDIR}/toil-LWref-work"
LOGFILE="${WORKDIR}/Pchalceus-LWref-cactus.log"

cd "${WORKDIR}"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Working directory: ${WORKDIR}"

# Check required programs.
[[ -x "${CACTUS}" ]] || {
    echo "ERROR: cactus-pangenome is not executable: ${CACTUS}" >&2
    exit 1
}

command -v singularity >/dev/null 2>&1 || {
    echo "ERROR: singularity is unavailable." >&2
    echo "Load the appropriate Apptainer/Singularity module in this script." >&2
    exit 1
}

# Check required input files.
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
    echo "ERROR: Chromosome list not found or empty: ${CHRFILE}" >&2
    exit 1
}

# Confirm that the chromosome list contains exactly 11 entries.
CHR_COUNT=$(awk 'NF {count++} END {print count+0}' "${CHRFILE}")

if [[ "${CHR_COUNT}" -ne 11 ]]; then
    echo "ERROR: ${CHRFILE} contains ${CHR_COUNT} nonempty entries; expected 11." >&2
    exit 1
fi

# Check for duplicate chromosome identifiers.
DUPLICATES=$(
    awk 'NF {print $1}' "${CHRFILE}" |
        sort |
        uniq -d
)

if [[ -n "${DUPLICATES}" ]]; then
    echo "ERROR: Duplicate identifiers found in ${CHRFILE}:" >&2
    echo "${DUPLICATES}" >&2
    exit 1
fi

# Confirm that every chromosome identifier occurs in the LW FASTA index.
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

# Do not accidentally mix this run with an earlier workflow.
if [[ -e "${JOBSTORE}" ]]; then
    echo "ERROR: Job store already exists: ${JOBSTORE}" >&2
    echo "Use a new JOBSTORE name or intentionally restart the existing workflow." >&2
    exit 1
fi

if [[ -e "${OUTDIR}" ]]; then
    echo "ERROR: Output directory already exists: ${OUTDIR}" >&2
    echo "Use a new OUTDIR name to avoid mixing results from different runs." >&2
    exit 1
fi

# Toil requires its working directory to exist before starting.
mkdir -p "${TOIL_WORKDIR}"

# Create the two-column Cactus sequence file with LW first.
printf "LW\t%s\nSW\t%s\n" \
    "${LW_FASTA}" \
    "${SW_FASTA}" \
    > "${SEQFILE}"

# Read the chromosome identifiers into a Bash array.
mapfile -t REF_CONTIGS < <(
    awk 'NF {print $1}' "${CHRFILE}"
)

echo
echo "Starting Minigraph-Cactus"
echo "Reference assembly: LW"
echo
echo "Sequence file:"
cat "${SEQFILE}"
echo
echo "LW reference chromosomes:"
printf '%s\n' "${REF_CONTIGS[@]}"
echo

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
    --panacus \
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
    echo "Inspect: ${LOGFILE}" >&2
    exit "${status}"
fi
