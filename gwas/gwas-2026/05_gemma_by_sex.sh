#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=gwas_sexsplit
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o gwas_sexsplit.%j.out

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh

VCF="final-vcfs/Pchal_Bar_SW.filtered.multiSplit.renamed.vcf.gz"
PHENO="phenotypes_with_behavior_sex.txt"

TRAITS=("relMRWS" "WL" "EL" "immediate3trials" "finalPosition3trials")
SEXES=("females" "males")

# Make keep files: FID IID
awk 'BEGIN{FS=OFS="\t"} NR>1 && $15=="F" {print $1,$2}' "$PHENO" > females.keep
awk 'BEGIN{FS=OFS="\t"} NR>1 && $15=="M" {print $1,$2}' "$PHENO" > males.keep

module load PLINK/1.9

for SEX in "${SEXES[@]}"; do

    KEEP="${SEX}.keep"

    for TRAIT in "${TRAITS[@]}"; do

        PREFIX="gwas_${SEX}_${TRAIT}"
        OUTDIR="gemma-results-${SEX}-${TRAIT}"
        KDIR="kinship_matrix_${SEX}_${TRAIT}"

        echo "====================================="
        echo "Running $TRAIT for $SEX"
        echo "====================================="

        conda activate variant_tools

        plink \
          --vcf "$VCF" \
          --keep "$KEEP" \
          --pheno "$PHENO" \
          --pheno-name "$TRAIT" \
          --double-id \
          --make-bed \
          --allow-extra-chr \
          --allow-no-sex \
          --out "$PREFIX"

        plink \
          --bfile "$PREFIX" \
          --missing \
          --allow-extra-chr \
          --allow-no-sex \
          --out "missing_${SEX}_${TRAIT}"

        conda activate gemma

        mkdir -p "$KDIR" "$OUTDIR"

        gemma \
          -bfile "$PREFIX" \
          -gk 1 \
          -outdir "$KDIR" \
          -o "$PREFIX"

        gemma \
          -bfile "$PREFIX" \
          -k "$KDIR/${PREFIX}.cXX.txt" \
          -lmm 4 \
          -outdir "$OUTDIR" \
          -o "gemma_lmm_results_${SEX}_${TRAIT}"

    done
done
