#!/bin/bash -l
#SBATCH --job-name=dl_eggnog_db
#SBATCH --cluster=genius
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/download_eggnog_%j.out
#SBATCH -e /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/download_eggnog_%j.err

set -euo pipefail

PROJECT=/scratch/leuven/357/vsc35707/annotation/func-annotation
DBDIR=$PROJECT/databases/eggnog

mkdir -p "$DBDIR"
cd "$DBDIR"

echo "======================================="
echo "eggNOG DB download job"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Target dir: $DBDIR"
echo "======================================="

download_if_missing () {
    local outfile="$1"
    local url="$2"

    if [[ -s "$outfile" ]]; then
        echo "[SKIP] $outfile already exists"
    else
        echo "[GET ] $outfile"
        wget -c -O "$outfile" "$url"
    fi
}

echo ""
echo "Step 1/3: eggnog.db.gz"
download_if_missing "eggnog.db.gz" \
    "http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz"

if [[ -s "eggnog.db.gz" && ! -s "eggnog.db" ]]; then
    echo "[UNZIP] eggnog.db.gz"
    gunzip -f eggnog.db.gz
fi

echo ""
echo "Step 2/3: eggnog.taxa.tar.gz"
download_if_missing "eggnog.taxa.tar.gz" \
    "http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz"

if [[ -s "eggnog.taxa.tar.gz" && ! -s "eggnog.taxa.db" ]]; then
    echo "[UNTAR] eggnog.taxa.tar.gz"
    tar -xzf eggnog.taxa.tar.gz
fi

echo ""
echo "Step 3/3: eggnog_proteins.dmnd.gz"
download_if_missing "eggnog_proteins.dmnd.gz" \
    "http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz"

if [[ -s "eggnog_proteins.dmnd.gz" && ! -s "eggnog_proteins.dmnd" ]]; then
    echo "[UNZIP] eggnog_proteins.dmnd.gz"
    gunzip -f eggnog_proteins.dmnd.gz
fi

echo ""
echo "Cleaning compressed archives if decompression succeeded..."
[[ -s "eggnog.db" ]] && rm -f eggnog.db.gz || true
[[ -s "eggnog.taxa.db" ]] && rm -f eggnog.taxa.tar.gz || true
[[ -s "eggnog_proteins.dmnd" ]] && rm -f eggnog_proteins.dmnd.gz || true

echo ""
echo "Final contents:"
ls -lh "$DBDIR"

echo ""
echo "Disk usage:"
du -sh "$DBDIR"

echo ""
echo "Done: $(date)"
