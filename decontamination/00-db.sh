#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name download_db
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o download_db.%j.out
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/blobtools/sorted_prim_dud/
mkdir -p nt
cd nt/

# Download files from nt.000.tar.gz to nt.192.tar.gz along with their .md5 files
for i in $(seq -w 0 192); do
    wget "https://ftp.ncbi.nlm.nih.gov/blast/db/nt.$i.tar.gz"
    wget "https://ftp.ncbi.nlm.nih.gov/blast/db/nt.$i.tar.gz.md5"
done

# Verify the integrity of each .tar.gz file using its .md5 file
for file in *.tar.gz; do
    md5file="$file.md5"
    # Extract the expected MD5 checksum from the .md5 file
    expected_md5=$(cut -d ' ' -f 1 "$md5file")
    # Calculate the MD5 checksum of the downloaded file
    actual_md5=$(md5sum "$file" | awk '{print $1}')
    
    # Compare the checksums
    if [[ "$expected_md5" == "$actual_md5" ]]; then
        echo "$file: MD5 checksum matches"
    else
        echo "$file: MD5 checksum does NOT match"
    fi
done

# Extract all downloaded .tar.gz files if they passed the MD5 check
for file in *.tar.gz; do
    tar -pxvzf "$file"
done


export BLASTDB=/scratch/leuven/357/vsc35707/blobtools/sorted_prim_dud/nt

cd /scratch/leuven/357/vsc35707/blobtools/sorted_prim_dud/
mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;
