#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name download_db
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH -o download_db.%j.out
#SBATCH -A lp_svbelleghem

mkdir -p nr
cd nr/

# Download files from nr.000.tar.gz to nr.118.tar.gz along with their .md5 files
for i in $(seq -w 0 118); do
    wget "https://ftp.ncbi.nlm.nih.gov/blast/db/nr.$i.tar.gz"
    wget "https://ftp.ncbi.nlm.nih.gov/blast/db/nr.$i.tar.gz.md5"
done

# Verify the integrity of each .tar.gz file using its .md5 file
for file in *.tar.gz; do
    md5file="$file.md5"
    expected_md5=$(cut -d ' ' -f 1 "$md5file")
    actual_md5=$(md5sum "$file" | awk '{print $1}')
    
    if [[ "$expected_md5" == "$actual_md5" ]]; then
        echo "$file: MD5 checksum matches"
    else
        echo "$file: MD5 checksum does NOT match"
    fi
done

# Extract all downloaded .tar.gz files
for file in *.tar.gz; do
    tar -pxvzf "$file"
done

# Set BLASTDB environment variable
export BLASTDB=./
