#! /bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name rename2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=8:00:00
#SBATCH -o rename2.%j.out
#SBATCH -A lp_svbelleghem

# Loop through directories BAR4-01 to BAR4-10
for dir in BAR4-0{1..9} BAR4-10; do
    echo "Processing $dir..."

    # Check if the directory contains fq.gz files
    fq_files=("$dir"/*.fq.gz)
    if [ ${#fq_files[@]} -eq 0 ]; then
        echo "  No fq.gz files found in $dir. Skipping..."
        continue
    fi

    # Count the number of fq.gz files
    fq_count=${#fq_files[@]}

    # Extract sample ID (e.g., "01" from "BAR4-01")
    sample_id=$(echo "$dir" | sed 's/BAR4-//')

    # If there are 4 fq.gz files, concatenate them
    if [ "$fq_count" -eq 4 ]; then
        echo "  Found 4 fq.gz files. Concatenating..."

        # Concatenate all "_1.fq.gz" files into a single "_1.fq.gz"
        cat "$dir"/*_1.fq.gz > "$dir"/Bar4_"$sample_id"_R1.fq.gz
        cat "$dir"/*_2.fq.gz > "$dir"/Bar4_"$sample_id"_R2.fq.gz

        # Remove original 4 files after merging
        rm "$dir"/*_1.fq.gz "$dir"/*_2.fq.gz
    elif [ "$fq_count" -eq 2 ]; then
        echo "  Found 2 fq.gz files. Renaming..."

        # Rename files directly
        mv "$dir"/*_1.fq.gz "$dir"/Bar4_"$sample_id"_R1.fq.gz
        mv "$dir"/*_2.fq.gz "$dir"/Bar4_"$sample_id"_R2.fq.gz
    else
        echo "  Unexpected number of fq.gz files in $dir. Skipping..."
        continue
    fi

    echo "  Done with $dir."
done

echo "All directories processed."
