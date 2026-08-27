#!/usr/bin/env bash
#SBATCH --cluster=wice
#SBATCH --job-name=extract_pogonus
#SBATCH --output=extract_pogonus_%j.out
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -A lp_svbelleghem

INPUT_DIR="$PWD"
OUTPUT_DIR="$INPUT_DIR/241213_Pogonus_PacBio"

ZIP_FILES=(
  "1737128976_241213_Pogonus_PacBio_241213_Pogonus_PacBio_part_1_of_2.zip"
  "1737133348_241213_Pogonus_PacBio_241213_Pogonus_PacBio_part_2_of_2.zip"
)

mkdir -p "$OUTPUT_DIR"

for zip_file in "${ZIP_FILES[@]}"; do
    zip_path="$INPUT_DIR/$zip_file"

    if [[ ! -f "$zip_path" ]]; then
        echo "ERROR: File not found: $zip_path" >&2
        exit 1
    fi

    echo "Checking $zip_file..."
    unzip -tq "$zip_path"

    echo "Extracting $zip_file..."
    unzip -q -o "$zip_path" -d "$OUTPUT_DIR"
done

echo "Extraction completed successfully."
echo "Output directory: $OUTPUT_DIR"
