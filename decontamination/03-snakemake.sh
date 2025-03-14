#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name run_smk
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o smk.%j.out
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/blobtools/blobtoolkit/insdc-pipeline/

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh

conda activate btk_env

# The working directory
WORKDIR=/scratch/leuven/357/vsc35707/blobtools/sorted_prim_dud       
# The assembly prefix
ASSEMBLY=sorted_prim_dud
# A directory to contain the Conda environments for individual Snakemake rules
CONDA_DIR=/scratch/leuven/357/vsc35707/blobtools/blobtoolkit/.conda
# The maximum number of parallel threads to run
THREADS=32
# The tool being used
TOOL=blobtoolkit
# The directory with the snakefiles
SNAKE_DIR=/scratch/leuven/357/vsc35707/blobtools/blobtoolkit/insdc-pipeline

snakemake -p \
          --use-conda \
          --conda-prefix $CONDA_DIR \
          --directory $WORKDIR/ \
          --configfile $WORKDIR/config.yaml \
          --stats $ASSEMBLY.snakemake.stats \
          -j $THREADS \
          -s $SNAKE_DIR/$TOOL.smk \
          --resources btk=1
