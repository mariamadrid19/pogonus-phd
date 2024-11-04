# Set name of directory in which config.yaml file can be found
# The file should be available at /path/to/datasets/$ACCESSION/config.yaml
ACCESSION=sorted_prim_dud
# Set maximum number of threads to use
THREADS=24
# Set TRANSFER=true to remove intermediate files and place results in a separate directory
# at /path/to/output/<PREFIX>, where PREFIX is taken from assembly.prefix in config.yaml
# The final BlobDir will be available as /path/to/output/$PREFIX/PREFIX.tar
TRANSFER=true
docker run --rm \
    --name btk-$ACCESSION \
    -e ASSEMBLY=$ACCESSION \
    -e THREADS=$THREADS \
    -e TRANSFER=$TRANSFER \
    -v /scratch/leuven/357/vsc35707/blobtools:/blobtoolkit/datasets \
    -v /scratch/leuven/357/vsc35707/blobtools:/blobtoolkit/databases \
    -v /scratch/leuven/357/vsc35707/blobtools:/blobtoolkit/output \
    genomehubs/blobtoolkit:$RELEASE
