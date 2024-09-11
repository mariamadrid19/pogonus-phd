#install braker

mamba create -n braker anaconda::perl anaconda::biopython bioconda::perl-app-cpanminus bioconda::perl-file-spec bioconda::perl-hash-merge bioconda::perl-list-util bioconda::perl-module-load-conditional bioconda::perl-posix bioconda::perl-file-homedir bioconda::perl-parallel-forkmanager bioconda::perl-scalar-util-numeric bioconda::perl-yaml bioconda::perl-class-data-inheritable bioconda::perl-exception-class bioconda::perl-test-pod bioconda::perl-file-which bioconda::perl-mce bioconda::perl-threaded bioconda::perl-list-util bioconda::perl-math-utils bioconda::cdbtools eumetsat::perl-yaml-xs bioconda::perl-data-dumper

mamba activate braker

#GeneMark-ETP (Data directory)
#AUGUSTUS
module load AUGUSTUS/3.5.0-foss-2022a

#Pyhton3
#BAMTOOLS
module load BamTools/2.5.2-GCC-11.3.0

#DIAMOND
module load DIAMOND/0.9.24-foss-2018a

#STRINGTIE2
#BEDTools
#GffRead
