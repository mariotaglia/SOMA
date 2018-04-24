module --force purge
module use /usr/local/software/jureca/OtherStages
module load Stages/2017a
module load PGI
module load MVAPICH2
module load CMake
module load HDF5
export CC=pgcc
