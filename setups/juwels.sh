module purge
ml use /gpfs/software/juwels/otherstages
ml Stages/2018b

module load CUDA
module load PGI
module load MVAPICH2
module load CMake
module load HDF5


export CC=pgcc
