module purge
ml use /gpfs/software/juwels/otherstages
ml Stages/2018b

module load CUDA
module load PGI
module load MVAPICH2
module load CMake
module load HDF5


export CC=pgcc


#For installation it may be necessary to set the following environment:
# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$EBROOTCUDA/lib64/stubs"
# Do not use for running simulations
