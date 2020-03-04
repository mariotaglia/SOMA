module purge

module load CUDA
module load PGI
module load MVAPICH2
module load CMake
module load HDF5

#export SOMA_C_FLAGS="-Mnollvm -L$EBROOTNUMACTL/lib/ -L$EBROOTGCCCORE/lib64/ -L$EBROOTHDF5/lib/ -L$EBROOTMVAPICH2/lib64/ -L$EBROOTSZIP/lib/ -L$EBROOTCUDA/lib64/"

export MV2_USE_GPUDIRECT_GDRCOPY=0

export CC=pgcc
