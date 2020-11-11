module purge

module load CUDA
module load NVHPC
module load ParaStationMPI
module load CMake
module load HDF5
module load mpi-settings/CUDA

#export SOMA_C_FLAGS="-Mnollvm -L$EBROOTNUMACTL/lib/ -L$EBROOTGCCCORE/lib64/ -L$EBROOTHDF5/lib/ -L$EBROOTMVAPICH2/lib64/ -L$EBROOTSZIP/lib/ -L$EBROOTCUDA/lib64/"


export CC=nvc
