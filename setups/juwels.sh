module purge


ml use $OTHERSTAGES
ml Stages/2020
ml NVHPC/21.1-GCC-9.3.0
ml ParaStationMPI
ml mpi-settings/CUDA
ml CMake
ml NCCL
ml h5py
ml HDF5

export CC=nvc
