module purge

ml use $OTHERSTAGES
ml Stages/2020
ml NVHPC
ml ParaStationMPI
ml mpi-settings/CUDA
ml HDF5
ml CMake
ml NCCL
ml h5py

export CC=nvc
