module purge

ml Stages/2020
ml NVHPC/21.1-GCC-9.3.0
ml ParaStationMPI
ml CMake
ml NCCL
ml h5py
ml HDF5

export CUDA_VISIBLE_DEVICES=0,1,2,3


export CC=nvc
