module purge

ml NVHPC
ml ParaStationMPI
ml HDF5
ml CMake
ml NCCL
ml h5py

export CUDA_VISIBLE_DEVICES=0,1,2,3

export CC=nvc
