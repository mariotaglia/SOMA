module --force purge
ml Stages/2022
ml NVHPC
ml ParaStationMPI
ml HDF5
ml CMake
ml NCCL
ml h5py
ml HDF5
export CUDA_VISIBLE_DEVICES=0,1,2,3
export CC=nvc
