module purge

ml use $OTHERSTAGES
ml Stages/2020
ml NVHPC/20.7-GCC-9.3.0
ml ParaStationMPI/5.4.7-1
ml HDF5/1.10.6
ml CMake/3.18.0

export CC=nvc
