module purge

module use $OTHERSTAGES
module load Stages/Devel-2019a
module load PGI/20.1-GCC-8.3.0 ParaStationMPI/5.4.6-1-CUDA HDF5
module load GCCcore/.8.3.0 CMake/3.16.5

#export SOMA_C_FLAGS="-Mnollvm -L$EBROOTNUMACTL/lib/ -L$EBROOTGCCCORE/lib64/ -L$EBROOTHDF5/lib/ -L$EBROOTMVAPICH2/lib64/ -L$EBROOTSZIP/lib/ -L$EBROOTCUDA/lib64/"

export MV2_USE_GPUDIRECT_GDRCOPY=0

export CC=pgcc
