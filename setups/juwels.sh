module purge
# ml use /gpfs/software/juwels/otherstages
# ml Stages/2018b

module load CUDA
module load PGI
module load MVAPICH2/2.3.2-GDR
module load CMake
module load HDF5

#export SOMA_C_FLAGS="-Mnollvm -L$EBROOTNUMACTL/lib/ -L$EBROOTGCCCORE/lib64/ -L$EBROOTHDF5/lib/ -L$EBROOTMVAPICH2/lib64/ -L$EBROOTSZIP/lib/ -L$EBROOTCUDA/lib64/"

export export LD_PRELOAD=$MV2_PATH/lib64/libmpi.so

export CC=pgcc
