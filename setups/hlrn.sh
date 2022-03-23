module purge
module load sw.skl
module load slurm
module load HLRNenv
module load gcc
module load openmpi/gcc.9/3.1.5
module load hdf5-parallel/ompi/gcc.9/1.10.6
module load cmake
module load anaconda3

#HLRN fails to set this correctly automatically, so do it manually
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LD_RUN_PATH
export CC=gcc
