module purge
module load sw.skl
module load slurm
module load HLRNenv
module load gcc
module load openmpi/3.1.2-gcc
module load hdf5/1.10.4_gcc.8.2_openmpi.3.1.2
module load cmake

#HLRN fails to set this correctly automatically, so do it manually
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LD_RUN_PATH
export CC=gcc
