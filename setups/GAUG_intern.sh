#module unload pgi
module load pgi
module load mpi
export OMPI_CC=pgcc
export OMPI_CFLAGS=""
export CC=pgcc
