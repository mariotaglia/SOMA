module purge

module use $OTHERSTAGES
module load Stages/Devel-2019a
module load PGI/20.1-GCC-8.3.0
module load ParaStationMPI/5.4.6-1-CUDA
module load HDF5/1.10.5
module load CMake/3.16.5

export CC=pgcc
