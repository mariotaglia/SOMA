module purge
module load compilers/intel
module load intelmpi
module load applications/python/3.8
module load hdf5

export CC=icc
alias cmake=cmake3
alias ccmake3=ccmake3
export HDF5_ROOT=$HDF5HOME
