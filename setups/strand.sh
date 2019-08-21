module purge
module load compilers/intel
module load intelmpi
module load applications/utils
module load applications/anaconda3
module load hdf5

export CC=icc
#export CMAKE_XTRA_HDF5="-DHDF5_CXX_COMPILER_EXECUTABLE=$HDF5HOME/bin/h5c++ -DHDF5_C_COMPILER_EXECUTABLE=$HDF5HOME/bin/h5cc -DHDF5_C_INCLUDE_DIR=$HDF5HOME/include -DHDF5_DIFF_EXECUTABLE=$HDF5HOME/bin/h5diff -DHDF5_hdf5_LIBRARY_RELEASE=$HDF5HOME/lib/libhdf5.a -DHDF5_DIR=$HDF5HOME/"
#export CMAKE_XTRA_HDF5="-DHDF5_hdf5_LIBRARY_RELEASE=$HDF5HOME/lib/libhdf5.a -DHDF5_C_INCLUDE_DIR=$HDF5HOME/include"

#Use cmake3 over cmake2 to use the ROOT way to find the correct HDF5 library
alias cmake=cmake3
alias ccmake3=ccmake3
export HDF5_ROOT=$HDF5HOME
