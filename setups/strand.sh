module purge
module load compilers/intel
module load intelmpi
module load applications/utils
module load applications/anaconda3
module load hdf5

export CC=icc
#export CMAKE_XTRA_HDF5="-DHDF5_CXX_COMPILER_EXECUTABLE=$HDF5HOME/bin/h5c++ -DHDF5_C_COMPILER_EXECUTABLE=$HDF5HOME/bin/h5cc -DHDF5_C_INCLUDE_DIR=$HDF5HOME/include -DHDF5_DIFF_EXECUTABLE=$HDF5HOME/bin/h5diff -DHDF5_hdf5_LIBRARY_RELEASE=$HDF5HOME/lib/libhdf5.a -DHDF5_DIR=$HDF5HOME/"
export CMAKE_XTRA_HDF5="-DHDF5_hdf5_LIBRARY_RELEASE=$HDF5HOME/lib/libhdf5.a -DHDF5_C_INCLUDE_DIR=$HDF5HOME/include"
