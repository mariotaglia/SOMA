module purge

module load git
module load cmake
module load gcc/9.1.0
module load hdf5/1.8.16+openmpi-2.0.1
module load python/cpython-3.8.5

echo "For all features an installation of h5py is needed. You may choose to install it via pip3"

export CC=gcc
