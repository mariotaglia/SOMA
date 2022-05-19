## It is necessary to change the omp compile flags to:
## -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a
module load GCC
module load CMake
module load OpenMPI
module load HDF5
export CC=/opt/rocm/bin/amdclang
