#cmake_minimum_required(VERSION 3.13.4 FATAL_ERROR)
cmake_minimum_required(VERSION 3.11.0 FATAL_ERROR)

include(CheckCCompilerFlag)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug" "RelWithDebInfo")
endif()
if( ${CMAKE_VERSION} VERSION_GREATER "3.7.0" )
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
    cmake_policy(SET CMP0074 NEW)
  endif()
endif()

project(SOMA C)

# bring in custom modules
add_subdirectory (CMake)
#set version information
include(GetGitRevisionDescription)
git_describe(SOMA_GIT_VERSION)
if (SOMA_GIT_VERSION)
  set(SOMA_VERSION_LONG "${SOMA_GIT_VERSION}")
else (SOMA_GIT_VERSION)
  set(SOMA_VERSION_LONG "${SOMA_VERSION}-unknown")
endif (SOMA_GIT_VERSION)

set(SOMA_SYSTEM_INFO "${CMAKE_SYSTEM} with ${CMAKE_C_COMPILER_ID} for ${CMAKE_SYSTEM_PROCESSOR} ")

#set extra flags from environment
set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} $ENV{SOMA_C_FLAGS}" CACHE STRING "" FORCE)

#Set proper release flags for PGI.
if( "${CMAKE_C_COMPILER_ID}" STREQUAL "PGI" )
  set( CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG" CACHE STRING "" FORCE)
  set( CMAKE_C_FLAGS_RELWITHDEBINFO "-g -O0 -DNDEBUG" CACHE STRING "" FORCE)
endif()


option(ENABLE_MPI "Enable distributed memory execution via MPI" ON)
if(ENABLE_MPI)
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} MPI")
else(ENABLE_MPI)
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} NO-MPI")
endif(ENABLE_MPI)

if(ENABLE_MPI)
  #Find dependency packages
  find_package(MPI REQUIRED)
  if (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  else (MPI_FOUND)
    message(SEND_ERROR "MPI not found. Install MPI or disable MPI support.")
  endif (MPI_FOUND)
endif(ENABLE_MPI)

#Cache the HDF5_ROOT if specified by user
if (NOT "$ENV{HDF5_ROOT}" STREQUAL "")
  set(HDF5_ROOT "$ENV{HDF5_ROOT}" CACHE INTERNAL "Stored user specified HDF5_ROOT")
endif()


find_package(HDF5 REQUIRED)
if (HDF5_FOUND)
  include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
else (HDF5_FOUND)
  message(SEND_ERROR "SOMA requires a hdf5 library.")
endif(HDF5_FOUND)

if(ENABLE_MPI)
  #Test HDF5 for parallel support.
  if(NOT HDF5_IS_PARALLEL)
    message(SEND_ERROR "SOMA requires a hdf5 library with parallel IO support.")
  endif(NOT HDF5_IS_PARALLEL)
endif(ENABLE_MPI)

#Check for supported compilers.
if( NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" AND NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "PGI" )
  message("SOMA is intended to be used with either GNU GCC or PGI. Proceed on your own risk.")
endif(NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" AND NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "PGI" )


#set up the compiler flags
set(soma_link_flags "-lm -ldl -lz")


if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
  set(CMAKE_C_FLAGS_DEBUG "-g -fstack-protector-all -Wstack-protector -fstack-protector-strong")
  option(ENABLE_SANITIZE_DEBUG "Activate gcc sanitize debug plugins (undefined and address)." OFF)
  if(ENABLE_SANITIZE_DEBUG)
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=undefined -fsanitize=address")
    set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} SANITIZE_DEBUG")

  endif(ENABLE_SANITIZE_DEBUG)
endif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")

if("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
  #sorry the PGI compiler messages are utterly useless. And I don't know how to shutoff the warnings in the code.
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99 -Minform=severe")
endif("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")

if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xHost")
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
endif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")

find_package(OpenACC)
#Set up options build options
if(OpenACC_C_FOUND)
  if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    option(ENABLE_OPENACC "Activate OpenACC for  accelerators." OFF)
  else("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    option(ENABLE_OPENACC "Activate OpenACC for  accelerators." ON)
  endif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
endif(OpenACC_C_FOUND)

if(NOT ENABLE_OPENACC)
find_package(OpenMP)
endif(NOT ENABLE_OPENACC)

if(OpenMP_C_FOUND)
  option(ENABLE_OPENMP "Activate OpenMP for shared memory system parallelism." ON)
endif(OpenMP_C_FOUND)

option(SINGLE_PRECISION  "Single precision floats for real values. Otherwise double precision." OFF)

if(SINGLE_PRECISION)
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} SINGLE")
else()
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} DOUBLE")
endif()

option(ENABLE_MIC "Activate compilation with minimum image convention for bond distances" OFF)
if(ENABLE_MIC)
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} MIC")
endif(ENABLE_MIC)

option(ENABLE_DOMAIN_DECOMPOSITION "Activate compilation with domain decomposition support." OFF)
if(ENABLE_DOMAIN_DECOMPOSITION)
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} DOMAIN_DECOMPOSITION")
endif(ENABLE_DOMAIN_DECOMPOSITION)

if( !ENABLE_MPI)
  if(ENABLE_DOMAIN_DECOMPOSITION)
    message(SEND_ERROR "Domain decomposition requires enabled MPI")
  endif(ENABLE_DOMAIN_DECOMPOSITION)
endif( !ENABLE_MPI )

option(ENABLE_MONOTYPE_CONVERSIONS "Activate monomer type conversion, where monomer types are stored individually for each polymer and can be switched." OFF)
if(ENABLE_MONOTYPE_CONVERSIONS)
    set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} MONOTYPE_CONVERSIONS")
endif(ENABLE_MONOTYPE_CONVERSIONS)

if( ENABLE_OPENACC AND ENABLE_OPENMP )
  message(SEND_ERROR "SOMA can be compiled with either OpenACC or OpenMP, not both.")
endif(ENABLE_OPENACC AND ENABLE_OPENMP )

if(ENABLE_OPENACC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenACC_C_FLAGS}")
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} OpenACC")

  if("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Minfo=accel")
  endif("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
  
  if(ENABLE_MPI)
    # now perform some more in-depth tests of whether the MPI library supports CUDA memory
    if ( MPI_mpi_LIBRARY MATCHES libmpi OR MPI_mpi_LIBRARY MATCHES "MVAPICH" )
      # find out if this is MVAPICH2
      get_filename_component(_mpi_library_dir ${MPI_mpi_LIBRARY} PATH)
      find_program(MPICH2_VERSION
	NAMES mpichversion mpich2version
	HINTS ${_mpi_library_dir} ${_mpi_library_dir}/../bin
	)
      if (MPICH2_VERSION)
	execute_process(COMMAND ${MPICH2_VERSION}
          OUTPUT_VARIABLE _output)
    if (_output MATCHES "--enable-cuda" OR _output MATCHES "--enable-psp-cuda-awareness")
          set(MPI_CUDA TRUE)
          message(STATUS "Found MVAPICH2 or ParaStationMPI with CUDA support.")
    endif (_output MATCHES "--enable-cuda" OR _output MATCHES "--enable-psp-cuda-awareness")
      endif(MPICH2_VERSION)
    elseif(MPI_mpi_LIBRARY MATCHES libmpi)
      # find out if this is OpenMPI
      get_filename_component(_mpi_library_dir ${MPI_mpi_LIBRARY} PATH)
      find_program(OMPI_INFO
	NAMES ompi_info
	HINTS ${_mpi_library_dir} ${_mpi_library_dir}/../bin
	)
      if (OMPI_INFO)
	execute_process(COMMAND ${OMPI_INFO}
          OUTPUT_VARIABLE _output)
	if (_output MATCHES "smcuda")
          set(MPI_CUDA TRUE)
          message(STATUS "Found OpenMPI with CUDA support.")
	endif(_output MATCHES "smcuda")
      endif(OMPI_INFO)
    endif()


    if (MPI_CUDA)
      option(ENABLE_MPI_CUDA "Enable MPI<->CUDA interoperability" off)
      if(ENABLE_MPI_CUDA)
	message(STATUS "MPI<->CUDA interoperability is available and is enabled.")
      else(ENABLE_MPI_CUDA)
        message(STATUS "MPI<->CUDA interoperability is available but is disabled.")
      endif(ENABLE_MPI_CUDA)
    else(MPI_CUDA)
      message(STATUS "MPI<->CUDA interoperability does not appear to be available.")
      option(ENABLE_MPI_CUDA "Enable MPI<->CUDA interoperability" off)
    endif(MPI_CUDA)
endif(ENABLE_MPI)

endif(ENABLE_OPENACC)

if(ENABLE_MPI_CUDA)
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} MPI_CUDA")
  option(ENABLE_NCCL "Enable GPU communications via CUDA NCCL." off)
  if(ENABLE_NCCL)
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} CUDA_NCCL")
  endif(ENABLE_NCCL)
endif(ENABLE_MPI_CUDA)

if(ENABLE_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} OpenMP")
endif(ENABLE_OPENMP)

add_subdirectory(c_src)
add_subdirectory(setups)

find_package(Python)

if(PYTHON_FOUND)
  set(install_python_text "Install the python scripts. Not needed for pure server installations.")
  #check for numpy and h5py
  set(find_numpy "
from __future__ import print_function; import sys, os; sys.stdout = open(os.devnull, 'w');
try: import numpy;
except ImportError: print('OFF',file=sys.stderr,end='');
else: print('ON',file=sys.stderr,end='');")
  execute_process(COMMAND ${Python_EXECUTABLE} -c "${find_numpy}" ERROR_VARIABLE numpy_FOUND)
  if(${numpy_FOUND})
    message(STATUS "Numpy found")

    set(find_h5py "
from __future__ import print_function; import sys, os; #sys.stdout = open(os.devnull, 'w');
try: import h5py;
except ImportError: print('OFF',file=sys.stderr,end='');
else: print('ON',file=sys.stderr,end='');")
    execute_process(COMMAND ${Python_EXECUTABLE} -c "${find_h5py}" ERROR_VARIABLE h5py_FOUND)

    #Extra resilient h5py detection
    string(FIND ${h5py_FOUND} "ON" ON_FOUND)
    if( ON_FOUND GREATER -1)
       string(SUBSTRING ${h5py_FOUND} ${ON_FOUND} 2 h5py_FOUND)
    endif()

    if(${h5py_FOUND})
      message(STATUS "h5py found")
      option(INSTALL_PYTHON "${install_python_text}" ON)
      set(SOMA_PYTHON_SHEBANG "#!${Python_EXECUTABLE}")
    else(${h5py_FOUND})
      message(STATUS "h5py NOT found. ${h5py_FOUND}")
      option(INSTALL_PYTHON "${install_python_text}" OFF)
    endif(${h5py_FOUND})
  else(${numpy_FOUND})
    message(STATUS "numpy NOT found. (not looking for h5py anymore)")
    option(INSTALL_PYTHON "${install_python_text}" OFF)
  endif(${numpy_FOUND})
else(PYTHON_FOUND)
  message(STATUS "no viable python found.")
  option(INSTALL_PYTHON "${install_python_text}" OFF)
endif(PYTHON_FOUND)

if(INSTALL_PYTHON)
  add_subdirectory(python_src)
  option(ENABLE_TESTING "Include the automated tests." ON)

  #Note: testing is only possible if the python scripts are installed and valid h5py/numpy is found
  if(ENABLE_TESTING)
    include(CTest)
    add_subdirectory(testing)
    enable_testing()
  endif(ENABLE_TESTING)
else(INSTALL_PYTHON)
  message(WARNING "No python script installation. There may be a problem with python environment (see earlier messages). It is necessary for testing and initialization from XML. YELL: \"you probably want to fix this (unless you know what you are doing.)\"")
endif(INSTALL_PYTHON)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  option(ENABLE_DOXYGEN "Enable the generation of API documentation." ON)
endif(DOXYGEN_FOUND)

if(ENABLE_DOXYGEN)
  add_subdirectory(documentation)
endif(ENABLE_DOXYGEN)

message(STATUS "Build configured as: ${SOMA_VERSION_LONG}")
