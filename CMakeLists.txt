cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

include(CheckCCompilerFlag)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug")
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

set(SOMA_SYSTEM_INFO "${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_C_COMPILER_ID}")

#set extra flags from environment
set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} $ENV{SOMA_C_FLAGS}" CACHE STRING "" FORCE)

#Set proper release flags for PGI.
if( "${CMAKE_C_COMPILER_ID}" STREQUAL "PGI" )
  set( CMAKE_C_FLAGS_RELEASE "-O0 -DNDEBUG" CACHE STRING "" FORCE)
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
  endif(ENABLE_SANITIZE_DEBUG)
endif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")

if("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99")
endif("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")

if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xHost")
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
endif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")


#Set up options build options
if("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
  option(ENABLE_NVIDIA_GPU "Activate NVIDIA GPUs as OpenACC accelerator." ON)
  option(ENABLE_OPENMP "Activate OpenMP for shared memory system parallelism." OFF)
else("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
  option(ENABLE_NVIDIA_GPU "Activate NVIDIA GPUs as OpenACC accelerator." OFF)
  option(ENABLE_OPENMP "Activate OpenMP for shared memory system parallelism." ON)
endif("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")


option(SINGLE_PRECISION  "Single precision floats for real values. Otherwise double precision." OFF)

if(SINGLE_PRECISION)
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} SINGLE")
else()
  set(SOMA_VERSION_LONG "${SOMA_VERSION_LONG} DOUBLE")
endif()

option(ENABLE_MIC "Activate compilation with minimum image convention for bond distances" OFF)
option(ENABLE_DOMAIN_DECOMPOSITION "Activate compilation with domain decomposition support." OFF)
if( !ENABLE_MPI)
  if(ENABLE_DOMAIN_DECOMPOSITION)
    message(SEND_ERROR "Domain decomposition requires enabled MPI")
  endif(ENABLE_DOMAIN_DECOMPOSITION)
endif( !ENABLE_MPI )

if( ENABLE_NVIDIA_GPU AND ENABLE_OPENMP )
  message(SEND_ERROR "SOMA can be compiled with either OpenACC or OpenMP, not both.")
endif(ENABLE_NVIDIA_GPU AND ENABLE_OPENMP )

if(ENABLE_NVIDIA_GPU)
  find_package(CUDA REQUIRED QUIET)

  if(CUDA_VERSION_MAJOR VERSION_LESS 8)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ta=nvidia,host")
  else(CUDA_VERSION_MAJOR VERSION_LESS 8)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ta=nvidia:cuda${CUDA_VERSION_MAJOR}.${CUDA_VERSION_MINOR},host")
  endif(CUDA_VERSION_MAJOR VERSION_LESS 8)

  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Minfo=accel")

  if(ENABLE_MPI)
    # now perform some more in-depth tests of whether the MPI library supports CUDA memory
    if ( MPI_LIBRARY MATCHES mpich OR MPI_LIBRARY MATCHES "MVAPICH" )
      # find out if this is MVAPICH2
      get_filename_component(_mpi_library_dir ${MPI_LIBRARY} PATH)
      find_program(MPICH2_VERSION
	NAMES mpichversion mpich2version
	HINTS ${_mpi_library_dir} ${_mpi_library_dir}/../bin
	)
      if (MPICH2_VERSION)
	execute_process(COMMAND ${MPICH2_VERSION}
          OUTPUT_VARIABLE _output)
	if (_output MATCHES "--enable-cuda")
          set(MPI_CUDA TRUE)
          message(STATUS "Found MVAPICH2 with CUDA support.")
	endif(_output MATCHES "--enable-cuda")
      endif(MPICH2_VERSION)
    elseif(MPI_LIBRARY MATCHES libmpi)
      # find out if this is OpenMPI
      get_filename_component(_mpi_library_dir ${MPI_LIBRARY} PATH)
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
      message(STATUS "MPI<->CUDA interoperability is available and disabled by default.")
      option(ENABLE_MPI_CUDA "Enable MPI<->CUDA interoperability" off)
    else(MPI_CUDA)
      message(STATUS "MPI<->CUDA interoperability does not appear to be available.")
      option(ENABLE_MPI_CUDA "Enable MPI<->CUDA interoperability" off)
    endif(MPI_CUDA)
  endif(ENABLE_MPI)

endif(ENABLE_NVIDIA_GPU)

if(ENABLE_OPENMP)
  if("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mp")
  endif("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
  if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
  endif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopenmp-offload=mic")
  endif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
endif(ENABLE_OPENMP)


add_subdirectory(c_src)
add_subdirectory(setups)

option(INSTALL_PYTHON "Install the python scripts. Not needed for pure server installations." ON)
if(INSTALL_PYTHON)
  add_subdirectory(python_src)
endif(INSTALL_PYTHON)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  option(ENABLE_DOXYGEN "Enable the generation of API documentation." ON)
endif(DOXYGEN_FOUND)

if(ENABLE_DOXYGEN)
  add_subdirectory(documentation)
endif(ENABLE_DOXYGEN)

option(ENABLE_TESTING "Include the automated tests." ON)
if(ENABLE_TESTING)
  include(CTest)
  add_subdirectory(testing)
  enable_testing()
endif(ENABLE_TESTING)
