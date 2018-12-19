#   Copyright (C) 2017-2018 Ludwig Schneider
#
# This file is part of SOMA.
#
# SOMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SOMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SOMA.  If not, see <http://www.gnu.org/licenses/>.

set(SOMA_PYTHON_DIR "\"${CMAKE_CURRENT_BINARY_DIR}/../python_src\"")

set(ACC_FLAG "")
if(ENABLE_OPENMP)
  set(ACC_FLAG " -n 2")
endif(ENABLE_OPENMP)
if(ENABLE_NVIDIA_GPU)
  if(ENABLE_MPI)
    set(ACC_FLAG " -g ${MPIEXEC_MAX_NUMPROCS}")
  else(ENABLE_MPI)
    set(ACC_FLAG " -o 0")
  endif(ENABLE_MPI)
endif(ENABLE_NVIDIA_GPU)

if(ENABLE_MPI)
  set(MPI_PREFIX "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2")
else(ENABLE_MPI)
    set(MPI_PREFIX "")
endif(ENABLE_MPI)


#copy relevant files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../example/coord.xml ${CMAKE_CURRENT_BINARY_DIR}/coord.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../example/complex.xml ${CMAKE_CURRENT_BINARY_DIR}/complex.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../example/coord.dat ${CMAKE_CURRENT_BINARY_DIR}/coord.dat COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/homo.xml ${CMAKE_CURRENT_BINARY_DIR}/homo.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/diblock.xml ${CMAKE_CURRENT_BINARY_DIR}/diblock.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scmf1.xml ${CMAKE_CURRENT_BINARY_DIR}/scmf1.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/domain.xml ${CMAKE_CURRENT_BINARY_DIR}/domain.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mic.xml ${CMAKE_CURRENT_BINARY_DIR}/mic.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mic_time.xml ${CMAKE_CURRENT_BINARY_DIR}/mic_time.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/diverse_polytype.xml ${CMAKE_CURRENT_BINARY_DIR}/diverse_polytype.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/async.xml ${CMAKE_CURRENT_BINARY_DIR}/async.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/structure_factor.xml ${CMAKE_CURRENT_BINARY_DIR}/structure_factor.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/external_field.xml ${CMAKE_CURRENT_BINARY_DIR}/external_field.xml COPYONLY)

#quick
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quick.sh.in ${CMAKE_CURRENT_BINARY_DIR}/quick.sh)
add_test(NAME quick WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./quick.sh)

#Walltime
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/walltime.sh.in ${CMAKE_CURRENT_BINARY_DIR}/walltime.sh)
add_test(NAME walltime WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./walltime.sh)

#convert
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/convert.sh.in ${CMAKE_CURRENT_BINARY_DIR}/convert.sh)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compareDatH5.py.in ${CMAKE_CURRENT_BINARY_DIR}/compareDatH5.py)
add_test(NAME convert WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./convert.sh)

#xml
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/xml.sh.in ${CMAKE_CURRENT_BINARY_DIR}/xml.sh)
add_test(NAME xml WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./xml.sh)

if(ENABLE_MPI)
  #async
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/async.py.in ${CMAKE_CURRENT_BINARY_DIR}/async.py)
  add_test(NAME async WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python async.py ${ACC_FLAG})
endif(ENABLE_MPI)

#ConfGen_AnaGen
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/AnaGen.sh.in ${CMAKE_CURRENT_BINARY_DIR}/AnaGen.sh)
add_test(NAME AnaGen WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./AnaGen.sh)

#Compare Bead data even if the order of chains has changed.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compare_mixed_bead_data.py.in ${CMAKE_CURRENT_BINARY_DIR}/compare_mixed_bead_data.py)

#Structure factor
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/structure_factor.py.in ${CMAKE_CURRENT_BINARY_DIR}/structure_factor.py)
add_test(NAME StructureFactor WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python structure_factor.py ${ACC_FLAG})

#External field
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/external_field.py.in ${CMAKE_CURRENT_BINARY_DIR}/external_field.py)
add_test(NAME ExternalField WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python external_field.py)

if(HDF5_DIFF_EXECUTABLE)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/RestartExact.sh.in ${CMAKE_CURRENT_BINARY_DIR}/RestartExact.sh)
  add_test(NAME Restart WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./RestartExact.sh)

  if(ENABLE_NVIDIA_GPU)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GPUexact.sh.in ${CMAKE_CURRENT_BINARY_DIR}/GPUexact.sh)
    add_test(NAME GPUexact WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./GPUexact.sh)
  endif(ENABLE_NVIDIA_GPU)

endif(HDF5_DIFF_EXECUTABLE)

if(ENABLE_MPI)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/LoadBalance.sh.in ${CMAKE_CURRENT_BINARY_DIR}/LoadBalance.sh)
  add_test(NAME Load-Balance WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./LoadBalance.sh)
endif(ENABLE_MPI)


#domain decomposition tests
if(ENABLE_DOMAIN_DECOMPOSITION)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/RestartDomain.sh.in ${CMAKE_CURRENT_BINARY_DIR}/RestartDomain.sh)
  add_test(NAME RestartDomain WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./RestartDomain.sh)
endif(ENABLE_DOMAIN_DECOMPOSITION)

if(ENABLE_MIC)
  #mic
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mic.py.in ${CMAKE_CURRENT_BINARY_DIR}/mic.py)
  add_test(NAME mic WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python mic.py)
endif(ENABLE_MIC)


#stat
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/statistics.py.in ${CMAKE_CURRENT_BINARY_DIR}/statistics.py)
add_test(NAME default-stat WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=${ACC_FLAG})#default
if(ENABLE_MPI)
  add_test(NAME MPI-stat WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --prefix=${MPI_PREFIX} --additional-flags=${ACC_FLAG})
endif(ENABLE_MPI)
#add_test(NAME PCG32-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=-pPCG32${ACC_FLAG}) #default case
add_test(NAME MT-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=-pMT${ACC_FLAG})
if(ENABLE_MIC)
  add_test(NAME MIC-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=--bond-minimum-image-convention${ACC_FLAG})
endif(ENABLE_MIC)
#add_test(NAME TT800-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=-pTT800${ACC_FLAG}) #not working prng
#add_test(NAME SMART-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=--move-type=SMART${ACC_FLAG}) #default
#add_test(NAME TRIAL-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=--move-type=TRIAL${ACC_FLAG}) #not working move type (possibly because of bad statistics)
#add_test(NAME POLYMER-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=--iteration-alg=POLYMER${ACC_FLAG}) #default
add_test(NAME SET-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=--iteration-alg=SET${ACC_FLAG})
#add_test(NAME SET-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=--iteration-alg=SET--set-generation-algorithm=SIMPLE${ACC_FLAG}) #default
add_test(NAME SET-FIXED-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py "--additional-flags=--iteration-alg=SET --set-generation-algorithm=FIXED${ACC_FLAG}") #default

#domain stat test
if(ENABLE_DOMAIN_DECOMPOSITION)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/domain_statistics.py.in ${CMAKE_CURRENT_BINARY_DIR}/domain_statistics.py)
  add_test(NAME domain-statistics-test WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python domain_statistics.py --additional-flags=${ACC_FLAG})
endif(ENABLE_DOMAIN_DECOMPOSITION)
