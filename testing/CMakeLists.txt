#   Copyright (C) 2017 Ludwig Schneider
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

#copy relevant files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../example/coord.xml ${CMAKE_CURRENT_BINARY_DIR}/coord.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../example/complex.xml ${CMAKE_CURRENT_BINARY_DIR}/complex.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../example/coord.dat ${CMAKE_CURRENT_BINARY_DIR}/coord.dat COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/homo.xml ${CMAKE_CURRENT_BINARY_DIR}/homo.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/diblock.xml ${CMAKE_CURRENT_BINARY_DIR}/diblock.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scmf1.xml ${CMAKE_CURRENT_BINARY_DIR}/scmf1.xml COPYONLY)

#quick
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quick.sh.in ${CMAKE_CURRENT_BINARY_DIR}/quick.sh)
add_test(NAME quick WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./quick.sh)

#convert
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/convert.sh.in ${CMAKE_CURRENT_BINARY_DIR}/convert.sh)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compareDatH5.py.in ${CMAKE_CURRENT_BINARY_DIR}/compareDatH5.py)
add_test(NAME convert WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./convert.sh)

#xml
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/xml.sh.in ${CMAKE_CURRENT_BINARY_DIR}/xml.sh)
add_test(NAME xml WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./xml.sh)

set(ACC_FLAG "")
if(ENABLE_OPENMP)
  set(ACC_FLAG " -n ${MPIEXEC_MAX_NUMPROCS}")
endif(ENABLE_OPENMP)
if(ENABLE_NVIDIA_GPU)
  set(ACC_FLAG " -g ${MPIEXEC_MAX_NUMPROCS}")
endif(ENABLE_NVIDIA_GPU)

set(MPI_PREFIX "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}")

if(HDF5_DIFF_EXECUTABLE)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/RestartExact.sh.in ${CMAKE_CURRENT_BINARY_DIR}/RestartExact.sh)
  add_test(NAME Restart WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./RestartExact.sh)
  
  if(ENABLE_NVIDIA_GPU)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GPUexact.sh.in ${CMAKE_CURRENT_BINARY_DIR}/GPUexact.sh)
    add_test(NAME GPUexact WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./GPUexact.sh)
  endif(ENABLE_NVIDIA_GPU)

endif(HDF5_DIFF_EXECUTABLE)


#stat
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/statistics.py.in ${CMAKE_CURRENT_BINARY_DIR}/statistics.py)
add_test(NAME simple-stat WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=${ACC_FLAG})
add_test(NAME MPI-stat WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --prefix=${MPI_PREFIX} --additional-flags=${ACC_FLAG})
#foreach(PRNG PCG32 MT TT800)
foreach(PRNG PCG32 MT)
  add_test(NAME ${PRNG}-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=-p${PRNG}${ACC_FLAG})
endforeach(PRNG PCG32 MT)
#endforeach(PRNG PCG32 MT TT800)
# foreach(MOVETYPE TRIAL SMART)
#   add_test(NAME ${MOVETYPE}-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=--move-type=${MOVETYPE}${ACC_FLAG})
# endforeach(MOVETYPE TRIAL SMART)
foreach(ITERALG POLY SET)
  add_test(NAME ${ITERALG}-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --additional-flags=--iteration-alg=${ITERALG}${ACC_FLAG})
endforeach(ITERALG POLY SET)
