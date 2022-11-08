#   Copyright (C) 2017-2021 Ludwig Schneider
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

if(NOT DEFINED SOMA_PYTHON_SHEBANG)
  message(WARNING "An invalid python/numpy/h5py environment may prevent the execution of certain tests")
endif()

if(ENABLE_OPENACC)
  set(NTEST_GPUS 1 CACHE STRING "Number of GPUs (per node) available for testing.")
endif(ENABLE_OPENACC)

set(ACC_FLAG "")
if(ENABLE_OPENMP)
  set(ACC_FLAG " -n 2")
endif(ENABLE_OPENMP)
if(ENABLE_OPENACC)
  if(NTEST_GPUS GREATER 1)
    set(ACC_FLAG " -g ${NTEST_GPUS}")
  else(NTEST_GPUS GREATER 1)
    set(ACC_FLAG " -o 0")
  endif(NTEST_GPUS GREATER 1)
endif(ENABLE_OPENACC)


set(MPI_PREFIX "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2")
set(MPI_PREFIX_SINGLE "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2")


#copy relevant files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../example/coord.xml ${CMAKE_CURRENT_BINARY_DIR}/coord.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../example/complex.xml ${CMAKE_CURRENT_BINARY_DIR}/complex.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../example/coord.dat ${CMAKE_CURRENT_BINARY_DIR}/coord.dat COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/homo.xml ${CMAKE_CURRENT_BINARY_DIR}/homo.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/diblock.xml ${CMAKE_CURRENT_BINARY_DIR}/diblock.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scmf1.xml ${CMAKE_CURRENT_BINARY_DIR}/scmf1.xml COPYONLY)
if(ENABLE_MONOTYPE_CONVERSIONS)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/domain_monoconv.xml ${CMAKE_CURRENT_BINARY_DIR}/domain.xml COPYONLY)
else(ENABLE_MONOTYPE_CONVERSIONS)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/domain.xml ${CMAKE_CURRENT_BINARY_DIR}/domain.xml COPYONLY)
endif(ENABLE_MONOTYPE_CONVERSIONS)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mobility_high_concentration_eq.xml ${CMAKE_CURRENT_BINARY_DIR}/mobility_high_concentration_eq.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mobility_high_concentration.xml ${CMAKE_CURRENT_BINARY_DIR}/mobility_high_concentration.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mobility_low_concentration_eq.xml ${CMAKE_CURRENT_BINARY_DIR}/mobility_low_concentration_eq.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mobility_low_concentration.xml ${CMAKE_CURRENT_BINARY_DIR}/mobility_low_concentration.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/area51.xml ${CMAKE_CURRENT_BINARY_DIR}/area51.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mic.xml ${CMAKE_CURRENT_BINARY_DIR}/mic.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mic_time.xml ${CMAKE_CURRENT_BINARY_DIR}/mic_time.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/diverse_polytype.xml ${CMAKE_CURRENT_BINARY_DIR}/diverse_polytype.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/async.xml ${CMAKE_CURRENT_BINARY_DIR}/async.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/structure_factor.xml ${CMAKE_CURRENT_BINARY_DIR}/structure_factor.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/external_field.xml ${CMAKE_CURRENT_BINARY_DIR}/external_field.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/weight.xml ${CMAKE_CURRENT_BINARY_DIR}/weight.xml COPYONLY)

#quick
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quick.sh.in ${CMAKE_CURRENT_BINARY_DIR}/quick.sh)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME quick WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./quick.sh)
endif()

#Walltime
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/walltime.sh.in ${CMAKE_CURRENT_BINARY_DIR}/walltime.sh)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME walltime WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./walltime.sh)
endif()

#convert
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/convert.sh.in ${CMAKE_CURRENT_BINARY_DIR}/convert.sh)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compareDatH5.py.in ${CMAKE_CURRENT_BINARY_DIR}/compareDatH5.py)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME convert WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./convert.sh)
endif()

#xml
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/xml.sh.in ${CMAKE_CURRENT_BINARY_DIR}/xml.sh)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME xml WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./xml.sh)
endif()

#self_documentation
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/self_documentation.sh.in ${CMAKE_CURRENT_BINARY_DIR}/self_documentation.sh)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME self-documentation WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./self_documentation.sh)
endif()

if(ENABLE_MPI)
  #async
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/async.py.in ${CMAKE_CURRENT_BINARY_DIR}/async.py)
  if(DEFINED SOMA_PYTHON_SHEBANG)
    add_test(NAME async WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} async.py ${ACC_FLAG})
  endif()
endif(ENABLE_MPI)

#ConfGen_AnaGen
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/AnaGen.sh.in ${CMAKE_CURRENT_BINARY_DIR}/AnaGen.sh)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME AnaGen WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./AnaGen.sh)
endif()

#Compare Bead data even if the order of chains has changed.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compare_mixed_bead_data.py.in ${CMAKE_CURRENT_BINARY_DIR}/compare_mixed_bead_data.py)

#Structure factor
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/structure_factor.py.in ${CMAKE_CURRENT_BINARY_DIR}/structure_factor.py)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME StructureFactor WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} structure_factor.py ${ACC_FLAG})
endif()

#External field
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/external_field.py.in ${CMAKE_CURRENT_BINARY_DIR}/external_field.py)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME ExternalField WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} external_field.py)
endif()

if(HDF5_DIFF_EXECUTABLE)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/RestartExact.sh.in ${CMAKE_CURRENT_BINARY_DIR}/RestartExact.sh)
  if(ENABLE_MONOTYPE_CONVERSIONS)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../example/monoconversion.xml ${CMAKE_CURRENT_BINARY_DIR}/restartexact.xml COPYONLY)
  else(ENABLE_MONOTYPE_CONVERSIONS)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../example/complex.xml ${CMAKE_CURRENT_BINARY_DIR}/restartexact.xml COPYONLY)
  endif(ENABLE_MONOTYPE_CONVERSIONS)

  if(DEFINED SOMA_PYTHON_SHEBANG)
    add_test(NAME Restart WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./RestartExact.sh)

    if(ENABLE_OPENACC)
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GPUexact.sh.in ${CMAKE_CURRENT_BINARY_DIR}/GPUexact.sh)
      add_test(NAME GPUexact WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./GPUexact.sh)
    endif(ENABLE_OPENACC)
  endif()

endif(HDF5_DIFF_EXECUTABLE)

if(ENABLE_MPI)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/LoadBalance.sh.in ${CMAKE_CURRENT_BINARY_DIR}/LoadBalance.sh)
  if(DEFINED SOMA_PYTHON_SHEBANG)
    add_test(NAME Load-Balance WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./LoadBalance.sh)
  endif()

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tagA.xml ${CMAKE_CURRENT_BINARY_DIR}/tagA.xml COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tagB.xml ${CMAKE_CURRENT_BINARY_DIR}/tagB.xml COPYONLY)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Tag.sh.in ${CMAKE_CURRENT_BINARY_DIR}/Tag.sh)
  if(DEFINED SOMA_PYTHON_SHEBANG)
    add_test(NAME tags-and-type-switching WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./Tag.sh)
  endif()

endif(ENABLE_MPI)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/TestArea51.sh.in ${CMAKE_CURRENT_BINARY_DIR}/TestArea51.sh)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME TestArea51 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./TestArea51.sh)
endif()


#domain decomposition tests
if(ENABLE_DOMAIN_DECOMPOSITION)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/RestartDomain.sh.in ${CMAKE_CURRENT_BINARY_DIR}/RestartDomain.sh)
  if(DEFINED SOMA_PYTHON_SHEBANG)
    add_test(NAME RestartDomain WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./RestartDomain.sh)
    endif()
endif(ENABLE_DOMAIN_DECOMPOSITION)

if(ENABLE_MIC)
  #mic
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mic.py.in ${CMAKE_CURRENT_BINARY_DIR}/mic.py)
  if(DEFINED SOMA_PYTHON_SHEBANG)
    add_test(NAME mic WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND  ${Python_EXECUTABLE} mic.py)
  endif()
endif(ENABLE_MIC)

if(ENABLE_MONOTYPE_CONVERSIONS)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/monoconversion.xml ${CMAKE_CURRENT_BINARY_DIR}/monoconversion.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/monoconversion_dd.xml ${CMAKE_CURRENT_BINARY_DIR}/monoconversion_dd.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/monoconversion_test.py.in ${CMAKE_CURRENT_BINARY_DIR}/monoconversion_test.py)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME MonotypeConversion WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} monoconversion_test.py --additional-flags=${ACC_FLAG})
endif()

else(ENABLE_MONOTYPE_CONVERSIONS)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conversion.xml ${CMAKE_CURRENT_BINARY_DIR}/conversion.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conversion_test.py.in ${CMAKE_CURRENT_BINARY_DIR}/conversion_test.py)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME PolytypeConversion WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} conversion_test.py --additional-flags=${ACC_FLAG})
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/partial_conversion.xml ${CMAKE_CURRENT_BINARY_DIR}/partial_conversion.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/partial_conversion_dd.xml ${CMAKE_CURRENT_BINARY_DIR}/partial_conversion_dd.xml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/partial_conversion_test.py.in ${CMAKE_CURRENT_BINARY_DIR}/partial_conversion_test.py)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME PolytypePartialConversion WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} partial_conversion_test.py --additional-flags=${ACC_FLAG})
endif()

endif(ENABLE_MONOTYPE_CONVERSIONS)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/weight.py.in ${CMAKE_CURRENT_BINARY_DIR}/weight.py)
if(DEFINED SOMA_PYTHON_SHEBANG)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_mobility.py.in ${CMAKE_CURRENT_BINARY_DIR}/test_mobility.py)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mobility.sh.in ${CMAKE_CURRENT_BINARY_DIR}/mobility.sh)
  add_test(NAME mobility-test WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ./mobility.sh)

  add_test(NAME weight-test WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} weight.py --prefix=${MPI_PREFIX_SINGLE} --additional-flags=${ACC_FLAG})
endif(DEFINED SOMA_PYTHON_SHEBANG)

#stat
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/statistics.py.in ${CMAKE_CURRENT_BINARY_DIR}/statistics.py)
if(DEFINED SOMA_PYTHON_SHEBANG)
  add_test(NAME default-stat WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} statistics.py --prefix=${MPI_PREFIX_SINGLE} --additional-flags=${ACC_FLAG})#default
  if(ENABLE_MPI)
    add_test(NAME MPI-stat WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} statistics.py --prefix=${MPI_PREFIX} --additional-flags=${ACC_FLAG})
  endif(ENABLE_MPI)
  #add_test(NAME PCG32-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND python statistics.py --prefix=${MPI_PREFIX_SINGLE} --additional-flags=-pPCG32${ACC_FLAG}) #default case
  add_test(NAME MT-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} statistics.py --prefix=${MPI_PREFIX_SINGLE} --additional-flags=-pMT${ACC_FLAG})
  if(ENABLE_MIC)
    add_test(NAME MIC-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} statistics.py --prefix=${MPI_PREFIX_SINGLE} --additional-flags=--bond-minimum-image-convention${ACC_FLAG})
  endif(ENABLE_MIC)
  #add_test(NAME TT800-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} statistics.py --additional-flags=-pTT800${ACC_FLAG}) #not working prng
  #add_test(NAME SMART-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} statistics.py --additional-flags=--move-type=SMART${ACC_FLAG}) #default
  #add_test(NAME TRIAL-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} statistics.py --additional-flags=--move-type=TRIAL${ACC_FLAG}) #not working move type (possibly because of bad statistics)
  #add_test(NAME POLYMER-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} statistics.py --prefix=${MPI_PREFIX_SINGLE} --additional-flags=--iteration-alg=POLYMER${ACC_FLAG}) #default
  add_test(NAME SET-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} statistics.py --prefix=${MPI_PREFIX_SINGLE} --additional-flags=--iteration-alg=SET${ACC_FLAG})
  #add_test(NAME SET-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} statistics.py --additional-flags=--iteration-alg=SET--set-generation-algorithm=SIMPLE${ACC_FLAG}) #default
  add_test(NAME SET-FIXED-stat- WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} statistics.py --prefix=${MPI_PREFIX_SINGLE} "--additional-flags=--iteration-alg=SET --set-generation-algorithm=FIXED${ACC_FLAG}") #default

#domain stat test
if(ENABLE_DOMAIN_DECOMPOSITION)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/domain_statistics.py.in ${CMAKE_CURRENT_BINARY_DIR}/domain_statistics.py)

  add_test(NAME domain-statistics-test WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python_EXECUTABLE} domain_statistics.py --additional-flags=${ACC_FLAG})


endif(ENABLE_DOMAIN_DECOMPOSITION)
endif()
