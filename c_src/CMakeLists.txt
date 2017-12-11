if(ENABLE_MPI_CUDA AND MPI_CUDA)
  add_definitions(-DENABLE_MPI_CUDA)
endif()

set(soma_source
  ana.c
  ana_info.c
  autotuner.c
  bond.c
  cmdline.c
  generate_positions.c
  independent_sets.c
  init.c
  io.c
  mc.c
  mesh.c
  monomer.c
  mpiroutines.c
  phase.c
  polymer.c
  rng.c
  signal.c
  soma_config.c
  soma_util.c
  test.c
  )

configure_file(soma_config.h.in soma_config.h)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/soma_config.h DESTINATION include)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(soma_header
  ana.h
  ana_info.h
  autotuner.h
  bond.h
  cmdline.h
  generate_positions.h
  independent_sets.h
  init.h
  io.h
  mc.h
  mesh.h
  monomer.h
  mpiroutines.h
  phase.h
  polymer.h
  rng.h
  signal.h
  soma_util.h
  test.h
  )

add_library(soma_lib ${soma_source})

add_executable(SOMA soma.c)
add_executable(CONVERT convert.c)

foreach(EXE SOMA CONVERT)
  target_link_libraries(${EXE} soma_lib)
  target_link_libraries(${EXE} ${MPI_LIBRARIES})
  target_link_libraries(${EXE} ${HDF5_LIBRARIES})
  target_link_libraries(${EXE} ${HDF5_C_LIBRARIES})
  target_link_libraries(${EXE} ${soma_link_flags})

  install(TARGETS ${EXE} DESTINATION bin)
endforeach(EXE)

foreach(HEADER ${soma_header})
  install(FILES ${HEADER} DESTINATION include)
endforeach(HEADER ${soma_header})

#SECTION FOR PLUGINS

#START INSTANTANEOUS_SHEAR_PLUGIN
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/instantaneous_shear_plugin)
  option(ENABLE_INSTANTANEOUS_SHEAR_PLUGIN "Compilation and installation of c-program to apply an instantaneous shear to a configuration." off)
  if(ENABLE_INSTANTANEOUS_SHEAR_PLUGIN)
    add_subdirectory(instantaneous_shear_plugin)
  endif(ENABLE_INSTANTANEOUS_SHEAR_PLUGIN)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/instantaneous_shear_plugin)
#END INSTANTANEOUS_SHEAR_PLUGIN

#START SPLIT_PLUGIN
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/split_hdf5)
  option(ENABLE_SPLIT_HDF5_PLUGIN "Compilation and installation of c-program that splits a dataset to multiple files." off)
  if(ENABLE_SPLIT_HDF5_PLUGIN)
    add_subdirectory(split_hdf5)
  endif(ENABLE_SPLIT_HDF5_PLUGIN)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/split_hdf5)
#END SPLIT_PLUGIN

#START DOUBLE_PLUGIN
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/double_plugin)
  option(ENABLE_DOUBLE_PLUGIN "Running two SOMA-instances on 4 ranks with 2 ranks assigned to each instance." off)
  if(ENABLE_DOUBLE_PLUGIN)
    add_subdirectory(double_plugin)
  endif(ENABLE_DOUBLE_PLUGIN)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/double_plugin)
#END DOUBLE_PLUGIN

