if(ENABLE_MPI_CUDA AND MPI_CUDA)
  add_definitions(-DENABLE_MPI_CUDA)
endif()
if(ENABLE_NCCL)
  add_definitions(-DENABLE_NCCL)
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
  io_old.c
  mc.c
  mesh.c
  monomer.c
  mpiroutines.c
  phase.c
  polymer.c
  poly_heavy.c
  rng.c
  self_documentation.c
  rng_alternative.c
  soma_config.c
  soma_memory.c
  soma_util.c
  walltime.c
  test.c
  polytype_conversion.c
  monotype_conversion.c
  mobility.c
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
  io_old.h
  mc.h
  mesh.h
  monomer.h
  mpiroutines.h
  phase.h
  polymer.h
  poly_heavy.h
  rng.h
  self_documentation.h
  rng_alternative.h
  soma_memory.h
  soma_util.h
  walltime.h
  test.h
  polytype_conversion.h
  monotype_conversion.h
  mobility.h
  )

add_library(soma_lib ${soma_source})

add_executable(SOMA soma.c)
add_executable(CONVERT convert.c)

option(BUILD_STATIC "Try to build static if possible. No guarantee your EXE is completely static" OFF)

foreach(EXE SOMA CONVERT)
  target_link_libraries(${EXE} soma_lib)
  target_link_libraries(${EXE} ${MPI_LIBRARIES})
  target_link_libraries(${EXE} ${HDF5_LIBRARIES})
  target_link_libraries(${EXE} ${HDF5_C_LIBRARIES})
  target_link_libraries(${EXE} ${soma_link_flags})

  if(BUILD_STATIC)
    if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
      target_link_libraries(${EXE} "-static -static-libgcc")
    endif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    if("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
      target_link_libraries(${EXE} "-Bstatic_pgi")
    endif("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
  endif(BUILD_STATIC)

  if(ENABLE_NCCL)
    target_link_libraries(${EXE} "-lnccl")
  endif(ENABLE_NCCL)

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


#START STRING_PLUGIN
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/string_plugin)
  option(ENABLE_STRING_PLUGIN "Running multiple SOMA-instances for the application of the string-method." off)
  if(ENABLE_STRING_PLUGIN)
    add_subdirectory(string_plugin)
  endif(ENABLE_STRING_PLUGIN)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/string_plugin)
#END STRING_PLUGIN

#START CROSSLINK_PLUGIN
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/crosslink)
  option(ENABLE_CROSSLINK_PLUGIN "Plugin for crosslinking." off)
  if(ENABLE_CROSSLINK_PLUGIN)
    add_subdirectory(crosslink)
  endif(ENABLE_CROSSLINK_PLUGIN)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/crosslink)
#END CROSSLINK_PLUGIN
