
include_directories( "${SOMA_SOURCE_DIR}/c_src/" )

set(crosslink_source
  #add .c files if any
  readIn_lib.c
  separateMolecule_lib.c
  perkolation_lib.c
  )
set(crosslink_header
  #add .h files if any
  readIn_lib.h
  separateMolecule_lib.h
  perkolation_lib.h
  )

add_library(crosslink_lib ${crosslink_source})

add_executable(readIn readIn.c)

foreach(EXE readIn)
  target_link_libraries(${EXE} crosslink_lib)
  target_link_libraries(${EXE} soma_lib)
  target_link_libraries(${EXE} ${MPI_LIBRARIES})
  target_link_libraries(${EXE} ${HDF5_LIBRARIES})
  target_link_libraries(${EXE} ${HDF5_C_LIBRARIES})
  target_link_libraries(${EXE} ${soma_link_flags})

  install(TARGETS ${EXE} DESTINATION bin)
endforeach(EXE)

foreach(HEADER ${crosslink_header})
  install(FILES ${HEADER} DESTINATION include)
endforeach(HEADER ${crosslink_header})
