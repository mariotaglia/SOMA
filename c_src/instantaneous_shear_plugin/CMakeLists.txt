include_directories( "${SOMA_SOURCE_DIR}/c_src/" )

#Modify the name of this list of source files to your plugin name
set(instantaneous_shear_plugin_source
  #add all .c files if any
  shear_lib.c
  )
#Modify the name of this list of header files to your plugin name
set(instantaneous_shear_plugin_header
  #add all .h files if any
  shear_lib.h
  )

#Use your source list to create a new library with your plugin name.
add_library(instantaneous_shear_plugin_lib ${instantaneous_shear_plugin_source})

#Add the executable of your plugin with its source file containing the "main()" function
add_executable(INSTANTANEOUS_SHEAR shear.c)

#Modify the names according to your plugin
foreach(EXE INSTANTANEOUS_SHEAR)
  target_link_libraries(${EXE} instantaneous_shear_plugin_lib)
  target_link_libraries(${EXE} soma_lib)
  target_link_libraries(${EXE} ${MPI_LIBRARIES})
  target_link_libraries(${EXE} ${HDF5_LIBRARIES})
  target_link_libraries(${EXE} ${HDF5_C_LIBRARIES})
  target_link_libraries(${EXE} ${soma_link_flags})

  install(TARGETS ${EXE} DESTINATION bin)
endforeach(EXE)

#Modify the names according to your plugin
foreach(HEADER ${instantaneous_shear_plugin_header})
  install(FILES ${HEADER} DESTINATION include)
endforeach(HEADER ${instantaneous_shear_plugin_header})
