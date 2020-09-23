set(soma_python_scripts
  AnaGen.py
  ConfGen.py
  create_xdmf.py
  handleAnaH5.py )

if( SINGLE_PRECISION )
  set(NUMPY_SOMA_SCALAR_TYPE "numpy.float32")
else()
  set(NUMPY_SOMA_SCALAR_TYPE "numpy.float64")
endif()

configure_file(soma_type.py.in soma_type.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/soma_type.py DESTINATION python-script)
set(SOMA_PYTHON_INSTALL_DIR "\"${CMAKE_INSTALL_PREFIX}/python-script\"")

foreach( SCRIPT ${soma_python_scripts})
  configure_file("${SCRIPT}.in" "${SCRIPT}")
  #Configure the files only as executable if they have a SHEBANG
  if(DEFINED SOMA_PYTHON_SHEBANG)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}" DESTINATION python-script PERMISSIONS OWNER_EXECUTE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_READ WORLD_EXECUTE)
  else()
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}" DESTINATION python-script PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
    endif()
endforeach(SCRIPT ${soma_python_scripts})
