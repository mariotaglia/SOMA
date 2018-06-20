if(DOXYGEN_FOUND)
  if(ENABLE_DOXYGEN)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

    add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen.stamp
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      COMMAND doxygen ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      COMMAND cmake -E touch ${CMAKE_CURRENT_BINARY_DIR}/doxygen.stamp
      COMMENT "Generating the Doxygen documentation"
      VERBATIM )
    add_custom_target( doc ALL
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen.stamp )
    
    option(INSTALL_DOC "Install the doxygen html documentation." OFF)
    
    if(INSTALL_DOC)
      
      install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ DESTINATION doc/ )
    endif(INSTALL_DOC)
endif(ENABLE_DOXYGEN)
endif(DOXYGEN_FOUND)
