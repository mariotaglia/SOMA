include_directories("${SOMA_SOURCE_DIR}/c_src/")

set(mpi_unity_source
        unity.c
        unity_capture.c
        unity_capture_test.c
        unity_fixture.c
        unity_mpi_addon.c
        unity_memory.c
        mesh_test.c
        ana_test.c
        )

set(mpi_unity_header
        unity.h
        unit_test_util.h
        unity_capture.h
        unity_capture_internals.h
        unity_fixture.h
        unity_memory.h
        unity_mpi_addon.h
        run_unit_tests.h
        )

add_library(mpi_unity_lib
        ${mpi_unity_source}
        )

add_executable(UNITTESTS
        unit_test_main.c
        ${CMAKE_CURRENT_BINARY_DIR}/run_unit_tests.c
        )

# extract_test_cases needs to run every time the test cases change, that is, any of the _test.c files.
# I do not know how to detect this, so I run the script on *every* invocation of make.
# add_custom_command without a DEPENDS-argument runs the script whenever one of the OUTPUTs do not exist, so if we specify a file as OUTPUT that does not exist and will not be created, we have the behaviour of running this every time.
set(EXTRACT_SCRIPT  "${CMAKE_CURRENT_SOURCE_DIR}/extract_test_cases.py")
add_custom_command(OUTPUT run_unit_tests.c non_existing_dummy.txt PRE_BUILD COMMAND ${PYTHON_EXECUTABLE} "${EXTRACT_SCRIPT}" COMMENT "locating unit tests...")

foreach(EXE UNITTESTS)
  target_link_libraries(${EXE} mpi_unity_lib)
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

endforeach(EXE)

foreach(HEADER ${mpi_unity_header})
  install(FILES ${HEADER} DESTINATION include)
endforeach()