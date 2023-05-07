
#message(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})


#if(WIN32)
#    configure_file(setup.py.win ${CMAKE_SOURCE_DIR}/pybind/setup.py)
#else()
#    configure_file(setup.py.linux ${CMAKE_SOURCE_DIR}/pybind/setup.py)
#endif()

find_package(pybind11)



add_custom_target(PyPackageInstall
        COMMAND "${PYTHON_EXECUTABLE}" -m pip install -e .
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/pybind"
        COMMENT "Installing wheel file"
        )
add_dependencies(PyPackageInstall HierAMuS_FEM HierAMuS_Plot HierAMuSPyWrapper)


add_custom_command(TARGET PyPackageInstall PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:HierAMuSPyWrapper> ${CMAKE_SOURCE_DIR}/pybind/HierAMuS/${targetfile}
)
add_custom_command(TARGET PyPackageInstall PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:HierAMuS_FEM> ${CMAKE_SOURCE_DIR}/pybind/HierAMuS/${targetfile}
)

