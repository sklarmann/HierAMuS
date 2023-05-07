############################################################################################
#
# Paraview check
#
############################################################################################


set(CATALYST_USE_MPI OFF)
find_package(ParaView QUIET
HINTS "${STAGED_INSTALL_PREFIX}/lib/cmake/paraview-5.10/")

if(NOT ParaView_FOUND)
    message(STATUS "Paraview not found. Downloading and installing it in the build directory.")
    configure_file(${CMAKE_SOURCE_DIR}/external/ParaviewCatalyst/CMakeLists.txt.in ${CMAKE_BINARY_DIR}/paraview-download/CMakeLists.txt)


    execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/paraview-download"
    )
    execute_process(COMMAND "${CMAKE_COMMAND}" --build . --config Release
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/paraview-download"
    )

    if(WIN32 AND NOT CYGWIN)
      set(DEF_Paraview_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/paraview-5.10)
    else()
      set(DEF_Paraview_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/paraview-5.10)
    endif()




    #file(TO_NATIVE_PATH "${DEF_Spectra_DIR}" DEF_Spectra_DIR)
    set(ParaView_DIR ${DEF_Paraview_DIR}
        CACHE PATH "Path to internally built paraview-config.cmake" FORCE)


endif()