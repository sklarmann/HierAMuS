############################################################################################
#
# Eigen check
#
############################################################################################

find_package(Eigen3 QUIET
HINTS "${STAGED_INSTALL_PREFIX}/share/eigen3/cmake" "${STAGED_INSTALL_PREFIX}/CMake"
)

if(NOT EIGEN3_FOUND)
    message(STATUS "Eigen3 not found. Downloading and installing it in the build directory.")
    configure_file(external/Eigen3/CMakeLists.txt.in eigen-download/CMakeLists.txt)


    execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/eigen-download"
    )
    execute_process(COMMAND "${CMAKE_COMMAND}" --build . --config Release
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/eigen-download"
    )

    if(WIN32 AND NOT CYGWIN)
      set(DEF_Eigen3_DIR ${STAGED_INSTALL_PREFIX}/CMake)
    else()
      set(DEF_Eigen3_DIR ${STAGED_INSTALL_PREFIX}/share/eigen3/cmake)
    endif()
      #file(TO_NATIVE_PATH "${DEF_Eigen3_DIR}" DEF_Eigen3_DIR)
    set(Eigen3_DIR ${DEF_Eigen3_DIR}
        CACHE PATH "Path to internally built Eigen3Config.cmake" FORCE)

endif()

############################################################################################
#
# Spectra check
#
############################################################################################
macro(find_dependency args)
        message(${args})
        find_package(${args})
endmacro()
find_package(Spectra QUIET
HINTS "${STAGED_INSTALL_PREFIX}/cmake/")

if(NOT Spectra_FOUND)
    message(STATUS "Spectra not found. Downloading and installing it in the build directory.")
    configure_file(external/Spectra/CMakeLists.txt.in spectra-download/CMakeLists.txt)


    execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/spectra-download"
    )
    execute_process(COMMAND "${CMAKE_COMMAND}" --build . --config Release
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/spectra-download"
    )

    if(WIN32 AND NOT CYGWIN)
      set(DEF_Spectra_DIR ${STAGED_INSTALL_PREFIX}/CMake)
    else()
      set(DEF_Spectra_DIR ${STAGED_INSTALL_PREFIX}/cmake)
    endif()




    #file(TO_NATIVE_PATH "${DEF_Spectra_DIR}" DEF_Spectra_DIR)
    set(spectra_DIR ${DEF_Spectra_DIR}
        CACHE PATH "Path to internally built spectra-config.cmake" FORCE)


endif()

############################################################################################
#
# spdlogger check
#
############################################################################################
find_package(spdlog QUIET
HINTS "${STAGED_INSTALL_PREFIX}/lib/cmake/")

if(NOT spdlog_FOUND)
  message(STATUS "Spdlogger not found. Downloading and installing it in the build directory.")
  configure_file(external/spdlogger/CMakeLists.txt.in spdlogger-download/CMakeLists.txt)


  execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/spdlogger-download"
  )
  execute_process(COMMAND "${CMAKE_COMMAND}" --build . --config Release
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/spdlogger-download"
  )

  if(WIN32 AND NOT CYGWIN)
    set(DEF_Spdlogger_DIR ${STAGED_INSTALL_PREFIX}/CMake)
  else()
    set(DEF_Spdlogger_DIR ${STAGED_INSTALL_PREFIX}/lib/cmake/spdlog)
  endif()
    #file(TO_NATIVE_PATH "${DEF_Eigen3_DIR}" DEF_Eigen3_DIR)
  set(spdlog_DIR ${DEF_Spdlogger_DIR}
      CACHE PATH "Path to internally built Eigen3Config.cmake" FORCE)

      message(${spdlog_DIR})
endif()