set(checkEnvironment $ENV{CPPFEM})

if(NOT checkEnvironment)
    message(FATAL_ERROR "
    Environment variable CPPFEM not found!!
    Please make sure to set the environment variable CPPFEM to the path where CPPFEM should be installed!
    Most commonly the path C:\\CPPFEM\\ will be used.
    Additionally add the CPPFEM variable to your PATH variable.

    Afterwards restart CMake, then the new environment variable should be recognized!
  ")
endif()

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "$ENV{CPPFEM}"
            CACHE PATH "default install path" FORCE)
endif()

# file(GLOB addFiles "ConfigData/*")
# message(${addFiles})
# foreach(CFIL ${addFiles})
#     get_filename_component(dummy ${CFIL} NAME)
#     message(${CFIL})
#         file(COPY ${CFIL}
#                 DESTINATION $ENV{CPPFEM})
#         # message("copy ${CFIL} to $ENV{FCFG}${dummy}")
#     # string(REGEX MATCH "(?!.*\\/).*" dummy ${CFIL})
#     # message(${dummy})
# endforeach()



macro(print_all_variables)
    message(STATUS "print_all_variables------------------------------------------{")
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
    message(STATUS "print_all_variables------------------------------------------}")
endmacro()




if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
	    "Chose the type of build, options are: Debug, Release, RelWithDebInfo, MinSizeRel."
    FORCE)
endif()