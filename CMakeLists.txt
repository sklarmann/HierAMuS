# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.16)
project (HierAMuS
VERSION 0.2 )


############################################################################################
#
# Basic CXX settings
#
############################################################################################
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON) 

############################################################################################
#
# Necessary variables for automatic dependencies installation
#
############################################################################################
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type set to ${CMAKE_BUILD_TYPE}")

set(STAGED_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/stage)
message(STATUS "${PROJECT_NAME} staged install: ${STAGED_INSTALL_PREFIX}")


include(GenerateExportHeader)
include(ExternalProject)

############################################################################################
#
# Environment check
#
############################################################################################
#include(${CMAKE_SOURCE_DIR}/cmake/CheckEnvironment.cmake)



############################################################################################
#
# MKL Settings
#
############################################################################################
include(${CMAKE_SOURCE_DIR}/cmake/MKLSettings.cmake)



############################################################################################
#
# OpenMP Settings
#
############################################################################################
set(USE_OPENMP OFF CACHE BOOL "Enable OpenMP")
if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    add_definitions(/DOPENMP)
endif()






############################################################################################
#
# Precision Settings
#
############################################################################################
include(${CMAKE_SOURCE_DIR}/cmake/PrecisionSettings.cmake)

############################################################################################
#
# Boost library
#
############################################################################################
if(BOOST_REQUIRED)
	find_package(Boost REQUIRED)
	include_directories(${Boost_INCLUDE_DIRS})
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/types)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

############################################################################################
#
# Eigen3 and Spectra Library
#
############################################################################################
include(external/installLibs.cmake)
find_package(Eigen3 REQUIRED
HINTS "${STAGED_INSTALL_PREFIX}/share/eigen3/cmake" "${STAGED_INSTALL_PREFIX}/CMake")
find_package(spectra REQUIRED
HINTS "${STAGED_INSTALL_PREFIX}/cmake/")
find_package(spdlog REQUIRED
HINTS "${STAGED_INSTALL_PREFIX}/lib/cmake/spdlog")


set(USE_SPECTRA ON)
set(USE_CATALYST ON)
add_definitions(/DUSE_SPECTRA)
add_definitions(/DUSE_VTK)


set(PREPOCESSOR_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR})


if(WIN32)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        #SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /std:c++17")
        
    else()
	#set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
      add_compile_options(/std:c++17 /MP /bigobj /Ob2)
      add_compile_options("$<$<CONFIG:DEBUG>:/MD>")
	SET(CMAKE_CXX_LINK_FLAGS   "${CMAKE_CXX_LINK_FLAGS} /INCREMENTAL")
    endif()
else()
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    SET(CMAKE_INSTALL_RPATH "\${ORIGIN}")
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wunused -Wuninitialized")
    
    #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined,memory -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer -fsanitize-blacklist=${CMAKE_CURRENT_SOURCE_DIR}/blacklist.txt -g -O0")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
    #SET(CMAKE_CXX_LINK_FLAGS  "${CMAKE_CXX_LINK_FLAGS} -fsanitize=undefined,memory")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        add_definitions(/DGNU)
        SET(CMAKE_CXX_LINK_FLAGS  "${CMAKE_CXX_LINK_FLAGS} -lquadmath")
    endif()
endif()


if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else()
    add_compile_options(-fPIC)
endif()



#
#  Adding the FEMSource Folder to the project
#
set(BUILD_FEMSOURCE ON CACHE BOOL "Link the simulator with Catalyst and enable postprocessing")
if(BUILD_FEMSOURCE)
    add_subdirectory(HierAMuSEquationsSource)
    add_subdirectory(HierAMuSLoadSource)
    add_subdirectory(HierAMuSPlotSource)
    add_subdirectory(HierAMuSGeometrySource)
    add_subdirectory(HierAMuSFEMSource)
endif()

#
#  Adding the Python Wrapper to the project
#
set(BUILD_PYBINDWRAPPER ON CACHE BOOL "Link the simulator with Catalyst and enable postprocessing")
if(BUILD_PYBINDWRAPPER)
    add_subdirectory(pybind)
    add_subdirectory(FEMPythonPackage)
endif()





# Get all propreties that cmake supports
if(NOT CMAKE_PROPERTY_LIST)
    execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)

    # Convert command output into a CMake list
    string(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
    string(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
endif()

function(print_properties)
    message("CMAKE_PROPERTY_LIST = ${CMAKE_PROPERTY_LIST}")
endfunction()

function(print_target_properties target)
    if(NOT TARGET ${target})
      message(STATUS "There is no target named '${target}'")
      return()
    endif()

    foreach(property ${CMAKE_PROPERTY_LIST})
        string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" property ${property})

        # Fix https://stackoverflow.com/questions/32197663/how-can-i-remove-the-the-location-property-may-not-be-read-from-target-error-i
        if(property STREQUAL "LOCATION" OR property MATCHES "^LOCATION_" OR property MATCHES "_LOCATION$")
            continue()
        endif()

        get_property(was_set TARGET ${target} PROPERTY ${property} SET)
        if(was_set)
            get_target_property(value ${target} ${property})
            message("${target} ${property} = ${value}")
        endif()
    endforeach()
endfunction()
#print_target_properties(HierAMuS_FEM)
#print_target_properties(Eigen3::Eigen)
#print_target_properties(CPPFEMPythonCppyy)


