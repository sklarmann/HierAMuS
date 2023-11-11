# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES)
set(HEADRS ${CMAKE_CURRENT_SOURCE_DIR}/HelperFunctions.h)

add_subdirectory(control)
add_subdirectory(elementFormulations)
add_subdirectory(finiteElements)
add_subdirectory(materials)
add_subdirectory(math)
add_subdirectory(plot)
add_subdirectory(pointercollection)
add_subdirectory(solver)
#add_subdirectory(stringoperations)

#
#
if(WIN32)
	set(HELPSOURCE ${SOURCEFILES} ${HEADERS})
	add_library(HierAMuS_FEM SHARED ${HELPSOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/../eigen.natvis)
else()
	add_library(HierAMuS_FEM SHARED ${SOURCEFILES} ${CMAKE_CURRENT_SOURCE_DIR}/../eigen.natvis)
endif()
#

target_link_libraries(HierAMuS_FEM PUBLIC Eigen3::Eigen HierAMuS_Plot HierAMuS_Equations HierAMuS_Loads HierAMuS_Geometry)

if(USE_OPENMP)
  target_link_libraries(HierAMuS_FEM PUBLIC OpenMP::OpenMP_CXX)
endif()

set_target_properties(HierAMuS_FEM PROPERTIES DEBUG_POSTFIX "d")

if(WIN32)

else()
    target_link_options(HierAMuS_FEM PRIVATE "-Wl,--exclude-libs,ALL")
endif()
target_compile_definitions(HierAMuS_FEM PRIVATE BUILDCPPFEMFEM)


if(WIN32)
foreach(FILE ${HELPSOURCE})
    # Get the directory of the source file
    get_filename_component(PARENT_DIR "${FILE}" DIRECTORY)

    # Remove common directory prefix to make the group
    string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}" "" GROUP "${PARENT_DIR}")

    # Make sure we are using windows slashes
    string(REPLACE "/" "\\" GROUP "${GROUP}")

    # Group into "Source Files" and "Header Files"
    if ("${FILE}" MATCHES ".*\\.cpp")
       set(GROUP "Source Files${GROUP}")
    elseif("${FILE}" MATCHES ".*\\.h")
       set(GROUP "Header Files${GROUP}")
    endif()

    source_group("${GROUP}" FILES "${FILE}")
endforeach()
endif()


target_link_libraries(HierAMuS_FEM PRIVATE Spectra::Spectra spdlog::spdlog)
target_include_directories(HierAMuS_FEM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_SOURCE_DIR}/../HierAMuSPlotSource)

if(USE_MKL)
target_link_libraries(HierAMuS_FEM PUBLIC MKL::MKL)
endif()
set(CPPYHEADERS ${HEADERS} PARENT_SCOPE)
set(SOURCEFILES ${SOURCEFILES} PARENT_SCOPE)


# install(TARGETS HierAMuS_FEM EXPORT HierAMuS_FEMTargets
#   LIBRARY DESTINATION .
#   ARCHIVE DESTINATION .
#   RUNTIME DESTINATION .
#   INCLUDES DESTINATION HierAMuS_FEM/include
# )

# install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/" # source directory
#         DESTINATION "HierAMuS_FEM/include" # target directory
#         FILES_MATCHING # install only matched files
#         PATTERN "*.h" # select header files
#         )
# install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../types" # source directory
#         DESTINATION "HierAMuS_FEM/include/" # target directory
#         FILES_MATCHING # install only matched files
#         PATTERN "*.h" # select header files
#         )

# install(
#   FILES
#         ${CMAKE_CURRENT_SOURCE_DIR}/../datatypes.h

# #    #"${CMAKE_CURRENT_BINARY_DIR}/testlib_export.h"
#   DESTINATION
#     HierAMuS_FEM/include
# #  COMPONENT
# #    Devel
# )


# include(CMakePackageConfigHelpers)
# write_basic_package_version_file(
#   "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_FEM/HierAMuS_FEMConfigVersion.cmake"
#   VERSION ${VERSION}
#   COMPATIBILITY AnyNewerVersion
# )

# export(EXPORT HierAMuS_FEMTargets
#   FILE "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_FEM/HierAMuS_FEMTargets.cmake"
#   NAMESPACE CPPFEM::
# )
# configure_file(cmake/HierAMuS_FEMConfig.cmake
#   "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_FEM/HierAMuS_FEMConfig.cmake"
#   COPYONLY
# )

# set(ConfigPackageLocation lib/cmake/HierAMuS_FEM)
# install(EXPORT HierAMuS_FEMTargets
#   FILE
#     HierAMuS_FEMTargets.cmake
#   NAMESPACE
#     CPPFEM::
#   DESTINATION
#     ${ConfigPackageLocation}
# )
# install(
#   FILES
#     cmake/HierAMuS_FEMConfig.cmake
#     "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_FEM/HierAMuS_FEMConfigVersion.cmake"
#   DESTINATION
#     ${ConfigPackageLocation}
#   COMPONENT
#     Devel
# )




