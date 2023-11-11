# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause




set(LOCALSOURCEFILES
	                 ${CMAKE_CURRENT_SOURCE_DIR}/LoadList.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/PropfunctionHandler.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/Timefunction.cpp
		)

set(LOCALHEADERS
	                  ${CMAKE_CURRENT_SOURCE_DIR}/LoadList.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/PropfunctionHandler.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/Timefunction.h
		)


if(WIN32)
	set(HELPSOURCE ${LOCALSOURCEFILES} ${LOCALHEADERS})
	add_library(HierAMuS_Loads SHARED ${HELPSOURCE})
else()
	add_library(HierAMuS_Loads SHARED ${LOCALSOURCEFILES})
endif()

target_link_libraries(HierAMuS_Loads PRIVATE spdlog::spdlog)

target_include_directories(HierAMuS_Loads PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(HierAMuS_Loads PROPERTIES DEBUG_POSTFIX "d")

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






# install(TARGETS HierAMuS_Loads EXPORT HierAMuS_LoadsTargets
#         LIBRARY DESTINATION .
#         ARCHIVE DESTINATION .
#         RUNTIME DESTINATION .
#         INCLUDES DESTINATION HierAMuS_Loads/include
#         )

# install( DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/" # source directory
#          DESTINATION "HierAMuS_Loads/include" # target directory
#          FILES_MATCHING # install only matched files
#          PATTERN "*.h" # select header files
# ) 

# include(CMakePackageConfigHelpers)
# write_basic_package_version_file(
#   "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_Loads/HierAMuS_LoadsConfig.cmake"
#   VERSION ${VERSION}
#   COMPATIBILITY AnyNewerVersion
# )


# export(EXPORT HierAMuS_LoadsTargets
#   FILE "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_Loads/HierAMuS_LoadsTargets.cmake"
#   NAMESPACE CPPFEM::
# )
# configure_file(cmake/HierAMuS_LoadsConfig.cmake
#   "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_Loads/HierAMuS_LoadsConfig.cmake"
#   COPYONLY
# )

# set(ConfigPackageLocation lib/cmake/HierAMuS_Loads)
# install(EXPORT HierAMuS_LoadsTargets
#   FILE
#     HierAMuS_Loadsargets.cmake
#   NAMESPACE
#     CPPFEM::
#   DESTINATION
#     ${ConfigPackageLocation}
# )
# install(
#   FILES
#     cmake/HierAMuS_LoadsConfig.cmake
#     "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_Loads/HierAMuS_LoadsConfigVersion.cmake"
#   DESTINATION
#     ${ConfigPackageLocation}
#   COMPONENT
#     Devel
# )