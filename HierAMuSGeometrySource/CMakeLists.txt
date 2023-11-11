# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause


set(SOURCEFILES)
set(HEADERS)
#set(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../HierAMuSFEMSource/HelperFunctions.h)

add_subdirectory(geometry)
add_subdirectory(shapefunctions)
add_subdirectory(math)


if(WIN32)
	set(HELPSOURCE ${SOURCEFILES} ${HEADERS})
	add_library(HierAMuS_Geometry SHARED ${HELPSOURCE})
else()
	add_library(HierAMuS_Geometry SHARED ${SOURCEFILES})
endif()

target_link_libraries(HierAMuS_Geometry PRIVATE spdlog::spdlog Eigen3::Eigen HierAMuS_Plot HierAMuS_Equations HierAMuS_Loads)

target_include_directories(HierAMuS_Geometry PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../HierAMuSPlotSource ${CMAKE_CURRENT_SOURCE_DIR}/../HierAMuSFEMSource)
target_include_directories(HierAMuS_Geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(HierAMuS_Geometry PROPERTIES DEBUG_POSTFIX "d")

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

if(USE_MKL)
target_link_libraries(HierAMuS_Geometry PRIVATE MKL::MKL)
endif()



install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/HierAMuS" # source directory
        DESTINATION . # target directory
        )


install(TARGETS HierAMuS_Equations
        LIBRARY DESTINATION .
        ARCHIVE DESTINATION .
        RUNTIME DESTINATION .
        )