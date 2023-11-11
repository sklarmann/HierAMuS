# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause




set(LOCALSOURCEFILES
	                 ${CMAKE_CURRENT_SOURCE_DIR}/EquationHandler.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/GenericNodes.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/NodeSet.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/NodeSetList.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/NodeSetNodeList.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/MeshIdNodeList.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/Nodetypes.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/DegreeOfFreedom.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/NodeSetManager.cpp
		)

set(LOCALHEADERS
	                  ${CMAKE_CURRENT_SOURCE_DIR}/DofStatus.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/EquationHandler.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/GenericNodes.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/NodeSet.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/NodeSetList.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/NodeSetNodeList.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/MeshIdNodeList.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/Nodetypes.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/DegreeOfFreedom.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/NodeSetManager.h
		)


if(WIN32)
	set(HELPSOURCE ${LOCALSOURCEFILES} ${LOCALHEADERS})
	add_library(HierAMuS_Equations SHARED ${HELPSOURCE})
else()
	add_library(HierAMuS_Equations SHARED ${LOCALSOURCEFILES})
endif()

target_link_libraries(HierAMuS_Equations PRIVATE spdlog::spdlog)

target_include_directories(HierAMuS_Equations PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(HierAMuS_Equations PROPERTIES DEBUG_POSTFIX "d")

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



install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/HierAMuS" # source directory
        DESTINATION . # target directory
        )


install(TARGETS HierAMuS_Equations
        LIBRARY DESTINATION .
        ARCHIVE DESTINATION .
        RUNTIME DESTINATION .
        )