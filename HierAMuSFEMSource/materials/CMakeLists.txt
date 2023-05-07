# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

add_subdirectory(2D)
add_subdirectory(3D)
add_subdirectory(Homogenization)

set(SOURCEFILES ${SOURCEFILES}
		${CMAKE_CURRENT_SOURCE_DIR}/Material.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/MaterialList.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/MaterialformulationList.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/ElementformulationList.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/GenericMaterialFormulation.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/MaterialFormulation2.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
		${CMAKE_CURRENT_SOURCE_DIR}/Material.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/MaterialList.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/MaterialformulationList.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/ElementformulationList.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/GenericMaterialFormulation.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/MaterialFormulation2.h
		PARENT_SCOPE)


#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR})
#
#if(WIN32)
#	add_library(CPPFEM_Materials STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_Materials SHARED ${SOURCEFILES} ${HEADERS})
#endif()
#
#if(USE_CATALYST)
#target_link_libraries(CPPFEM_Materials PRIVATE ParaView::Catalyst ParaView::RemotingServerManager ParaView::RemotingAnimation ParaView::RemotingLive VTK::CommonDataModel)
#endif()