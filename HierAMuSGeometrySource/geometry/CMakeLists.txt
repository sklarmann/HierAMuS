# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause
add_subdirectory(Edges)
add_subdirectory(Faces)
add_subdirectory(Volumes)
add_subdirectory(Special)

set(SOURCEFILES  ${SOURCEFILES}
				# General Code
				${CMAKE_CURRENT_SOURCE_DIR}/GeometryShape.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/GeometryBaseData.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/GeometryBaseRuntime.cpp
				#${CMAKE_CURRENT_SOURCE_DIR}/GeometryTypes.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/GeometryData.cpp
				# Data Lists
				#${CMAKE_CURRENT_SOURCE_DIR}/datalists/compressedIndexToIndexList.cpp
				# Vertices
				${CMAKE_CURRENT_SOURCE_DIR}/VertexData.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/VertexRuntime.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
			# General Code
			${CMAKE_CURRENT_SOURCE_DIR}/GeometryTypes.h
			${CMAKE_CURRENT_SOURCE_DIR}/GeometryShape.h
			${CMAKE_CURRENT_SOURCE_DIR}/GeometryBaseData.h
			${CMAKE_CURRENT_SOURCE_DIR}/GeometryBaseRuntime.h
			${CMAKE_CURRENT_SOURCE_DIR}/GeometryData.h
			# Data Lists
			#${CMAKE_CURRENT_SOURCE_DIR}/datalists/compressedIndexToIndexList.h
		        ${CMAKE_CURRENT_SOURCE_DIR}/datalists/GeometryList.h
		        ${CMAKE_CURRENT_SOURCE_DIR}/datalists/GeometryListSingle.h
			# Vertices
			${CMAKE_CURRENT_SOURCE_DIR}/VertexData.h
			${CMAKE_CURRENT_SOURCE_DIR}/VertexRuntime.h
		PARENT_SCOPE)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR})
#
#
#if(WIN32)
#	add_library(CPPFEM_Geometry STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_Geometry SHARED ${SOURCEFILES} ${HEADERS})
#endif()
#
#if(USE_CATALYST)
#target_link_libraries(CPPFEM_Geometry PRIVATE ParaView::Catalyst ParaView::RemotingServerManager ParaView::RemotingAnimation ParaView::RemotingLive VTK::CommonDataModel)
#endif()