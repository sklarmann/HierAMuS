# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES  ${SOURCEFILES}
				# Edges
				${CMAKE_CURRENT_SOURCE_DIR}/EdgesData.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/EdgesRuntime.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearEdgeData.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearEdgeRuntime.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/QuadraticEdge.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
			# Edges
			${CMAKE_CURRENT_SOURCE_DIR}/EdgesData.h
			${CMAKE_CURRENT_SOURCE_DIR}/EdgesDataInterface.h
			${CMAKE_CURRENT_SOURCE_DIR}/EdgesRuntime.h
			${CMAKE_CURRENT_SOURCE_DIR}/EdgesRuntimeDataInterface.h
			${CMAKE_CURRENT_SOURCE_DIR}/EdgeH1ShapesInterface.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearEdgeData.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearEdgeRuntime.h
			${CMAKE_CURRENT_SOURCE_DIR}/QuadraticEdge.h
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