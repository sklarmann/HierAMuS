# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES  ${SOURCEFILES}
				# Faces
				${CMAKE_CURRENT_SOURCE_DIR}/FacesData.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/FacesRuntime.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearTriangleData.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearTriangleRuntime.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearQuadrilateralData.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearQuadrilateralRuntime.cpp
				# Scaled Boundary elements
				${CMAKE_CURRENT_SOURCE_DIR}/ScaledBoundary2DData.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/ScaledBoundary2DRuntime.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
			# Faces
			${CMAKE_CURRENT_SOURCE_DIR}/FaceOrientationFlags.h
			${CMAKE_CURRENT_SOURCE_DIR}/FacesData.h
			${CMAKE_CURRENT_SOURCE_DIR}/FacesDataInterface.h
			${CMAKE_CURRENT_SOURCE_DIR}/FacesRuntime.h
			${CMAKE_CURRENT_SOURCE_DIR}/FacesRuntimeDataInterface.h
			${CMAKE_CURRENT_SOURCE_DIR}/FacesH1Interface.h
			${CMAKE_CURRENT_SOURCE_DIR}/FacesHDivInterface.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearTriangleData.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearTriangleRuntime.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearQuadrilateralData.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearQuadrilateralRuntime.h
			# Scaled Boundary elements
			${CMAKE_CURRENT_SOURCE_DIR}/ScaledBoundary2DData.h
			${CMAKE_CURRENT_SOURCE_DIR}/ScaledBoundary2DRuntime.h
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