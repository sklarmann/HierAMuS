# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES  ${SOURCEFILES}
				# General Code
				${CMAKE_CURRENT_SOURCE_DIR}/Base.cpp
				#${CMAKE_CURRENT_SOURCE_DIR}/GeometryTypes.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/GeometryData.cpp
				# Data Lists
				#${CMAKE_CURRENT_SOURCE_DIR}/datalists/compressedIndexToIndexList.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/datalists/vertexlist.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/datalists/edgelist.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/datalists/facelist.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/datalists/volumelist.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/datalists/speciallist.cpp
				# Vertices
				${CMAKE_CURRENT_SOURCE_DIR}/Vertex.cpp
				# Edges
				${CMAKE_CURRENT_SOURCE_DIR}/Edges.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearEdge.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/QuadraticEdge.cpp
				# Faces
				${CMAKE_CURRENT_SOURCE_DIR}/Faces.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearTriangle.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearQuadrilateral.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/QuadraticQuadrilateral.cpp
				# Volumes
				${CMAKE_CURRENT_SOURCE_DIR}/Volumes.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearBrick.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearPrism.cpp
				# Special
				${CMAKE_CURRENT_SOURCE_DIR}/Special.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/BeamInterface2D.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/BeamInterface3D.cpp
				# Scaled Boundary elements
				${CMAKE_CURRENT_SOURCE_DIR}/ScaledBoundary2D.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
			# General Code
			${CMAKE_CURRENT_SOURCE_DIR}/GeometryTypes.h
			${CMAKE_CURRENT_SOURCE_DIR}/Base.h
			${CMAKE_CURRENT_SOURCE_DIR}/GeometryData.h
			# Data Lists
			#${CMAKE_CURRENT_SOURCE_DIR}/datalists/compressedIndexToIndexList.h
			${CMAKE_CURRENT_SOURCE_DIR}/datalists/vertexlist.h
			${CMAKE_CURRENT_SOURCE_DIR}/datalists/edgelist.h
			${CMAKE_CURRENT_SOURCE_DIR}/datalists/facelist.h
			${CMAKE_CURRENT_SOURCE_DIR}/datalists/volumelist.h
			${CMAKE_CURRENT_SOURCE_DIR}/datalists/speciallist.h
			${CMAKE_CURRENT_SOURCE_DIR}/datalists/datalistiterator.h
			# Vertices
			${CMAKE_CURRENT_SOURCE_DIR}/Vertex.h
			# Edges
			${CMAKE_CURRENT_SOURCE_DIR}/Edges.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearEdge.h
			${CMAKE_CURRENT_SOURCE_DIR}/QuadraticEdge.h
			# Faces
			${CMAKE_CURRENT_SOURCE_DIR}/Faces.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearTriangle.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearQuadrilateral.h
			${CMAKE_CURRENT_SOURCE_DIR}/QuadraticQuadrilateral.h
			# Volumes
			${CMAKE_CURRENT_SOURCE_DIR}/Volumes.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearBrick.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearPrism.h
			# Special
			${CMAKE_CURRENT_SOURCE_DIR}/Special.h
			${CMAKE_CURRENT_SOURCE_DIR}/BeamInterface2D.h
			${CMAKE_CURRENT_SOURCE_DIR}/BeamInterface3D.h
			# Scaled Boundary elements
			${CMAKE_CURRENT_SOURCE_DIR}/ScaledBoundary2D.h
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