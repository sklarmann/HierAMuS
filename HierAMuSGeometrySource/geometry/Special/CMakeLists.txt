# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES  ${SOURCEFILES}
				# Special
				${CMAKE_CURRENT_SOURCE_DIR}/Special.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/BeamInterface2D.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/BeamInterface3D.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
			# Special
			${CMAKE_CURRENT_SOURCE_DIR}/Special.h
			${CMAKE_CURRENT_SOURCE_DIR}/BeamInterface2D.h
			${CMAKE_CURRENT_SOURCE_DIR}/BeamInterface3D.h
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