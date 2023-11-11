# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES  ${SOURCEFILES}
				# Volumes
				${CMAKE_CURRENT_SOURCE_DIR}/VolumesData.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/VolumesRuntime.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearBrickData.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearBrickRuntime.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearPrismData.cpp
				${CMAKE_CURRENT_SOURCE_DIR}/LinearPrismRuntime.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
			# Volumes
			${CMAKE_CURRENT_SOURCE_DIR}/VolumesData.h
			${CMAKE_CURRENT_SOURCE_DIR}/VolumesDataInterface.h
			${CMAKE_CURRENT_SOURCE_DIR}/VolumesRuntime.h
			${CMAKE_CURRENT_SOURCE_DIR}/VolumesRuntimeDataInterface.h
			${CMAKE_CURRENT_SOURCE_DIR}/VolumesH1Interface.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearBrickData.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearBrickRuntime.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearPrismData.h
			${CMAKE_CURRENT_SOURCE_DIR}/LinearPrismRuntime.h
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