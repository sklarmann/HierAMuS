# Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES}
	${CMAKE_CURRENT_SOURCE_DIR}/HomogeniztaionBase.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/HomogenizationBeam.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Homogenization2DSolid.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Homogenization3DSolid.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/HomogenizationShell.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Homogenization3DThermoMechBeam.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
	${CMAKE_CURRENT_SOURCE_DIR}/HomogeniztaionBase.h
	${CMAKE_CURRENT_SOURCE_DIR}/HomogenizationBeam.h
	${CMAKE_CURRENT_SOURCE_DIR}/Homogenization2DSolid.h
	${CMAKE_CURRENT_SOURCE_DIR}/Homogenization3DSolid.h
	${CMAKE_CURRENT_SOURCE_DIR}/HomogenizationShell.h
	${CMAKE_CURRENT_SOURCE_DIR}/Homogenization3DThermoMechBeam.h
		PARENT_SCOPE)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR} ${SpectraDir})
#
#
#if(WIN32)
#	add_library(CPPFEM_solver STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_solver SHARED ${SOURCEFILES} ${HEADERS})
#endif()
#
#if(USE_CATALYST)
#target_link_libraries(CPPFEM_solver PRIVATE ParaView::Catalyst ParaView::RemotingServerManager ParaView::RemotingAnimation ParaView::RemotingLive VTK::CommonDataModel)
#endif()