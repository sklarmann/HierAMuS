# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES}
	         ${CMAKE_CURRENT_SOURCE_DIR}/pointercollection.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
	          ${CMAKE_CURRENT_SOURCE_DIR}/pointercollection.h
		PARENT_SCOPE)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR})
#
#if(WIN32)
#	add_library(CPPFEM_PointerCollection STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_PointerCollection SHARED ${SOURCEFILES} ${HEADERS})
#endif()
#
#if(USE_CATALYST)
#target_link_libraries(CPPFEM_PointerCollection PRIVATE ParaView::Catalyst ParaView::RemotingServerManager ParaView::RemotingAnimation ParaView::RemotingLive VTK::CommonDataModel)
#endif()