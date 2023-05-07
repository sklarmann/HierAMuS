# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES}
	         ${CMAKE_CURRENT_SOURCE_DIR}/plotCommands.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
					 ${CMAKE_CURRENT_SOURCE_DIR}/plotCommands.h
		PARENT_SCOPE)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE})
#
#if(WIN32)
#	add_library(CPPFEM_Control STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_Control SHARED ${SOURCEFILES} ${HEADERS})
#endif()
#
#
#if(USE_CATALYST)
#target_link_libraries(CPPFEM_Control PRIVATE ParaView::Catalyst VTK::CommonDataModel)
#endif()