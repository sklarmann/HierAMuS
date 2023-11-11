# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES} 
					 ${CMAKE_CURRENT_SOURCE_DIR}/Plane.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS} 
					 ${CMAKE_CURRENT_SOURCE_DIR}/Plane.h
		             ${CMAKE_CURRENT_SOURCE_DIR}/geometry.h
		PARENT_SCOPE)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
#
#
#if(WIN32)
#	add_library(CPPFEM_Math STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_Math SHARED ${SOURCEFILES} ${HEADERS})
#endif()
