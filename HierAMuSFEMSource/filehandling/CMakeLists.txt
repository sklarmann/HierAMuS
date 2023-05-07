# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES}
	    ${CMAKE_CURRENT_SOURCE_DIR}/filehandling.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
            ${CMAKE_CURRENT_SOURCE_DIR}/filehandling.h
					  PARENT_SCOPE)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE})
#
#
#if(WIN32)
#	add_library(CPPFEM_FileHandling STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_FileHandling SHARED ${SOURCEFILES} ${HEADERS})
#endif()