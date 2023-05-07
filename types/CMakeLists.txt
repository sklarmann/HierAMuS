# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(CDIR ${CMAKE_CURRENT_SOURCE_DIR}/types)

set(HEADERS 
			 ${CMAKE_CURRENT_SOURCE_DIR}/MatrixTypes.h
					 )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR})

if(WIN32)
	add_library(CPPFEM_Types STATIC ${SOURCEFILES} ${HEADERS})
else()
	add_library(CPPFEM_Types SHARED ${SOURCEFILES} ${HEADERS})
endif()

target_link_libraries(CPPFEM_Types)