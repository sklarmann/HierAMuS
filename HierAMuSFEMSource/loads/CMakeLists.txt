# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES}
		${CMAKE_CURRENT_SOURCE_DIR}/PropfunctionHandler.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/Timefunction.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/GenericLoad.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/DofLoad.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/LoadList.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
		${CMAKE_CURRENT_SOURCE_DIR}/PropfunctionHandler.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/Timefunction.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/GenericLoad.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/DofLoad.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/LoadList.h
					 PARENT_SCOPE)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR})
#
#
#if(WIN32)
#	add_library(CPPFEM_Loads STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_Loads SHARED ${SOURCEFILES} ${HEADERS})
#endif()
#
#if(USE_CATALYST)
#target_link_libraries(CPPFEM_Loads PUBLIC ParaView::Catalyst VTK::CommonDataModel)
#endif()