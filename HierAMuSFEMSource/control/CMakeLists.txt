# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(CDIR ${CMAKE_CURRENT_SOURCE_DIR}/control)


set(SOURCEFILES ${SOURCEFILES}					 
		${CMAKE_CURRENT_SOURCE_DIR}/HandlingStructs.cpp
	         #${CMAKE_CURRENT_SOURCE_DIR}/stringCommandHandler.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/OutputHandler.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/FEMExceptions.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/ParameterList.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
					 ${CMAKE_CURRENT_SOURCE_DIR}/HandlingStructs.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/OutputHandler.h
					 #${CMAKE_CURRENT_SOURCE_DIR}/stringCommandHandler.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/FEMExceptions.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/BinaryWrite.h
		${CMAKE_CURRENT_SOURCE_DIR}/ParameterList.h
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