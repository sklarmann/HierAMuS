# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(CDIR ${CMAKE_CURRENT_SOURCE_DIR}/equations)


set(SOURCEFILES ${SOURCEFILES}
	         ${CMAKE_CURRENT_SOURCE_DIR}/EquationHandler.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/GenericNodes.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/NodeSet.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/Nodetypes.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/DegreeOfFreedom.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
	          ${CMAKE_CURRENT_SOURCE_DIR}/DofStatus.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/EquationHandler.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/GenericNodes.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/NodeSet.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/Nodetypes.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/DegreeOfFreedom.h
		PARENT_SCOPE)



#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR})
#
#
#if(WIN32)
#	add_library(CPPFEM_Equations STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_Equations SHARED ${SOURCEFILES} ${HEADERS})
#endif()
#
#if(USE_CATALYST)
#target_link_libraries(CPPFEM_Equations PUBLIC ParaView::Catalyst VTK::CommonDataModel)
#endif()