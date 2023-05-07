# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES}
		${CMAKE_CURRENT_SOURCE_DIR}/vtkPlotClass.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/vtkPlotClassBase.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
		${CMAKE_CURRENT_SOURCE_DIR}/vtkplotClassBase.h
			${CMAKE_CURRENT_SOURCE_DIR}/vtkplotClass.h
		PARENT_SCOPE)


if(USE_CATALYST)


set(HEADERS ${HEADERS} 
			${CMAKE_CURRENT_SOURCE_DIR}/catalystcxxinterface/gridHandler.h
			${CMAKE_CURRENT_SOURCE_DIR}/catalystcxxinterface/management.h
			${CMAKE_CURRENT_SOURCE_DIR}/paraviewManagementClass.h
		PARENT_SCOPE)
endif()

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR})
#
#
#if(WIN32)
#	add_library(HierAMuS_Plot STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(HierAMuS_Plot SHARED ${SOURCEFILES} ${HEADERS})
#endif()
#
#if(USE_CATALYST)
#target_link_libraries(HierAMuS_Plot PRIVATE ParaView::Catalyst ParaView::RemotingServerManager ParaView::RemotingAnimation ParaView::RemotingLive VTK::CommonDataModel)
#endif()