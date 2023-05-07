# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES}
		#${CMAKE_CURRENT_SOURCE_DIR}/GaussPoints.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/KernelShapes.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/LagrangeShape.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/helperClasses/IntegrationPointsManagement.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/IntegrationPoints.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/IntegrationPointsIterator.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/dataClasses/IntegrationPointsBase.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/dataClasses/GaussPoints1D.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/dataClasses/GaussPoints2D.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/dataClasses/GaussPoints2DTriangle.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/dataClasses/GaussPoints3D.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/LobattoShapes.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/LegendreShapes.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
		#${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/GaussPoints.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/KernelShapes.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/LagrangeShape.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/helperClasses/IntegrationPointsManagement.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/IntegrationPoints.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/IntegrationPointsIterator.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/dataClasses/IntegrationPointsBase.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/dataClasses/GaussPoints1D.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/dataClasses/GaussPoints2D.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/dataClasses/GaussPoints2DTriangle.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/IntegrationsPoints/dataClasses/GaussPoints3D.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/LobattoShapes.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/LegendreShapes.h
		PARENT_SCOPE)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR})
#
#
#if(WIN32)
#	add_library(CPPFEM_ShapeFunctions STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_ShapeFunctions SHARED ${SOURCEFILES} ${HEADERS})
#endif()