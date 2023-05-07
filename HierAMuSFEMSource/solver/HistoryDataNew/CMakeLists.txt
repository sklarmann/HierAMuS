# Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES}
					 ${CMAKE_CURRENT_SOURCE_DIR}/HistoryDataManager.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/HistoryDataIterator.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/HistoryDataSetup.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/HistoryDataStructure.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/HistoryDataStructureCombinded.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
					 ${CMAKE_CURRENT_SOURCE_DIR}/HistoryDataManager.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/HistoryDataIterator.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/HistoryDataSetup.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/HistoryDataStructure.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/HistoryDataStructureCombinded.h
		PARENT_SCOPE)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR} ${SpectraDir})
#
#
#if(WIN32)
#	add_library(CPPFEM_solver STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_solver SHARED ${SOURCEFILES} ${HEADERS})
#endif()
#
#if(USE_CATALYST)
#target_link_libraries(CPPFEM_solver PRIVATE ParaView::Catalyst ParaView::RemotingServerManager ParaView::RemotingAnimation ParaView::RemotingLive VTK::CommonDataModel)
#endif()