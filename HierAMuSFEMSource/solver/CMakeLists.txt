# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

add_subdirectory(HistoryDataNew)
add_subdirectory(Constraints)
add_subdirectory(Homogenization)

set(SOURCEFILES ${SOURCEFILES}
	         ${CMAKE_CURRENT_SOURCE_DIR}/GenericSolver.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenSimplicialLDLT.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenSimplicialLLT.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenSparseLU.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenSparseQR.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenPardisoLDLT.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenPardisoLLT.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenPardisoLU.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/StaticSolutionState.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/StaticSolutionStateHomogenization.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/GenericSolutionState.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/TransientSolutionNewmark.cpp


		PARENT_SCOPE)

set(HEADERS ${HEADERS}
		${CMAKE_CURRENT_SOURCE_DIR}/GenericSolver.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenSimplicialLDLT.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenSimplicialLLT.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenSparseLU.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenSparseQR.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenPardisoLDLT.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenPardisoLLT.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/EigenPardisoLU.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/SolverTypes.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/StaticSolutionState.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/StaticSolutionStateHomogenization.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/GenericSolutionState.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/SolutionTypes.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/TransientSolution.h
					 ${CMAKE_CURRENT_SOURCE_DIR}/TransientSolutionNewmark.h

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