# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(CDIR ${CMAKE_CURRENT_SOURCE_DIR}/elementFormulations)


set(SOURCEFILES ${SOURCEFILES}
	         ${CMAKE_CURRENT_SOURCE_DIR}/GenericElementFormulation.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL101_Bernoulli2D.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL102_Timoshenko2D.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL103_Timoshenko3D.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL104_TimoshenkoPrism.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL201_2DShell.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL202_Piezo2D.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL203_BeamInterface2D.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL204_BeamInterface2D.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL205_HDivTest.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL206_Plate.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL207_FaceConstraint.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL300_Solid3DLinear.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL301_Piezo3DLinear.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL302_BeamCoupling3D.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL303_ThermoMechanikSolid3D.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL304_QPVolumeElement.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL307_VolumeConstraint.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/LSFEMBernoulli.cpp
			 ${CMAKE_CURRENT_SOURCE_DIR}/EL290_PythonElement.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
	          ${CMAKE_CURRENT_SOURCE_DIR}/GenericElementFormulation.h
	          ${CMAKE_CURRENT_SOURCE_DIR}/GenericElementFormulationInterface.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL101_Bernoulli2D.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL102_Timoshenko2D.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL103_Timoshenko3D.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL104_TimoshenkoPrism.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL201_2DShell.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL202_Piezo2D.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL203_BeamInterface2D.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL204_BeamInterface2D.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL205_HDivTest.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL206_Plate.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL207_FaceConstraint.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL300_Solid3DLinear.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL301_Piezo3DLinear.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL302_BeamCoupling3D.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL303_ThermoMechanikSolid3D.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL304_QPVolumeElement.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL307_VolumeConstraint.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/LSFEMBernoulli.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/EL290_PythonElement.h
		PARENT_SCOPE)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.. ${PREPOCESSOR_INCLUDE} ${Eigen3_DIR})
#
#
#if(WIN32)
#	add_library(CPPFEM_Elementformulations STATIC ${SOURCEFILES} ${HEADERS})
#else()
#	add_library(CPPFEM_Elementformulations SHARED ${SOURCEFILES} ${HEADERS})
#endif()
#
#if(USE_CATALYST)
#target_link_libraries(CPPFEM_Elementformulations PRIVATE ParaView::Catalyst ParaView::RemotingServerManager ParaView::RemotingAnimation ParaView::RemotingLive VTK::CommonDataModel)
#endif()