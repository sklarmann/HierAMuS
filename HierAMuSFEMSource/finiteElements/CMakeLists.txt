# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES}
	         ${CMAKE_CURRENT_SOURCE_DIR}/ElementList.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/ElementTypes.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Edge.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Face.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Volume.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/GenericFiniteElement.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/beamInterfaceElement2D.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/beamInterfaceElement3D.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/LinearPrism.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/FaceConstraint.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/VolumeConstraint.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
		${CMAKE_CURRENT_SOURCE_DIR}/ElementList.h
		${CMAKE_CURRENT_SOURCE_DIR}/ElementTypes.h
		${CMAKE_CURRENT_SOURCE_DIR}/Edge.h
		${CMAKE_CURRENT_SOURCE_DIR}/Face.h
		${CMAKE_CURRENT_SOURCE_DIR}/Volume.h
		${CMAKE_CURRENT_SOURCE_DIR}/NormTypes.h
		${CMAKE_CURRENT_SOURCE_DIR}/GenericFiniteElement.h
		${CMAKE_CURRENT_SOURCE_DIR}/beamInterfaceElement2D.h
		${CMAKE_CURRENT_SOURCE_DIR}/beamInterfaceElement3D.h
		${CMAKE_CURRENT_SOURCE_DIR}/LinearPrism.h
		${CMAKE_CURRENT_SOURCE_DIR}/FaceConstraint.h
		${CMAKE_CURRENT_SOURCE_DIR}/VolumeConstraint.h
		PARENT_SCOPE)
 
