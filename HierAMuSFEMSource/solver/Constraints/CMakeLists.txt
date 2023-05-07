# Copyright 2022 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(CDIR ${CMAKE_CURRENT_SOURCE_DIR}/equations)


set(SOURCEFILES ${SOURCEFILES}
${CMAKE_CURRENT_SOURCE_DIR}/ConstraintHandler.cpp
${CMAKE_CURRENT_SOURCE_DIR}/BaseConstraint.cpp
${CMAKE_CURRENT_SOURCE_DIR}/GeneralLink.cpp
		PARENT_SCOPE)

set(HEADERS ${HEADERS}
${CMAKE_CURRENT_SOURCE_DIR}/ConstraintHandler.h
${CMAKE_CURRENT_SOURCE_DIR}/BaseConstraint.h
${CMAKE_CURRENT_SOURCE_DIR}/GeneralLink.h
		PARENT_SCOPE)

