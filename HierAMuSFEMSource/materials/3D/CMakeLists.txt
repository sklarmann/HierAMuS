# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES}
					 ${CMAKE_CURRENT_SOURCE_DIR}/MA1_LinearElastic_Isotrop.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/MA2_NeoHook.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/MA3_SmallStrainPlasticity.cpp
					 PARENT_SCOPE)

set(HEADERS ${HEADERS} ${CMAKE_CURRENT_SOURCE_DIR}/MA1_LinearElastic_Isotrop.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/MA2_NeoHook.h
					  ${CMAKE_CURRENT_SOURCE_DIR}/MA3_SmallStrainPlasticity.h
					 PARENT_SCOPE)

