# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES ${SOURCEFILES}
					 ${CMAKE_CURRENT_SOURCE_DIR}/MA1_2D_PlainStrain_3D.cpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/MA3_2D_LinearElastic_Isotrop.cpp
					 PARENT_SCOPE)

set(HEADERS ${HEADERS}
					   ${CMAKE_CURRENT_SOURCE_DIR}/MA1_2D_PlainStrain_3D.h
					   ${CMAKE_CURRENT_SOURCE_DIR}/MA3_2D_LinearElastic_Isotrop.h
					 PARENT_SCOPE)

