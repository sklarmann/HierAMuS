# Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
#
# SPDX-License-Identifier: BSD-3-Clause

set(SOURCEFILES)
set(HEADRS)

add_subdirectory(plot)

option(BUILD_PLOT "Enable build of plotting library, requires Paraview" OFF)


if(BUILD_PLOT)
#
#
include(${CMAKE_SOURCE_DIR}/external/installParaviewLib.cmake)
find_package(ParaView  REQUIRED
HINTS "${STAGED_INSTALL_PREFIX}/lib/cmake/paraview-5.10/")

if(WIN32)
	set(HELPSOURCE ${SOURCEFILES} ${HEADERS})
	add_library(HierAMuS_Plot SHARED ${HELPSOURCE})
else()
	add_library(HierAMuS_Plot SHARED ${SOURCEFILES})
endif()
#

target_link_libraries(HierAMuS_Plot PUBLIC Eigen3::Eigen)
set_target_properties(HierAMuS_Plot PROPERTIES DEBUG_POSTFIX "d")

if(WIN32)

else()
    target_link_options(HierAMuS_Plot PRIVATE "-Wl,--exclude-libs,ALL")
endif()
target_compile_definitions(HierAMuS_Plot PRIVATE BUILDCPPFEMFEM)


if(WIN32)
foreach(FILE ${HELPSOURCE})
    # Get the directory of the source file
    get_filename_component(PARENT_DIR "${FILE}" DIRECTORY)

    # Remove common directory prefix to make the group
    string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}" "" GROUP "${PARENT_DIR}")

    # Make sure we are using windows slashes
    string(REPLACE "/" "\\" GROUP "${GROUP}")

    # Group into "Source Files" and "Header Files"
    if ("${FILE}" MATCHES ".*\\.cpp")
       set(GROUP "Source Files${GROUP}")
    elseif("${FILE}" MATCHES ".*\\.h")
       set(GROUP "Header Files${GROUP}")
    endif()

    source_group("${GROUP}" FILES "${FILE}")
endforeach()
endif()



if(USE_CATALYST)
target_link_libraries(HierAMuS_Plot PRIVATE ParaView::Catalyst ParaView::RemotingServerManager ParaView::RemotingLive ParaView::RemotingAnimation VTK::CommonDataModel)
endif()
target_link_libraries(HierAMuS_Plot PRIVATE Eigen3::Eigen Spectra::Spectra)



target_include_directories(HierAMuS_Plot PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

if(USE_MKL)
target_link_libraries(HierAMuS_Plot PRIVATE MKL::MKL)

endif()
set(CPPYHEADERS ${HEADERS} PARENT_SCOPE)
set(SOURCEFILES ${SOURCEFILES} PARENT_SCOPE)



install(TARGETS HierAMuS_Plot EXPORT HierAMuS_PlotTargets
  LIBRARY DESTINATION .
  ARCHIVE DESTINATION .
  RUNTIME DESTINATION .
  INCLUDES DESTINATION HierAMuS_Plot/include
)

#install(
#  FILES
#    ${HEADERS}
#    #"${CMAKE_CURRENT_BINARY_DIR}/testlib_export.h"
#  DESTINATION
#    include
#  COMPONENT
#    Devel
#)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/" # source directory
        DESTINATION "HierAMuS_Plot/include" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.h" # select header files
        )
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../types" # source directory
        DESTINATION "HierAMuS_Plot/include/" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.h" # select header files
        )

install(
  FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/../datatypes.h

#    #"${CMAKE_CURRENT_BINARY_DIR}/testlib_export.h"
  DESTINATION
    HierAMuS_Plot/include
#  COMPONENT
#    Devel
)


include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_Plot/HierAMuS_PlotConfigVersion.cmake"
  VERSION ${VERSION}
  COMPATIBILITY AnyNewerVersion
)

export(EXPORT HierAMuS_PlotTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_Plot/HierAMuS_PlotTargets.cmake"
  NAMESPACE CPPFEM::
)
configure_file(cmake/HierAMuS_PlotConfig.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_Plot/HierAMuS_PlotConfig.cmake"
  COPYONLY
)

set(ConfigPackageLocation lib/cmake/HierAMuS_Plot)
install(EXPORT HierAMuS_PlotTargets
  FILE
    HierAMuS_PlotTargets.cmake
  NAMESPACE
    CPPFEM::
  DESTINATION
    ${ConfigPackageLocation}
)
install(
  FILES
    cmake/HierAMuS_PlotConfig.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/HierAMuS_Plot/HierAMuS_PlotConfigVersion.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)

add_custom_command(TARGET HierAMuS_Plot POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:HierAMuS_Plot> ${CMAKE_SOURCE_DIR}/pybind/HierAMuS/${targetfile}
)


else()
  
  #find_package(HierAMuS_Plot REQUIRED)
  add_library(HierAMuS_Plot SHARED IMPORTED GLOBAL)
  if(WIN32)
    set_property(TARGET HierAMuS_Plot PROPERTY IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/precompiled/plot/windows/HierAMuS_Plot.dll")
    set_property(TARGET HierAMuS_Plot PROPERTY IMPORTED_IMPLIB "${CMAKE_SOURCE_DIR}/precompiled/plot/windows/HierAMuS_Plot.lib")
    #file(COPY
    #${CMAKE_SOURCE_DIR}/precompiled/plot/windows/HierAMuS_Plot.dll
    #DESTINATION ${CMAKE_INSTALL_PREFIX}/)
    #file(COPY
    #${CMAKE_SOURCE_DIR}/precompiled/plot/windows/catalyst.dll
    #DESTINATION ${CMAKE_INSTALL_PREFIX}/)
    
    file(COPY
    ${CMAKE_SOURCE_DIR}/precompiled/plot/windows/HierAMuS_Plot.dll
    DESTINATION ${CMAKE_SOURCE_DIR}/pybind/HierAMuS)
    file(COPY
    ${CMAKE_SOURCE_DIR}/precompiled/plot/windows/catalyst.dll
    DESTINATION ${CMAKE_SOURCE_DIR}/pybind/HierAMuS)
    
  elseif(UNIX)
    set_property(TARGET HierAMuS_Plot PROPERTY IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/precompiled/plot/linux/libHierAMuS_Plot.so")
    #file(COPY
    #${CMAKE_SOURCE_DIR}/precompiled/plot/linux/libHierAMuS_Plot.so
    #DESTINATION ${CMAKE_INSTALL_PREFIX}/)
    #file(COPY
    #${CMAKE_SOURCE_DIR}/precompiled/plot/linux/libcatalyst.so.3
    #DESTINATION ${CMAKE_INSTALL_PREFIX}/)
    #file(COPY
    #${CMAKE_SOURCE_DIR}/precompiled/plot/linux/libcatalyst.so
    #DESTINATION ${CMAKE_INSTALL_PREFIX}/)

    
    file(COPY
    ${CMAKE_SOURCE_DIR}/precompiled/plot/linux/libHierAMuS_Plot.so
    DESTINATION ${CMAKE_SOURCE_DIR}/pybind/HierAMuS)
    file(COPY
    ${CMAKE_SOURCE_DIR}/precompiled/plot/linux/libcatalyst.so.3
    DESTINATION ${CMAKE_SOURCE_DIR}/pybind/HierAMuS)
    file(COPY
    ${CMAKE_SOURCE_DIR}/precompiled/plot/linux/libcatalyst.so
    DESTINATION ${CMAKE_SOURCE_DIR}/pybind/HierAMuS)
  endif()
  #message(${CMAKE_CURRENT_SOURCE_DIR}/plot)
  #set_property(TARGET HierAMuS_Plot PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/plot)
  #set_property(TARGET HierAMuS_Plot PROPERTY BINARY_DIR ${CMAKE_SOURCE_DIR}/precompiled/plot/windows/)
  
endif()

