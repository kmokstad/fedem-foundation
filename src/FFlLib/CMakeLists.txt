# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID FFlLib )
set ( LIB_ID_LIST FFlFEParts FFlIOAdaptors )
if ( "${APPLICATION_ID}" STREQUAL "fedemKernel" AND BUILD_TESTS )
  list ( APPEND LIB_ID_LIST FFlTests )
endif ( "${APPLICATION_ID}" STREQUAL "fedemKernel" AND BUILD_TESTS )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

set ( USE_STRAINCOAT TRUE )

if ( USE_VTFAPI )
  find_package ( VTFAPI )
endif ( USE_VTFAPI )

if ( USE_REFCOUNT )
  string ( APPEND CMAKE_CXX_FLAGS " -DFFL_REFCOUNT" )
endif ( USE_REFCOUNT )

if ( USE_PROFILER )
  string ( APPEND CMAKE_CXX_FLAGS " -DFFL_TIMER" )
endif ( USE_PROFILER )

if ( USE_STRAINCOAT )
  string ( APPEND CMAKE_CXX_FLAGS " -DFT_USE_STRAINCOAT" )
endif ( USE_STRAINCOAT )

if ( USE_VERTEXOBJ )
  list ( INSERT LIB_ID_LIST 2 FFlVisualization )
  string ( APPEND CMAKE_CXX_FLAGS " -DFT_USE_VERTEX" )
endif ( USE_VERTEXOBJ )

if ( USE_VISUALS )
  string ( APPEND CMAKE_CXX_FLAGS " -DFT_USE_VISUALS" )
endif ( USE_VISUALS )

if ( USE_CONNECTORS )
  string ( APPEND CMAKE_CXX_FLAGS " -DFT_USE_CONNECTORS" )
endif ( USE_CONNECTORS )

if ( USE_MEMPOOL )
  string ( APPEND CMAKE_CXX_FLAGS " -DFT_USE_MEMPOOL" )
endif ( USE_MEMPOOL )

if ( FT_DEBUG GREATER 12 )
  math ( EXPR FFL_DEBUG "${FT_DEBUG} - 10" )
else ( FT_DEBUG GREATER 12 )
  set ( FFL_DEBUG 2 )
endif ( FT_DEBUG GREATER 12 )
string ( APPEND CMAKE_CXX_FLAGS_DEBUG " -DFFL_DEBUG=${FFL_DEBUG}" )


# Subfolder handling

foreach ( FOLDER ${LIB_ID_LIST} )
  add_subdirectory ( ${FOLDER} )
endforeach ( FOLDER ${LIB_ID_LIST} )


## Files with header and source with same name
set ( COMPONENT_FILE_LIST FFlAttributeBase FFlConnectorItems FFlElementBase
                          FFlFEAttributeRefs FFlFEAttributeSpec
                          FFlFEElementTopSpec FFlFENodeRefs FFlFieldBase
                          FFlGroup FFlLinkHandler FFlLoadBase FFlMemPool
                          FFlNamedPartBase FFlPartBase FFlUtils
)
if ( USE_VERTEXOBJ )
  list ( APPEND COMPONENT_FILE_LIST FFlVertex )
endif ( USE_VERTEXOBJ )
if ( USE_VISUALS )
  list ( APPEND COMPONENT_FILE_LIST FFlVisualBase FFlVisualRefs )
endif ( USE_VISUALS )

## Pure header files, i.e., header files without a corresponding source file
set ( HEADER_FILE_LIST FFlFEResultBase FFlField FFlInit
                       FFlLinkCSMask FFlReference FFlTypeInfoSpec
)

## Pure implementation files, i.e., source files without corresponding header
set ( SOURCE_FILE_LIST )
if ( "${APPLICATION_ID}" STREQUAL "fedemKernel" AND USE_FORTRAN )
  set ( SOURCE_FILE_LIST FFlLinkHandler_F )
endif ( "${APPLICATION_ID}" STREQUAL "fedemKernel" AND USE_FORTRAN )

## Fortran 90 source files
set ( F90_FILE_LIST FFlLinkHandlerInterface )


foreach ( FILE ${COMPONENT_FILE_LIST} )
  list ( APPEND CPP_SOURCE_FILES ${FILE}.C )
  list ( APPEND CPP_HEADER_FILES ${FILE}.H )
endforeach ( FILE ${COMPONENT_FILE_LIST} )

foreach ( FILE ${HEADER_FILE_LIST} )
  list ( APPEND CPP_HEADER_FILES ${FILE}.H )
endforeach ( FILE ${HEADER_FILE_LIST} )

foreach ( FILE ${SOURCE_FILE_LIST} )
  list ( APPEND CPP_SOURCE_FILES ${FILE}.C )
endforeach ( FILE ${SOURCE_FILE_LIST} )

foreach ( FILE ${F90_FILE_LIST} )
  list ( APPEND F90_SOURCE_FILES ${FILE}.f90 )
endforeach ( FILE ${F90_FILE_LIST} )


set ( DEPENDENCY_LIST FFlFEParts FFaAlgebra FFaDefinitions )

if ( "${APPLICATION_ID}" STREQUAL "fedemKernel" AND USE_FORTRAN )
  list ( INSERT DEPENDENCY_LIST 0 FFlIOAdaptors )
endif ( "${APPLICATION_ID}" STREQUAL "fedemKernel" AND USE_FORTRAN )

if ( USE_MEMPOOL )
  list ( APPEND DEPENDENCY_LIST FFaPatterns )
endif ( USE_MEMPOOL )

if ( USE_PROFILER )
  list ( APPEND DEPENDENCY_LIST FFaProfiler )
endif ( USE_PROFILER )

if ( VTFAPI_FOUND )
  list ( APPEND DEPENDENCY_LIST ${VTF_LIBRARY} )
endif ( VTFAPI_FOUND )

add_library ( ${LIB_ID} ${CPP_SOURCE_FILES} ${CPP_HEADER_FILES} )
target_link_libraries ( ${LIB_ID} ${DEPENDENCY_LIST} )

if ( "${APPLICATION_ID}" STREQUAL "fedemKernel" AND USE_FORTRAN )
  add_library ( ${LIB_ID}_F ${F90_SOURCE_FILES} )
  target_link_libraries ( ${LIB_ID}_F ${LIB_ID} FFaCmdLineArg )
endif ( "${APPLICATION_ID}" STREQUAL "fedemKernel" AND USE_FORTRAN )
