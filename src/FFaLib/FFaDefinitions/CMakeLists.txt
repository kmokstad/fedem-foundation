# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID FFaDefinitions )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )


## Files with header and source with same name
set ( COMPONENT_FILE_LIST FFaAppInfo FFaListViewItem FFaMsg
                          FFaResultDescription FFaVersionNumber )

## Pure header files, i.e., header files without a corresponding source file
set ( HEADER_FILE_LIST FFaViewItem )

## Pure implementation files, i.e., source files without corresponding header
set ( SOURCE_FILE_LIST FFaViewItemUtils )
if ( USE_FORTRAN )
  list ( APPEND SOURCE_FILE_LIST FFaMsg_F )
endif ( USE_FORTRAN )


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


add_library ( ${LIB_ID} ${CPP_SOURCE_FILES} ${CPP_HEADER_FILES} )
target_link_libraries ( ${LIB_ID} FFaString Admin )
