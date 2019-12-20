# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Build setup

set ( LIB_ID FFaFunctionTests )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

add_executable ( ${LIB_ID}_C funcTest.C )
target_link_libraries ( ${LIB_ID}_C FFaFunctionLib )
if ( USE_FORTRAN )
  add_executable ( ${LIB_ID} funcTest.f90 dummies.f90 )
  target_link_libraries ( ${LIB_ID} FFaFunctionLib_F90 )
  if ( pFUnit_FOUND )
    if ( LINUX )
#     Some of the pFUnit source files might have longer lines than 132 characters
      set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-256" )
    endif ( LINUX )
    add_fortran_test ( test_spline testSpline.pf FFaFunctionLib_F90 )
  endif ( pFUnit_FOUND )
endif ( USE_FORTRAN )
