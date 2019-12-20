// $Id$
////////////////////////////////////////////////////////////////////////////////
//
//    F E D E M    T E C H N O L O G Y    A S
//
//    Copyright (C)
//    1999 - 2018
//    FEDEM Technology AS
//    all rights reserved
//
//    This is UNPUBLISHED PROPRIETARY SOURCE CODE of FEDEM Technology AS;
//    the contents of this file may not be disclosed to third parties,
//    copied or duplicated in any form, in whole or in part, without
//    the prior written permission of FEDEM Technology AS.
//
////////////////////////////////////////////////////////////////////////////////
/*!
  \file FFaSplineHandler_F.C
  \brief Fortran wrapper for the FFaSpline class.
*/

#include <cstring>

#include "FFaLib/FFaAlgebra/FFaSpline.H"
#include "FFaLib/FFaOS/FFaFortran.H"

static std::vector<FFaSpline*> ourSplines; //!< Spline object container


/*!
  \brief Static helper checking that the \a splineIndex is in valid range.
*/

static bool checkSplineIndex (int splineIndex, int& ierr)
{
  int nSpline = ourSplines.size();
  if (splineIndex >= 0 && splineIndex < nSpline)
  {
    ierr = 0;
    return true;
  }
  else
  {
    std::cerr <<" *** Spline index "<< splineIndex <<" out of range [0,"
              << nSpline-1 <<"]."<< std::endl;
    ierr = -1;
    return false;
  }
}


/*!
  \brief Evaluates a spline object at a specified point.
*/

SUBROUTINE(ffa_eval_spline,FFA_EVAL_SPLINE) (const int& splineIndex,
                                             const int& derivOrder,
                                             const double& splinePrm,
                                             double* Xpt, int& ierr)
{
  if (checkSplineIndex(splineIndex,ierr))
  {
    FaVec3 Xs = ourSplines[splineIndex]->evaluate(splinePrm,derivOrder);
    memcpy(Xpt,Xs.getPt(),3*sizeof(double));
#if FFA_DEBUG > 2
    std::cout <<"ffa_eval_spline("<< splinePrm
              <<","<< derivOrder <<"): "<< Xs << std::endl;
#endif
  }
}


/*!
  \brief Erases all heap-allocated spline object.
*/

SUBROUTINE(ffa_add_spline,FFA_ADD_SPLINE) (const int& npts,
                                           const double* points, int& splineIdx)
{
  if (npts < 1)
  {
    std::cerr <<" *** Too few points "<< npts
              <<" in spline object."<< std::endl;
    splineIdx = -2;
    return;
  }

  std::vector<FaVec3> pointVec;
  pointVec.reserve(npts);
  for (int i = 0; i < npts; i++)
    pointVec.push_back(FaVec3(points+3*i));

  splineIdx = ourSplines.size();
  ourSplines.push_back(new FFaSpline(pointVec));
  std::cout <<"Parameter range ["<< ourSplines.back()->start()
            <<","<< ourSplines.back()->stop() <<"]"<< std::endl;
}


/*!
  \brief Erases all heap-allocated spline objects.
*/

SUBROUTINE(ffa_erase_splines,FFA_ERASE_SPLINES) ()
{
  for (FFaSpline* spline : ourSplines) delete spline;
  ourSplines.clear();
}
