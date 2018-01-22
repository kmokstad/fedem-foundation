// $Id$
////////////////////////////////////////////////////////////////////////////////
//
//    F E D E M    T E C H N O L O G Y   A S
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
  \file FFaSpline.H
  \brief Spline curve representations.
*/

#ifndef FFA_SPLINE_H
#define FFA_SPLINE_H

#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include <vector>
#include <array>


/*!
  \brief Class representing a cubic Catmull-Rom spline curve.
*/

class FFaSpline
{
public:
  //! \brief The constructor create a spline passing through a list of points.
  //! \param[in] points Spatial coordinates of the spline control points
  FFaSpline(const std::vector<FaVec3>& points);
  //! \brief Evaluates the spline at a parametric coordinate
  //! \param[in] xi Parameter of evaluation point, in range [t[1],t[t.size()-2]]
  //! \param[in] deriv Derivative order (0, 1 or 2)
  FaVec3 evaluate(double xi, int deriv = 0) const;
  //! \brief Evaluation operator.
  FaVec3 operator()(double xi) const { return this->evaluate(xi); }
  //! \brief Returns the start parameter of the spline curve.
  double start() const { return t.size() > 2 ? t[1] : 0.0; }
  //! \brief Returns the end parameter of the spline curve.
  double stop() const { return t.size() > 2 ? t[t.size()-2] : 0.0; }
  //! \brief Returns the tangent vector at the specified control point.
  const FaVec3& getTangent(size_t i) const { return tangent[i]; }
  //! \brief Returns the normal (up) vector at the specified control point.
  const FaVec3& getNormal(size_t i) const { return normal[i]; }

private:
  std::vector<double> t; //!< knot vector
  std::vector<FaVec3> P; //!< control points

  std::vector<FaVec3> tangent; //!< tangent vectors at control points
  std::vector<FaVec3> normal ; //!< normal vectors at control points

  mutable std::array<FaVec3,3> A; //!< temporary used in evaluation
  mutable std::array<FaVec3,2> B; //!< temporary used in evaluation
  mutable std::array<FaVec3,3> C; //!< temporary used in evaluation
  mutable std::array<FaVec3,2> D; //!< temporary used in evaluation
  mutable std::array<FaVec3,2> E; //!< temporary used in evaluation
};

#endif
