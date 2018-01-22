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
  \file FFaSpline.C
  \brief Spline curve representations.
*/

#include "FFaLib/FFaAlgebra/FFaSpline.H"


/*!
  See https://en.wikipedia.org/wiki/Centripetal_Catmull-Rom_spline
*/

FFaSpline::FFaSpline (const std::vector<FaVec3>& points)
{
  // Add an extra point at each end since the end segments cannot be evaluated
  P.reserve(2+points.size());
  if (points.size() > 3 && points.front().equals(points.back(),1.0e-10))
  {
    // Looping spline
    P.push_back(points[points.size()-2]);
    P.insert(P.end(),points.begin(),points.end());
    P.push_back(points[1]);
  }
  else if (points.size() > 1)
  {
    P.push_back(2.0*points.front() - points[1]);
    P.insert(P.end(),points.begin(),points.end());
    P.push_back(2.0*points.back() - points[points.size()-2]);
  }
  else
    return; // Too few points

  // Calculate the knots
  size_t i;
  t.resize(P.size(),0.0);
  for (i = 2; i < P.size(); i++)
    t[i] = t[i-1] + sqrt((P[i]-P[i-1]).length()); // centripetal Catmull-Rom
  t.front() = -t[2];

  // Calculate the tangent and normal vector at each point
  const double eps = 1.0e-8;

  size_t nCurve = 0;
  std::vector<bool> straight(points.size(),false);
  tangent.resize(points.size());
  normal.resize(points.size());
  for (i = 0; i < normal.size(); i++)
  {
    tangent[i] = this->evaluate(t[1+i],1);
    tangent[i].normalize();
    double dt = (t[i+2]-t[i])*eps;
    normal[i] = (P[i+1]-P[i]) ^ (P[i+2]-P[i+1]);
    if (normal[i].sqrLength() < dt*dt)
      straight[i] = true;
    else
    {
      normal[i] = tangent[i] ^ normal[i];
      normal[i].normalize();
      nCurve++;
    }
  }

  if (nCurve == 0)
  {
    // The spline is a straight line, choose some arbitrary normal direction
    FaVec3 tan1 = P[2] - P[1];
    if (tan1.isParallell(FaVec3(0.0,0.0,1.0)))
      normal.front() = FaVec3(1.0,0.0,0.0);
    else
    {
      normal.front() = FaVec3(tan1[VY],-tan1[VX],0.0); // tangent x Z-axis
      normal.front().normalize();
    }
    for (i = 1; i < normal.size(); i++)
      normal[i] = normal.front();
  }
  else if (nCurve < normal.size())
  {
    // Some (but not all) segments are straight. Now interpolate the other ones.
    int i0 = -1;
    for (i = 0; i < straight.size(); i++)
      if (straight[i])
      {
        if (i+1 == straight.size() && i0 >= 0) // At last point
          normal[i] = normal[i0];
        else for (size_t i1 = i+1; i1 < straight.size(); i1++)
          if (!straight[i1])
          {
            if (i0 < 0) // At first point
              normal[i] = normal[i1];
            else // Interpolate from the nearest non-straight segments
              normal[i] = (normal[i0]*(t[i+1]-t[i0+1]) +
                           normal[i1]*(t[i1+i]-t[i+1])) / (t[i1+1]-t[i0+1]);
          }

        FaVec3 tangent = this->evaluate(t[1],1);
        FaVec3 binormal = tangent ^ normal[i];
        normal[i] = binormal ^ tangent;
        normal[i].normalize();
      }
      else
        i0 = i;
  }
}


FaVec3 FFaSpline::evaluate (double xi, int deriv) const
{
  if (t.size() < 4 || deriv > 2)
    return FaVec3(); // empty spline

  size_t i = 1, j = t.size()-2, k = 0;
  if (xi < t[i] || xi > t[j])
    return FaVec3(); // outside domain

  // Binary search to find the segment t_i < xi < t_i+1
  while (j-i > 1)
  {
    k = (i+j)/2;
    if (t[k] > xi)
      j = k;
    else if (t[k] < xi)
      i = k;
    else if (deriv <= 0)
      return P[k];
    else
      break;
  }

  // Evaluate the Catmull-Rom spline

  for (j = 0, k = i; j < 3; j++, k++)
    A[j] = ((t[k]-xi)*P[k-1] + (xi-t[k-1])*P[k])/(t[k]-t[k-1]);

  for (j = 0, k = i; j < 2; j++, k++)
    B[j] = ((t[k+1]-xi)*A[j] + (xi-t[k-1])*A[j+1])/(t[k+1]-t[k-1]);

  if (deriv <= 0) // Return the value
    return ((t[i+1]-xi)*B[0] + (xi-t[i])*B[1])/(t[i+1]-t[i]);

  // First derivatives of A and B

  for (j = 0, k = i; j < 3; j++, k++)
    C[j] = (P[k]-P[k-1])/(t[k]-t[k-1]);

  for (j = 0, k = i; j < 2; j++, k++)
    D[j] = ((t[k+1]-xi)*C[j]-A[j] + (xi-t[k-1])*C[j+1]+A[j+1])/(t[k+1]-t[k-1]);

  if (deriv == 1) // Return the first derivative
    return ((t[i+1]-xi)*D[0]-B[0] + (xi-t[i])*D[1]+B[1])/(t[i+1]-t[i]);

  // Second derivatives

  for (j = 0, k = i; j < 2; j++, k++)
    E[j] = 2.0*(C[j+1] + C[j])/(t[k+1]-t[k-1]);

  // Return the second derivative
  return ((t[i+1]-xi)*E[0]-2.0*D[0] + (xi-t[i])*E[1]+2.0*D[1])/(t[i+1]-t[i]);
}
