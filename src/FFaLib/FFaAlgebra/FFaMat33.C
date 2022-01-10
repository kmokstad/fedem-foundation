// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////
/*!
  \file FFaMat33.C
  \brief Point transformations in 3D space.
*/

#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"
#include "FFaLib/FFaOS/FFaFortran.H"


/*!
  \class FaMat33

  Matrix layout :
  \code
    [0][0]  [1][0]  [2][0]

    [0][1]  [1][1]  [2][1]

    [0][2]  [1][2]  [2][2]
  \endcode
*/

//! \brief Convenience macro for easy matrix element access.
#define THIS(i,j) this->operator()(i,j)


////////////////////////////////////////////////////////////////////////////////
//
// Constructors
//
////////////////////////////////////////////////////////////////////////////////

FaMat33::FaMat33 (const float* mat)
{
  v[0] = FaVec3(mat);
  v[1] = FaVec3(mat+3);
  v[2] = FaVec3(mat+6);
}


FaMat33::FaMat33 (const double* mat)
{
  v[0] = FaVec3(mat);
  v[1] = FaVec3(mat+3);
  v[2] = FaVec3(mat+6);
}


FaMat33::FaMat33 (const FaVec3& v0, const FaVec3& v1, const FaVec3& v2)
{
  v[0] = v0;
  v[1] = v1;
  v[2] = v2;
}


////////////////////////////////////////////////////////////////////////////////
//
// Local operators
//
////////////////////////////////////////////////////////////////////////////////

FaMat33& FaMat33::operator+= (const FaMat33& m)
{
  v[0] += m.v[0];
  v[1] += m.v[1];
  v[2] += m.v[2];
  return *this;
}


FaMat33& FaMat33::operator-= (const FaMat33& m)
{
  v[0] -= m.v[0];
  v[1] -= m.v[1];
  v[2] -= m.v[2];
  return *this;
}


FaMat33& FaMat33::operator*= (double d)
{
  v[0] *= d;
  v[1] *= d;
  v[2] *= d;
  return *this;
}


FaMat33& FaMat33::operator/= (double d)
{
  if (fabs(d) < EPS_ZERO)
  {
#ifdef FFA_DEBUG
    std::cerr <<"  ** FaMat33::operator/=(double): Division by zero ("
              << d <<")"<< std::endl;
#endif
    v[0] = v[1] = v[2] = FaVec3(HUGE_VAL, HUGE_VAL, HUGE_VAL);
  }
  else
  {
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
  }

  return *this;
}


////////////////////////////////////////////////////////////////////////////////
//
// Special functions
//
////////////////////////////////////////////////////////////////////////////////

FaMat33 FaMat33::inverse (double eps) const
{
  static FaMat33 b;

  double det
    = v[0][0]*(v[1][1]*v[2][2] - v[2][1]*v[1][2])
    - v[0][1]*(v[1][0]*v[2][2] - v[2][0]*v[1][2])
    + v[0][2]*(v[1][0]*v[2][1] - v[2][0]*v[1][1]);

  if (fabs(det) >= eps)
  {
    b.v[0][0] =  (v[1][1]*v[2][2] - v[2][1]*v[1][2]) / det;
    b.v[0][1] = -(v[0][1]*v[2][2] - v[2][1]*v[0][2]) / det;
    b.v[0][2] =  (v[0][1]*v[1][2] - v[1][1]*v[0][2]) / det;
    b.v[1][0] = -(v[1][0]*v[2][2] - v[2][0]*v[1][2]) / det;
    b.v[1][1] =  (v[0][0]*v[2][2] - v[2][0]*v[0][2]) / det;
    b.v[1][2] = -(v[0][0]*v[1][2] - v[1][0]*v[0][2]) / det;
    b.v[2][0] =  (v[1][0]*v[2][1] - v[2][0]*v[1][1]) / det;
    b.v[2][1] = -(v[0][0]*v[2][1] - v[2][0]*v[0][1]) / det;
    b.v[2][2] =  (v[0][0]*v[1][1] - v[1][0]*v[0][1]) / det;
  }
#ifdef FFA_DEBUG
  else
    std::cerr <<"  ** FaMat33::inverse(): Singular matrix, det = "
              << det << std::endl;
#endif

  return b;
}


void FaMat33::setIdentity ()
{
  v[0] = FaVec3(1.0, 0.0, 0.0);
  v[1] = FaVec3(0.0, 1.0, 0.0);
  v[2] = FaVec3(0.0, 0.0, 1.0);
}


FaMat33 FaMat33::transpose () const
{
  return FaMat33(FaVec3(v[0][0], v[1][0], v[2][0]),
                 FaVec3(v[0][1], v[1][1], v[2][1]),
                 FaVec3(v[0][2], v[1][2], v[2][2]));
}


FaMat33& FaMat33::shift (int delta)
{
  if (delta < -2 || delta%3 == 0) return *this;

  // Perform a cyclic permutation of the matrix columns
  FaVec3 v1 = v[0];
  FaVec3 v2 = v[1];
  FaVec3 v3 = v[2];

  v[(3+delta)%3] = v1;
  v[(4+delta)%3] = v2;
  v[(5+delta)%3] = v3;

  return *this;
}


bool FaMat33::isCoincident (const FaMat33& m, double tolerance) const
{
  if (v[0].isParallell(m[0],tolerance) != 1) return false;
  if (v[1].isParallell(m[1],tolerance) != 1) return false;
  if (v[2].isParallell(m[2],tolerance) != 1) return false;
  return true;
}


/*!
  The X-axis of the globalized system is set parallel to the vector \a v1,
  and the two other axes are then set as close as possible to the corresponding
  global coordinate axis directions.
*/

FaMat33& FaMat33::makeGlobalizedCS (const FaVec3& v1)
{
  FaVec3& eX = v[0];
  FaVec3& eY = v[1];
  FaVec3& eZ = v[2];

  eX = v1;
  eX.normalize();

  if (fabs(eX(3)) > fabs(eX(2)))
  {
    // Define eY by projecting the global Y-axis onto the plane defined by eX
    // eY = v1 x (Y x v1)
    eY(1) = -eX(2)*eX(1);
    eY(2) =  eX(1)*eX(1) + eX(3)*eX(3);
    eY(3) = -eX(2)*eX(3);

    eY.normalize();
    eZ = eX ^ eY;
  }
  else
  {
    // Define eZ by projecting the global Z-axis onto the plane defined by eX
    // eZ = v1 x (Z x v1)
    eZ(1) = -eX(3)*eX(1);
    eZ(2) = -eX(3)*eX(2);
    eZ(3) =  eX(1)*eX(1) + eX(2)*eX(2);

    eZ.normalize();
    eY = eZ ^ eX;
  }

  return *this;
}


/*!
  The XY-plane of the globalized system is defined by the two vectors \a v1
  and \a v2, such that its local Z-axis is parallel to the normal vector of
  that plane, and the two other axes are chosen as close as possible to the
  corresponding global coordinate axis directions.
*/

FaMat33& FaMat33::makeGlobalizedCS (const FaVec3& v1, const FaVec3& v2)
{
  FaVec3& eX = v[0];
  FaVec3& eY = v[1];
  FaVec3& eZ = v[2];

  eZ = v1 ^ v2;
  eZ.normalize();

  if (fabs(eZ(1)) < fabs(eZ(2)))
  {
    // Define eX by projecting the global X-axis onto the plane defined by eZ
    // eX = eZ x (X x eZ)
    eX(1) =  eZ(2)*eZ(2) + eZ(3)*eZ(3);
    eX(2) = -eZ(1)*eZ(2);
    eX(3) = -eZ(1)*eZ(3);

    eX.normalize();
    eY = eZ ^ eX;
  }
  else
  {
    // Define eY by projecting the global Y-axis onto the plane defined by eZ
    // eY = eZ x (Y x eZ)
    eY(1) = -eZ(2)*eZ(1);
    eY(2) =  eZ(1)*eZ(1) + eZ(3)*eZ(3);
    eY(3) = -eZ(2)*eZ(3);

    eY.normalize();
    eX = eY ^ eZ;
  }

  return *this;
}

FaMat33& FaMat33::makeGlobalizedCS (const FaVec3& v0,
                                    const FaVec3& v1, const FaVec3& v2)
{
  return makeGlobalizedCS(v1-v0,v2-v0);
}

FaMat33& FaMat33::makeGlobalizedCS (const FaVec3& v1, const FaVec3& v2,
                                    const FaVec3& v3, const FaVec3& v4)
{
  return makeGlobalizedCS(v3-v1,v4-v2);
}


/*!
  This method assumes the rotation angles are applied in the order
  X, Y and Z, i.e., R(x,y,z) = R(z)*R(y)*R(x).
*/

FaMat33& FaMat33::eulerRotateZYX (const FaVec3& angles)
{
  double cx = cos(angles.x());
  double cy = cos(angles.y());
  double cz = cos(angles.z());

  double sx = sin(angles.x());
  double sy = sin(angles.y());
  double sz = sin(angles.z());

  THIS(1,1) =  cz*cy;
  THIS(1,2) =  cz*sy*sx - sz*cx;
  THIS(1,3) =  cz*sy*cx + sz*sx;
  THIS(2,1) =  sz*cy;
  THIS(2,2) =  sz*sy*sx + cz*cx;
  THIS(2,3) =  sz*sy*cx - cz*sx;
  THIS(3,1) = -sy;
  THIS(3,2) =  cy*sx;
  THIS(3,3) =  cy*cx;

  return *this;
}


/*!
  This method is supposed to perform the inverse operation of eulerRotateZYX().
*/

FaVec3 FaMat33::getEulerZYX () const
{
  // Principles from Computer-Aided modelling, dynamic simulation and
  // dimensioning of Mechanisms, Ole Ivar Sivertsen, NTH 1995, pp. 49-51.
  /*
  double  aZ =  atan3(THIS(2,1),THIS(1,1));
  FaMat33 Tb = *this * makeZrotation(aZ);
  double  aY = -atan3(  Tb(3,1),  Tb(1,1));
  FaMat33 Tc =   Tb  * makeYrotation(aY);
  double  aX =  atan3(  Tc(3,2),  Tc(2,2));
  */
#ifdef FFA_DEBUG
  const char* func = "FaMat33::getEulerZYX";
#else
  const char* func = NULL;
#endif
  // New calculation, based on http://www.udel.edu/HNES/HESC427/EULER.DOC
  // See also http://en.wikipedia.org/wiki/Euler_angles
  // See also https://stackoverflow.com/questions/15022630/how-to-calculate-the-angle-from-rotation-matrix
  double aX, aY, aZ, R31 = THIS(3,1);
  if (fabs(R31) < 1.0-EPS_ZERO)
  {
    aZ = atan3(THIS(2,1),THIS(1,1),func);
    aY = atan3(-R31,hypot(THIS(1,1),THIS(2,1)),func);
    aX = atan3(THIS(3,2),THIS(3,3),func);
  }
  else
  {
    aX = atan3(-THIS(2,3),THIS(2,2));
    aY = copysign(0.5*M_PI,R31);
    aZ = 0.0;
  }

  return FaVec3(aX,aY,aZ);
}


/*!
  This method assumes the rotation angles are applied in the order
  Z, Y and X, i.e., R(z,y,x) = R(x)*R(y)*R(z).
*/

FaMat33& FaMat33::eulerRotateXYZ (const FaVec3& angles)
{
  double cx = cos(angles.x());
  double cy = cos(angles.y());
  double cz = cos(angles.z());

  double sx = sin(angles.x());
  double sy = sin(angles.y());
  double sz = sin(angles.z());

  THIS(1,1) =  cy*cz;
  THIS(1,2) = -cy*sz;
  THIS(1,3) =  sy;
  THIS(2,1) =  cx*sz + sx*sy*cz;
  THIS(2,2) =  cx*cz - sx*sy*sz;
  THIS(2,3) = -sx*cy;
  THIS(3,1) = -cx*sy*cz + sx*sz;
  THIS(3,2) =  cx*sy*sz + sx*cz;
  THIS(3,3) =  cx*cy;

  return *this;
}


/*!
  This method is supposed to perform the inverse operation of eulerRotateXYZ().
*/

FaVec3 FaMat33::getEulerXYZ () const
{
#ifdef FFA_DEBUG
  const char* func = "FaMat33::getEulerXYZ";
#else
  const char* func = NULL;
#endif
  double aX, aY, aZ, R13 = THIS(1,3);
  if (fabs(R13) < 1.0-EPS_ZERO)
  {
    aZ = atan3(-THIS(1,2),THIS(1,1),func);
    aY = asin(R13);
    aX = atan3(-THIS(2,3),THIS(3,3),func);
  }
  else
  {
    aX = 0.0;
    aY = copysign(0.5*M_PI,R13);
    aZ = atan3(THIS(2,1),THIS(2,2));
  }

  return FaVec3(aX,aY,aZ);
}


/*!
  This method uses a quaternion representation of the rotation, and is
  equivalent to the Fortran subroutine rotationmodule::vec_to_mat().

  The angles provided are those related to a Rodrigues parameterization.
  Rotation axis, with length equal to the angle to rotate about that axis.
*/

FaMat33& FaMat33::incRotate (const FaVec3& angles)
{
  double theta = angles.length();
  double quat0 = cos(0.5*theta);
  FaVec3 quatr = angles * (theta < EPS_ZERO ? 0.5 : sin(0.5*theta)/theta);
  double quatl = sqrt(quat0*quat0 + quatr.sqrLength());
  quat0 /= quatl;
  quatr /= quatl;

  THIS(1,1) = 2.0*(quatr(1)*quatr(1) + quat0*quat0) - 1.0;
  THIS(2,2) = 2.0*(quatr(2)*quatr(2) + quat0*quat0) - 1.0;
  THIS(3,3) = 2.0*(quatr(3)*quatr(3) + quat0*quat0) - 1.0;

  THIS(1,2) = 2.0*(quatr(1)*quatr(2) - quatr(3)*quat0);
  THIS(1,3) = 2.0*(quatr(1)*quatr(3) + quatr(2)*quat0);
  THIS(2,3) = 2.0*(quatr(2)*quatr(3) - quatr(1)*quat0);

  THIS(2,1) = 2.0*(quatr(2)*quatr(1) + quatr(3)*quat0);
  THIS(3,1) = 2.0*(quatr(3)*quatr(1) - quatr(2)*quat0);
  THIS(3,2) = 2.0*(quatr(3)*quatr(2) + quatr(1)*quat0);

  return *this;
}


/*!
  This method uses a quaternion representation of the rotation, and is
  equivalent to the Fortran subroutine rotationmodule::mat_to_vec().
*/

FaVec3 FaMat33::getRotation () const
{
  int i = 1;
  if (THIS(2,2) > THIS(i,i)) i = 2;
  if (THIS(3,3) > THIS(i,i)) i = 3;

  double quat0;
  FaVec3 quatr;
  double trace = THIS(1,1) + THIS(2,2) + THIS(3,3);
  if (trace > THIS(i,i))
  {
    quat0    = 0.5*sqrt(1.0+trace);
    quatr(1) = (THIS(3,2) - THIS(2,3)) * 0.25/quat0;
    quatr(2) = (THIS(1,3) - THIS(3,1)) * 0.25/quat0;
    quatr(3) = (THIS(2,1) - THIS(1,2)) * 0.25/quat0;
  }
  else
  {
    int j = 1 + i%3;
    int k = 1 + j%3;
    quatr(i) = sqrt(0.5*THIS(i,i) + 0.25*(1.0-trace));
    quat0    = (THIS(k,j) - THIS(j,k)) * 0.25/quatr(i);
    quatr(j) = (THIS(j,i) + THIS(i,j)) * 0.25/quatr(i);
    quatr(k) = (THIS(k,i) + THIS(i,k)) * 0.25/quatr(i);
  }

  double quatl = sqrt(quat0*quat0 + quatr.sqrLength());
  double sthh  = quatr.length() / quatl;
  double cthh  = quat0 / quatl;
  double theta = sthh < 0.7 ? 2.0*asin(sthh) : 2.0*acos(cthh);
  if (theta < EPS_ZERO)
    return quatr * 2.0;
  else if (sthh < 1.0)
    return quatr * theta/sthh;
  else
    return quatr * theta;
}


FaMat33 FaMat33::makeZrotation (double angle)
{
  static FaMat33 r;

  double c = cos(angle);
  double s = sin(angle);

  r(1,1) = c;
  r(2,1) = s;
  r(1,2) = -s;
  r(2,2) = c;

  return r;
}

FaMat33 FaMat33::makeYrotation (double angle)
{
  static FaMat33 r;

  double c = cos(angle);
  double s = sin(angle);

  r(1,1) = c;
  r(3,1) = -s;
  r(1,3) = s;
  r(3,3) = c;

  return r;
}

FaMat33 FaMat33::makeXrotation (double angle)
{
  static FaMat33 r;

  double c = cos(angle);
  double s = sin(angle);

  r(2,2) = c;
  r(3,2) = s;
  r(2,3) = -s;
  r(3,3) = c;

  return r;
}


////////////////////////////////////////////////////////////////////////////////
//
// Global operators
//
////////////////////////////////////////////////////////////////////////////////
//! \cond DO_NOT_DOCUMENT because already documented in the header file

FaMat33 operator- (const FaMat33& a)
{
  return FaMat33(-a.v[0], -a.v[1], -a.v[2]);
}


FaMat33 operator+ (const FaMat33& a, const FaMat33& b)
{
  return FaMat33(a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2]);
}


FaMat33 operator- (const FaMat33& a, const FaMat33& b)
{
  return FaMat33(a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]);
}


FaMat33 operator* (const FaMat33& a, const FaMat33& b)
{
  return FaMat33(a*b.v[0], a*b.v[1], a*b.v[2]);
}


FaMat33 operator* (const FaMat33& a, double d)
{
  return FaMat33(a.v[0]*d, a.v[1]*d, a.v[2]*d);
}


FaMat33 operator* (double d, const FaMat33& a)
{
  return FaMat33(a.v[0]*d, a.v[1]*d, a.v[2]*d);
}


FaVec3 operator* (const FaMat33& m, const FaVec3& v1)
{
  return m.v[0]*v1[0] + m.v[1]*v1[1] + m.v[2]*v1[2];
}


FaMat33 operator/ (const FaMat33& a, double d)
{
  if (fabs(d) < EPS_ZERO)
  {
#ifdef FFA_DEBUG
    std::cerr <<"  ** FaMat33 operator/(FaMat33&,double): Division by zero ("
              << d <<")"<< std::endl;
#endif
    FaVec3 huge_vec(HUGE_VAL, HUGE_VAL, HUGE_VAL);
    return FaMat33(huge_vec, huge_vec, huge_vec);
  }

  return FaMat33(a.v[0] / d, a.v[1] / d, a.v[2] / d);
}


bool operator== (const FaMat33& a, const FaMat33& b)
{
  return (a.v[0] == b.v[0]) && (a.v[1] == b.v[1]) && (a.v[2] == b.v[2]);
}

bool operator!= (const FaMat33& a, const FaMat33& b)
{
  return !(a == b);
}


std::ostream& operator<< (std::ostream& s, const FaMat33& m)
{
  return s <<'\n'<< m.v[0][0] <<' '<< m.v[1][0] <<' '<< m.v[2][0]
	   <<'\n'<< m.v[0][1] <<' '<< m.v[1][1] <<' '<< m.v[2][1]
	   <<'\n'<< m.v[0][2] <<' '<< m.v[1][2] <<' '<< m.v[2][2];
}

std::istream& operator>> (std::istream& s, FaMat33& m)
{
  FaMat33 m_tmp;
  s >> m_tmp.v[0][0] >> m_tmp.v[1][0] >> m_tmp.v[2][0]
    >> m_tmp.v[0][1] >> m_tmp.v[1][1] >> m_tmp.v[2][1]
    >> m_tmp.v[0][2] >> m_tmp.v[1][2] >> m_tmp.v[2][2];
  if (s) m = m_tmp;
  return s;
}


////////////////////////////////////////////////////////////////////////////////
//
// FORTRAN interface to selected functions
//
////////////////////////////////////////////////////////////////////////////////

SUBROUTINE(ffa_eulerzyx,FFA_EULERZYX) (double* a, double* b, double* angles)
{
  FaVec3 euler(FaMat33(FaMat33(a).transpose() * FaMat33(b)).getEulerZYX());
  angles[0] = euler[0];
  angles[1] = euler[1];
  angles[2] = euler[2];
}


SUBROUTINE(ffa_glbeulerzyx,FFA_GLBEULERZYX) (double* a, double* angles)
{
  FaVec3 euler(FaMat33(a).getEulerZYX());
  angles[0] = euler[0];
  angles[1] = euler[1];
  angles[2] = euler[2];
}

//! \endcond
