// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////
/*!
  \file FFaMat33.H
  \brief Point transformations in 3D space.
*/

#ifndef FFA_MAT33_H
#define FFA_MAT33_H

#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include <array>


/*!
  \brief Class for point transformations in 3D space.
  \details This class is used to represent point transformations in geometry
  calculations as 3&times;3 matrices. The class is equipped with methods for
  creating such objects and for the most common operations on them.
*/

class FaMat33
{
public:
  //! \brief Default constructor creating the identity transformation.
  FaMat33() { this->setIdentity(); }
  //! \brief Constructor creating a transformation from an array of floats.
  FaMat33(const float* mat);
  //! \brief Constructor creating a transformation from an array of doubles.
  FaMat33(const double* mat);
  //! \brief Constructor creating a transformation from three (unit) vectors.
  FaMat33(const FaVec3& v0, const FaVec3& v1, const FaVec3& v2);
  //! \brief Copy constructor.
  FaMat33(const FaMat33& m) { v = m.v; }

  //! \brief Assignment operator.
  FaMat33& operator=  (const FaMat33& m) { v = m.v; return *this; }
  //! \brief Addition operator.
  FaMat33& operator+= (const FaMat33& m);
  //! \brief Subtraction operator.
  FaMat33& operator-= (const FaMat33& m);
  //! \brief Scaling operator (multiplication by a scalar).
  FaMat33& operator*= (double d);
  //! \brief Scaling operator (division by a scalar).
  FaMat33& operator/= (double d);

  //! \brief Zero-based basis vector indexing operator (read only).
  const FaVec3& operator[] (int i) const;
  //! \brief Zero-based basis vector indexing operator (for manipulation).
  FaVec3&       operator[] (int i);

  //! \brief One-based indexing operator (read only).
  //! \param[in] i Row index
  //! \param[in] j Column index
  const double& operator() (int i, int j) const;
  //! \brief One-based indexing operator (for manipulation).
  //! \param[in] i Row index
  //! \param[in] j Column index
  double& operator() (int i, int j);

  //! \brief Defines the identity transformation.
  void setIdentity();
  //! \brief Returns the inverse.
  FaMat33 inverse(double eps = 1.0e-16) const;
  //! \brief Returns the transpose.
  FaMat33 transpose() const;
  //! \brief Performs a cyclic permuation of the three basis vectors.
  FaMat33& shift(int delta);

  //! \brief Checks if this transformation equals another one.
  bool isCoincident(const FaMat33& m, double tolerance = 1.0e-10) const;

  //! \brief Computes a globalized coordinate system from X-axis direction.
  FaMat33& makeGlobalizedCS(const FaVec3& v1);
  //! \brief Computes a globalized coordinate system from in-plane directions.
  FaMat33& makeGlobalizedCS(const FaVec3& v1, const FaVec3& v2);
  //! \brief Computes a globalized coordinate system from three points.
  FaMat33& makeGlobalizedCS(const FaVec3& v0,
                            const FaVec3& v1, const FaVec3& v2);
  //! \brief Computes a globalized coordinate system from four points.
  FaMat33& makeGlobalizedCS(const FaVec3& v1, const FaVec3& v2,
                            const FaVec3& v3, const FaVec3& v4);

  //! \brief Computes an incremental rotation tensor from given Euler angles.
  FaMat33& eulerRotateXYZ(const FaVec3& angles);
  //! \brief Computes an incremental rotation tensor from given Euler angles.
  FaMat33& eulerRotateZYX(const FaVec3& angles);
  //! \brief Computes an incremental rotation tensor from given rotation angles.
  FaMat33& incRotate(const FaVec3& angles);

  //! \brief Returns Euler angles corresponding to an incremental rotation.
  FaVec3 getEulerXYZ() const;
  //! \brief Returns Euler angles corresponding to an incremental rotation.
  FaVec3 getEulerZYX() const;
  //! \brief Returns rotation angles corresponding to an incremental rotation.
  FaVec3 getRotation() const;

  //! \brief Returns a transformation defining rotation about the global Z-axis.
  static FaMat33 makeZrotation(double angle);
  //! \brief Returns a transformation defining rotation about the global Y-axis.
  static FaMat33 makeYrotation(double angle);
  //! \brief Returns a transformation defining rotation about the global X-axis.
  static FaMat33 makeXrotation(double angle);

  //! \brief Unary negation operator.
  friend FaMat33  operator-  (const FaMat33& a);
  //! \brief Global addition operator.
  friend FaMat33  operator+  (const FaMat33& a, const FaMat33& b);
  //! \brief Global subtraction operator.
  friend FaMat33  operator-  (const FaMat33& a, const FaMat33& b);

  //! \brief Matrix-matrix multiplication.
  friend FaMat33  operator*  (const FaMat33& a, const FaMat33& b);
  //! \brief Matrix-vector multiplication.
  friend FaVec3   operator*  (const FaMat33& a, const FaVec3&  b);
  //! \brief Scaling by scalar.
  friend FaMat33  operator*  (const FaMat33& a, double d);
  //! \brief Scaling by scalar.
  friend FaMat33  operator*  (double d, const FaMat33& a);
  //! \brief Division by scalar.
  friend FaMat33  operator/  (const FaMat33& a, double d);

  //! \brief Equality operator.
  friend bool     operator== (const FaMat33& a, const FaMat33& b);
  //! \brief Unequality operator.
  friend bool     operator!= (const FaMat33& a, const FaMat33& b);

  //! \brief Printing operator.
  friend std::ostream& operator<< (std::ostream& s, const FaMat33& m);
  //! \brief Reading operator.
  friend std::istream& operator>> (std::istream& s, FaMat33& m);

private:
  std::array<FaVec3,3> v; //!< Three basis vectors defining the transformation
};


// --- inline functions ---

inline const FaVec3& FaMat33::operator[] (int i) const
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 2)
    std::cerr <<"FaMat33::operator[]: index i="<< i <<" is out of range [0,2]"
              << std::endl;
#endif
  return v[i];
}

inline FaVec3& FaMat33::operator[] (int i)
{
#ifdef FFA_INDEXCHECK
  if (i < 0 || i > 2)
    std::cerr <<"FaMat33::operator[]: index i="<< i <<" is out of range [0,2]"
              << std::endl;
#endif
  return v[i];
}


/*!
  This indexing operator accesses the matrix in a Fortran-like manner, where
  the first index is the row number whereas the second index is the column number.
*/

inline const double& FaMat33::operator() (int i, int j) const
{
#ifdef FFA_INDEXCHECK
  if (i < 1 || i > 3)
    std::cerr <<"FaMat33::operator(): index i="<< i <<" is out of range [1,3]"
              << std::endl;
  if (j < 1 || j > 3)
    std::cerr <<"FaMat33::operator(): index j="<< j <<" is out of range [1,3]"
              << std::endl;
#endif
  return v[j-1][i-1];
}

inline double& FaMat33::operator() (int i, int j)
{
#ifdef FFA_INDEXCHECK
   if (i < 1 || i > 3)
     std::cerr <<"FaMat33::operator(): index i="<< i <<" is out of range [1,3]"
               << std::endl;
   if (j < 1 || j > 3)
     std::cerr <<"FaMat33::operator(): index j="<< j <<" is out of range [1,3]"
               << std::endl;
#endif
  return v[j-1][i-1];
}

#endif
