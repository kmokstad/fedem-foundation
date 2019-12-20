!! $Id$
!! Copyright 2019 SAP SE. All rights reserved.
!!==============================================================================
!> @file splineModule.f90
!>
!> @brief Subroutines for spline curve evaluations.

!!==============================================================================
!> @brief Module with subroutines for spline curve evaluations.
!>
!> @details This module contains subroutines for managing spline curve objects,
!> i.e., adding curves and evaluating them. The actual spline objects are kept
!> in a C++ class FFaSpline with a Fortran wrapper in FFaSplineHandler_F.C,
!> such that the same code can also be used by the GUI.
!>
!> @author Knut Morten Okstad
!>
!> @date 20 Dec 2019

module SplineModule

  implicit none

  integer, parameter :: dp = kind(1.0D0) !< 8-byte real (double)

  !! Global variables used by the call-back subroutine for LMDER
  real(dp), save :: X0(3)   !< Spatial point to project onto the spline curve
  integer , save :: iSpline !< Spline object handle
  integer , save :: iSerr   !< Error flag
  integer , save :: lpUnit  !< Logical print unit

  interface

     !> @brief Creates a new spline curve object.
     !> @param[in] npts Number of control points in the spline curve
     !> @param[in] CtrlPts Spatial coordinates of the spline control points
     !> @param[out] ispl Spline object handle (negative on error)
     subroutine ffa_add_spline (npts,CtrlPts,ispl)
       integer , parameter   :: dp = kind(1.0D0)
       integer , intent(in)  :: npts
       real(dp), intent(in)  :: CtrlPts(3,npts)
       integer , intent(out) :: ispl
     end subroutine ffa_add_spline

     !> @brief Evaluate a spline curve at a parametric point.
     !> @param[in] iSpline Spline object handle
     !> @param[in] iDer Derivative order (0 or 1)
     !> @param[in] s Parameter value of point to evaluate at
     !> @param[out] X Spatial coordinates of point or tangent vector
     !> @param[out] ierr Error flag (negative on error)
     subroutine ffa_eval_spline (iSpline,iDer,s,X,ierr)
       integer , parameter   :: dp = kind(1.0D0)
       integer , intent(in)  :: iSpline, iDer
       real(dp), intent(in)  :: s
       real(dp), intent(out) :: X(3)
       integer , intent(out) :: ierr
     end subroutine ffa_eval_spline

     !> @brief Erases all heap-allocated spline objects.
     subroutine ffa_erase_splines ()
     end subroutine ffa_erase_splines

  end interface


contains

  !!============================================================================
  !> @brief Creates a new spline curve with the given list of control points.
  !>
  !> @param[in] points Spatial coordinates of the spline control points
  !> @return Spline object handle

  function addSpline (points) result(ispl)

    real(dp), intent(in) :: points(:,:)
    integer              :: ispl

    !! --- Logic section ---

    if (size(points,1) /= 3) then
       ispl = -1
    else
       call ffa_add_spline (size(points,2),points,ispl)
    end if

  end function addSpline


  !!============================================================================
  !> @brief Projects a point in space onto a spline curve.
  !>
  !> @param[in] ispl Index (handle) of the spline object to project onto
  !> @param[in] X1 Spatial coordinates of the point to project
  !> @param s Spline parameter of the point (initial guess on input)
  !> @param[out] X2 Spatial coordinates of the projected point
  !> @param[in] nprint Print switch for LMDER
  !> @param[in] lpu Logical print unit
  !> @param[out] ierr Error flag
  !>
  !> @details The point on the spline curve that is closest to the given point
  !> is detected by solving the nonlinear equations F(x) = 0
  !> where the function F is defined by the subroutine SPLINE_FCN.
  !> The Minpack subroutine LMDER from Netlib is used to solve the equation.
  !>
  !> References:
  !> http://www.netlib.org/minpack/
  !> http://www.mcs.anl.gov/~more/ANL8074a.pdf

  subroutine projectOntoCurve (ispl, X1, s, X2, nprint, lpu, ierr)

    integer , intent(in)    :: ispl
    real(dp), intent(in)    :: X1(3)
    real(dp), intent(inout) :: s
    real(dp), intent(out)   :: X2(3)
    integer , intent(in)    :: nprint, lpu
    integer , intent(out)   :: ierr

    !! Local variables
    integer  :: maxfev, mode, ipvt, nfev, njev
    real(dp) :: xtol, factor
    real(dp) :: Fvec, Fjac, diag, qtf, wa1, wa2, wa3, wa4

    !! Objective function
    external :: SPLINE_FCN

    !! --- Logic section ---

    iSpline = ispl
    X0      = X1
    lpUnit  = max(6,lpu)

    !! Set some input parameters for LMDER
    xtol    = 1.0e-8_dp
    maxfev  = 1000
    mode    = 1
    factor  = 100.0_dp

    !! Solve the nonlinear minimization problem using Minpack::LMDER
    call LMDER (SPLINE_FCN,1,1,s,fvec,fjac,1,xtol,xtol,xtol, &
         &      maxfev,diag,mode,factor,nprint,ierr, &
         &      nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4)
    if (iSerr < 0) goto 900

    select case (ierr)
    case (-9:-1)
       if (lpu > 0) write(lpu,*) '*** LMDER: User termination.'
       return
    case (0)
       if (lpu > 0) write(lpu,*) '*** LMDER: Improper input arguments.'
       ierr = -99
       return
    case (1:4)
       if (lpu > 0) write(lpu,*) '  * LMDER has converged.'
    case (5)
       if (lpu > 0) write(lpu,*) ' ** LMDER reached',maxfev,' evaluations.'
    case (6)
       if (lpu > 0) write(lpu,*) ' ** LMDER says ftol is too small.'
    case (7)
       if (lpu > 0) write(lpu,*) ' ** LMDER says xtol is too small.'
    case (8)
       if (lpu > 0) write(lpu,*) ' ** LMDER says gtol is too small.'
    end select

    call ffa_eval_spline (iSpline,0,s,X2,iSerr)
    if (iSerr == 0) return

900 if (lpu > 0) write(lpu,*) '*** Spline evaluation failure',iSerr
    ierr = iSerr

  end subroutine projectOntoCurve

end module SplineModule


!!==============================================================================
!> @brief Evaluates the nonlinear function, F(x), and its Jacobian.
!>
!> @param[in] m Dimension of the funcrtion F(x) (should equal 1 here)
!> @param[in] n Number of unknown variables in @a x (should equal 1 here)
!> @param[in] x Current value of the unknowns to evaluate the function at
!> @param[out] Fvec Function value at @a x
!> @param[out] Fjac Jacobian (gradient) of the function value at @a x
!> @param[in] ldFjac Leading dimension of array @ Fjac
!> @param[in] iflag Evaluation option (1=function value, 2=jacobian)
!>
!> @details The function F(x) defines the distance between a
!> given point in space and a point on a parametric (spline) curve.
!> This subroutine is defined in global scope as it is used as a call-back
!> invoked by the Minpack subroutine LMDER.

subroutine SPLINE_FCN (m, n, x, Fvec, Fjac, ldFjac, iflag)

  use splineModule, only : dp, X0, iSpline, iSerr, lpUnit
  use splineModule, only : ffa_eval_spline

  implicit none

  integer , intent(in)  :: m, n, ldFjac, iflag
  real(dp), intent(in)  :: x(n)
  real(dp), intent(out) :: Fvec(m), Fjac(ldFjac,n)

  !! Local variables
  real(dp) :: Xpnt(3), Xder(3), dist

  !! --- Logic section ---

  select case (iflag)
  case (0) ! Print for current iteration

     write(lpUnit,600) 'X',x
     write(lpUnit,600) 'F',FVec
600  format (A1,':',1P6E12.4 / (2X,1P6E12.4))

  case (1) ! Evaluate the function value

     call ffa_eval_spline (iSpline,0,x(1),Xpnt,iSerr)
     if (iSerr < 0) then
        Fvec = 0.0_dp
     else
        Xpnt = Xpnt - X0
        Fvec = sqrt(dot_product(Xpnt,Xpnt))
     end if

  case (2) ! Evaluate the Jacobian

     call ffa_eval_spline (iSpline,0,x(1),Xpnt,iSerr)
     if (iSerr < 0) then
        Fjac = 0.0_dp
     else
        Xpnt = Xpnt - X0
        dist = sqrt(dot_product(Xpnt,Xpnt))
        call ffa_eval_spline (iSpline,1,x(1),Xder,iSerr)
        if (iSerr < 0 .or. dist < 1.0e-16_dp) then
           Fjac = 0.0_dp
        else
           Fjac = 2.0_dp*dot_product(Xpnt,Xder)/dist
        end if
     end if

  end select

end subroutine SPLINE_FCN
