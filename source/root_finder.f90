!# Copyright (C) 2018,2019 Luiz Felippe S. Rodrigues, Luke Chamandy
!#
!# This file is part of Magnetizer.
!#
!# Magnetizer is free software: you can redistribute it and/or modify
!# it under the terms of the GNU General Public License as published by
!# the Free Software Foundation, either version 3 of the License, or
!# (at your option) any later version.
!#
!# Magnetizer is distributed in the hope that it will be useful,
!# but WITHOUT ANY WARRANTY; without even the implied warranty of
!# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!# GNU General Public License for more details.
!#
!# You should have received a copy of the GNU General Public License
!# along with Magnetizer.  If not, see <http://www.gnu.org/licenses/>.
!#
!# Significant parts of this file are based on the examples in the
!# documentation of the FGSL library (licensed under GPL2.
!# https://www.lrz.de/services/software/mathematik/gsl/fortran/
!# https://github.com/reinh-bader/fgsl
!#
! Interface to compute roots using (F)GSL

module root_finder
  use FGSL
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: CubicRootClose, FindRoot, FindRoot_deriv, template_function

  integer, save :: error_state = FGSL_success

  ! Auxiliary global module variables
  procedure (template_function), pointer :: Faux => null ()
  real(c_double) :: step_aux

contains

  function dfunc(x, params) bind(c)
    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: dfunc
    integer(FGSL_int) :: status
    type(FGSL_function) :: stdfunc
    real(FGSL_double) :: abserr

    stdfunc = FGSL_function_init(Faux, params)
    status = FGSL_deriv_central(stdfunc, x, step_aux, dfunc, abserr)
  end function

  subroutine fdfunc(x, params, y, dy) bind(c)
    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double), intent(out) :: y, dy

    y = Faux(x, params)
    dy = dfunc(x, params)
  end subroutine fdfunc

  function template_function(x, params) bind(c)
    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: template_function
    template_function = x**2-4-100d0/sin(x+1d-9)
    return
  end function template_function

  double precision function CubicRootClose(a3, a2, a1, a0, refpoint,all_roots)
    ! Finds the roots of the cubic equation: a3*x^3 + a2*x^2 + a1*x + a0 = 0
    ! Returns the root closest to refpoint
    implicit none

    double precision, intent(in) :: a1, a2, a3, a0, refpoint
    double precision, dimension(3) :: roots
    double precision, dimension(3), intent(out), optional :: all_roots
    double precision :: distance, previous_distance
    real(FGSL_double) :: A, B, C
    integer :: number_of_solutions, i

    A = a2/a3
    B = a1/a3
    C = a0/a3

    ! This GSL routine will solve:  x^3 + A*x^2 + B*x + C = 0
    number_of_solutions = FGSL_poly_solve_cubic(A, B, C,   &
                                                roots(1),roots(2),roots(3))
    ! If required, outputs all the roots
    if (present(all_roots)) then
      all_roots=0d0
      do i=1,number_of_solutions
        all_roots(i)=roots(i)
      enddo
    endif

    previous_distance = huge(previous_distance)
    ! Initialises the output to an invalid value
    CubicRootClose = -1d9
    ! Loops over solutions and use the one closest to the reference point
    do i=1,number_of_solutions
      distance = (roots(i)-refpoint)**2
      if (distance < previous_distance .and. roots(i)>0) then
        previous_distance = distance
        CubicRootClose = roots(i)
      endif
    end do
    return
  end function CubicRootClose

  function FindRoot(func, params_ptr, interval, success, max_it) result(root)
    ! Wrapper to FGSL's root finder, based on the FSGL's example roots.f90
    ! Uses FGSL_root_fsolver_brent
    !
    ! Parameters
    ! ----------
    ! func : function
    !     a function with two arguments: a scalar and a pointer to a
    !     C-array of parameters.
    ! params_ptr : c_ptr
    !     C-pointer to an array of parameters (consistent with func)
    ! interval : real(FGSL_double), dimension(2)
    !     Interval containing root (must contain a sign change!)
    ! max_it, optional :
    !     Maximum number of iterations. Default: 50
    ! success, optional : logical
    !     Changes to True if FindRoot converged, False otherwise.
    !
    ! Returns
    ! -------
    ! root : real(FGSL_double)
    !     the root!
    !
    procedure(template_function) :: func
    type(c_ptr), intent(in) :: params_ptr
    real(FGSL_double), dimension(2), intent(in) :: interval
    integer, optional :: max_it
    real(FGSL_double), parameter :: ABS_TOL = 0_FGSL_double
    real(FGSL_double), parameter :: REL_TOL = 1.0e-7_FGSL_double
    integer(FGSL_int) :: itmax = 200
    real(FGSL_double) :: root, xlo, xhi
    character(kind=FGSL_char,len=FGSL_strmax) :: name
    integer :: i
    integer(FGSL_int) :: status
    logical, intent(out), optional :: success
    type(FGSL_root_fsolver) :: root_fslv
    type(FGSL_function) :: FGSL_func

    type(FGSL_error_handler_t) :: default_errh, my_errh

    my_errh = FGSL_error_handler_init(err_sub)
    default_errh = FGSL_set_error_handler(my_errh)

    if (present(success)) success = .true.
    if (present(max_it)) itmax = max_it

    root = -huge(root)
    ! Allocates the solver and initializes the function
    root_fslv = FGSL_root_fsolver_alloc(FGSL_root_fsolver_brent)
    FGSL_func = FGSL_function_init(func, params_ptr)

    ! Iterates to find the root
    if (FGSL_well_defined(root_fslv)) then
      status = FGSL_root_fsolver_set(root_fslv, FGSL_func, interval(1), interval(2))
      name = FGSL_root_fsolver_name (root_fslv)

      if (status /= FGSL_success) then
        if (present(success)) success = .false.
        return
      endif
      i = 0
      do
          i = i + 1
          status = FGSL_root_fsolver_iterate(root_fslv)
          if (status /= FGSL_success .or. i > itmax) then
            if (present(success)) success = .false.
            write(6, *) 'FindRoot: Failed to converge or iterate'
            exit
          end if
          root = FGSL_root_fsolver_root(root_fslv)
          xlo = FGSL_root_fsolver_x_lower(root_fslv)
          xhi = FGSL_root_fsolver_x_upper(root_fslv)
          status = FGSL_root_test_interval(xlo, xhi, ABS_TOL, REL_TOL)
          if (status == FGSL_success) exit
      end do
    end if

    call FGSL_root_fsolver_free(root_fslv)
    call FGSL_function_free(FGSL_func)
  end function FindRoot

  function FindRoot_deriv(func, params_ptr, guess, step, success, max_it) result(root)
    ! Wrapper to FGSL's root finder, based on the FSGL's example roots.f90
    ! Uses gsl_root_fdfsolver_secant and FGSL_deriv_central to estimate the
    ! derivative (which gsl_root_fdfsolver_secant uses in the first step).
    !
    ! Parameters
    ! ----------
    ! func : function
    !     a function with two arguments: a scalar and a pointer to a
    !     C-array of parameters.
    ! params_ptr : c_ptr
    !     C-pointer to an array of parameters (consistent with func)
    ! guess : real(FGSL_double)
    !     Initial guess for finding the root
    ! h : real(FGSL_double)
    !     Step used
    ! max_it, optional :
    !     Maximum number of iterations. Default: 50
    ! success, optional : logical
    !     Changes to True if FindRoot converged, False otherwise.
    !
    ! Returns
    ! -------
    ! root : real(FGSL_double)
    !     the root!
    !
    procedure(template_function) :: func
    type(c_ptr), intent(in) :: params_ptr
    real(FGSL_double), intent(in) :: guess
    real(FGSL_double), intent(in) :: step
    real(FGSL_double), dimension(2) :: interval
    integer, optional :: max_it
    real(FGSL_double), parameter :: ABS_TOL = 0_FGSL_double
    real(FGSL_double), parameter :: REL_TOL = 1.0e-5_FGSL_double
    integer(FGSL_int) :: itmax = 200
    real(FGSL_double) :: ri, root, xlo, xhi

    character(kind=FGSL_char,len=FGSL_strmax) :: name
    integer :: i
    integer(FGSL_int) :: status
    logical, intent(out), optional :: success
    type(FGSL_root_fdfsolver) :: root_fdfslv
    type(FGSL_function_fdf) :: stdfunc_fdf
    type(FGSL_error_handler_t) :: default_errh, my_errh

    Faux => func
    step_aux = step

    my_errh = FGSL_error_handler_init(err_sub)
    default_errh = FGSL_set_error_handler(my_errh)
    if (present(success)) success = .true.
    if (present(max_it)) itmax = max_it

    root = -huge(root)
    ! Allocates the solver and initializes the function
    root_fdfslv = FGSL_root_fdfsolver_alloc(FGSL_root_fdfsolver_secant)
    stdfunc_fdf = FGSL_function_fdf_init(func, dfunc, fdfunc, params_ptr)

    ! Iterates to find the root
    if (FGSL_well_defined(root_fdfslv)) then
      status = FGSL_root_fdfsolver_set(root_fdfslv, stdfunc_fdf, guess)
      name = FGSL_root_fdfsolver_name (root_fdfslv)

      if (status /= FGSL_success) then
        if (present(success)) success = .false.
        return
      endif
      i = 0
      do
        i = i + 1
        status = FGSL_root_fdfsolver_iterate(root_fdfslv)
        if (status /= FGSL_success .or. i > itmax) then
          if (present(success)) success = .false.
          write(6, *) 'FindRoot_deriv: Failed to converge or iterate'
          exit
        end if
        ri = root
        root = FGSL_root_fdfsolver_root(root_fdfslv)
        status = FGSL_root_test_delta (root, ri, ABS_TOL, REL_TOL)
        if (status == FGSL_success) exit
      end do
    end if

    call FGSL_root_fdfsolver_free(root_fdfslv)
    call FGSL_function_fdf_free(stdfunc_fdf)
  end function FindRoot_deriv


  subroutine err_sub (reason, file, line, errno) bind(c)
    ! Here we define our own error handler. It is not thread-safe
    ! since global state is being maintained.
    type(c_ptr), value :: reason, file
    integer(c_int), value :: line, errno
    error_state = errno
  end subroutine err_sub
end module root_finder
