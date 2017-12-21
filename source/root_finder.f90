!# Copyright (C) 2018  Luiz Felippe S. Rodrigues, Luke Chamandy
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
! Interface to compute roots using (F)GSL

module root_finder
  use FGSL
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: CubicRootClose, FindRoot

  integer, save :: error_state = fgsl_success

contains
  double precision function CubicRootClose(a3, a2, a1, a0, refpoint,all_roots)
    ! Finds the roots of the cubic equation: a3*x^3 + a2*x^2 + a1*x + a0 = 0
    ! Returns the root closest to refpoint
    implicit none

    double precision, intent(in) :: a1, a2, a3, a0, refpoint
    double precision, dimension(3) :: roots
    double precision, dimension(3), intent(out), optional :: all_roots
    double precision :: distance, previous_distance
    real(fgsl_double) :: A, B, C
    integer :: number_of_solutions, i

    A = a2/a3
    B = a1/a3
    C = a0/a3
    

    ! This GSL routine will solve:  x^3 + A*x^2 + B*x + C = 0
    number_of_solutions = fgsl_poly_solve_cubic(A, B, C,   &
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
    ! Input: func -> a function with two arguments: a scalar and a pointer to a C-array
    !                of parameters.
    !        params -> an array of parameters (consistent with func)
    !        interval -> a 2-array containing the allowed interval
    !        max_it, optional -> maximum number of iterations
    ! Output: the root!
    !
    real(fgsl_double), external :: func
    type(c_ptr), intent(in) :: params_ptr
    real(fgsl_double), dimension(2), intent(in) :: interval
    integer, optional :: max_it
    real(fgsl_double), parameter :: ABS_TOL = 0_fgsl_double
    real(fgsl_double), parameter :: REL_TOL = 1.0e-7_fgsl_double
    integer(fgsl_int) :: itmax = 40
    real(fgsl_double) :: root, xlo, xhi
    character(kind=fgsl_char,len=fgsl_strmax) :: name
    integer :: i
    integer(fgsl_int) :: status
    logical, intent(out), optional :: success
    type(fgsl_root_fsolver) :: root_fslv
    type(fgsl_function) :: fgsl_func

    type(fgsl_error_handler_t) :: default_errh, my_errh

    my_errh = fgsl_error_handler_init(err_sub)
    default_errh = fgsl_set_error_handler(my_errh)

    if (present(success)) success = .true.
    if (present(max_it)) itmax = max_it

    root = -huge(root)
    ! Allocates the solver and initializes the function
    root_fslv = fgsl_root_fsolver_alloc(fgsl_root_fsolver_brent)
    fgsl_func = fgsl_function_init(func, params_ptr)

    ! Iterates to find the root
    if (fgsl_well_defined(root_fslv)) then
      status = fgsl_root_fsolver_set(root_fslv, fgsl_func, interval(1), interval(2))
      name = fgsl_root_fsolver_name (root_fslv)

      if (status /= fgsl_success) then
        if (present(success)) success = .false.
        return
      endif
      i = 0
      do
          i = i + 1
          status = fgsl_root_fsolver_iterate(root_fslv)
          if (status /= fgsl_success .or. i > itmax) then
            if (present(success)) success = .false.
            write(6, *) 'FindRoot: Failed to converge or iterate'
            exit
          end if
          root = fgsl_root_fsolver_root(root_fslv)
          xlo = fgsl_root_fsolver_x_lower(root_fslv)
          xhi = fgsl_root_fsolver_x_upper(root_fslv)
          status = fgsl_root_test_interval(xlo, xhi, ABS_TOL, REL_TOL)
          if (status == fgsl_success) exit
      end do
    end if

    call fgsl_root_fsolver_free(root_fslv)
    call fgsl_function_free(fgsl_func)
  end function FindRoot

  subroutine err_sub (reason, file, line, errno) bind(c)
    ! Here we define our own error handler. It is not thread-safe
    ! since global state is being maintained.
    type(c_ptr), value :: reason, file
    integer(c_int), value :: line, errno
    error_state = errno
  end subroutine err_sub
end module root_finder
