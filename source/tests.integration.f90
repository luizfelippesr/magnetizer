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
!! (Slightly) adapted from Galacticus (https://bitbucket.org/abensonca/galacticus/)
!!
!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Contains a program to test integration routines.

module Test_Integration_Functions
  use fgsl
  use, intrinsic :: iso_c_binding

  type   (fgsl_function             ), save, private :: integrandFunction
  type   (fgsl_integration_workspace), save, private :: integrationWorkspace
  real(fgsl_double), parameter :: exact = 1.3932039296856768_fgsl_double
  integer(fgsl_size_t) :: calls = 1000
contains
  function g(v_c, n, params) bind(c)
    integer(c_size_t), value :: n
    type(c_ptr), value :: v_c, params
    real(c_double) :: g
    real(c_double), dimension(:), pointer :: v
    call c_f_pointer(v_c, v, [n])
    g = 1.0_c_double/(m_pi**3)/ &
         (1.0_fgsl_double - cos(v(1))*cos(v(2))*cos(v(3)))
  end function g
  subroutine display_results(title,result,error)
    character(kind=fgsl_char,len=*), intent(in) :: title
    real(fgsl_double), intent(in) :: result, error
    write(6, '(A,'' ================='')') trim(title)
    write(6, '(''result = '',F10.6)') result
    write(6, '(''sigma  = '',F10.6)') error
    write(6, '(''exact  = '',F10.6)') exact
    write(6, '(''error  = '',F10.6,'' = '',F5.1,'' sigma'')') &
         result - exact, abs(result - exact) / error
  end subroutine display_results
end module
program Test_Integration
  use Integration
  use Test_Integration_Functions
  implicit none
  real(fgsl_double) :: chisq, xl(3), xu(3), res, err, y, yy, yyy
  integer(fgsl_int) :: status, stage, mode, verbose
  integer(fgsl_size_t) :: its
  type(fgsl_file) :: file
  type(fgsl_rng) :: r
  type(fgsl_rng_type) :: t
  type(fgsl_monte_plain_state) :: s
  type(fgsl_monte_miser_state) :: m
  type(fgsl_monte_vegas_state) :: v
  type(fgsl_monte_function) :: gfun
  type(c_ptr) :: ptr
!
  xl = 0.0_fgsl_double ; xu = m_pi
  t = fgsl_rng_env_setup()
  r = fgsl_rng_alloc(t)
  gfun = fgsl_monte_function_init(g, 3_fgsl_size_t, ptr)
  s = fgsl_monte_plain_alloc(3_fgsl_size_t)
  status = fgsl_monte_plain_integrate(gfun, xl, xu, 3_fgsl_size_t, &
       calls, r, s, res, err)
  print *, calls
  call fgsl_monte_plain_free(s)
  call display_results('plain',res,err)
!
  m = fgsl_monte_miser_alloc(3_fgsl_size_t)
  status = fgsl_monte_miser_integrate(gfun, xl, xu, 3_fgsl_size_t, &
       calls, r, m, res, err)
  call fgsl_monte_miser_free(m)
  call display_results('miser',res,err)
!
  calls = 1500
  v = fgsl_monte_vegas_alloc(3_fgsl_size_t)
  status = fgsl_monte_vegas_integrate(gfun, xl, xu, 3_fgsl_size_t, &
       calls, r, v, res, err)
  call display_results('vegas warm-up',res,err)
  write(6, *) 'converging ...'
  calls =2000
  do
    print *, '?'
     status = fgsl_monte_vegas_integrate(gfun, xl, xu, 3_fgsl_size_t, &
          calls/5_fgsl_size_t, r, v, res, err)
     call fgsl_monte_vegas_getparams(v, y, yy, chisq, yyy, &
       its, stage, mode, verbose, file)
     write(6, '(''result = '',F10.6,'' sigma = '',F10.6, &
          & '' chisq/dof = '',F6.1)') res,err,chisq
     if (abs(chisq - 1.0_fgsl_double) <= 0.5_fgsl_double) exit
  end do
  call display_results('vegas converged',res,err)
  call fgsl_monte_vegas_free(v)
  call fgsl_monte_function_free(gfun)
  call fgsl_rng_free(r)

!   !% Tests that numerical integration routines work.
!   use Integration
!   use Test_Integration_Functions
!   implicit none
!   type            (fgsl_function             ) :: integrandFunction
!   type            (fgsl_integration_workspace) :: integrationWorkspace
!   logical                                      :: integrationReset
!   double precision                             :: integral, expected, rel_err
!   integer i
!   double precision :: cpu_time_start, cpu_time_finish
!
!   integrationReset=.true.
! !   ! Test simple integrations.
!   integral=Integrate(0.0d0, 1.0d0, Integrand1, integrandFunction&
!        &,integrationWorkspace,toleranceRelative=1.0d-6, reset=integrationReset)
!   print *, "integrate f(x)=x          from x=0……1"
!   expected = 0.5d0
!   rel_err = abs((integral-expected)/expected)
!   print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
!   print *,
!
! !   integrationReset=.false.
!   integral=Integrate(0.0d0,2.0d0*m_pi,Integrand2,integrandFunction&
!        &,integrationWorkspace,toleranceRelative=1.0d-6, reset=integrationReset)
!   print *, "integrate f(x)=sin(x)     from x=0…2π"
!   expected = 0.0d0
!   print *, '  result=',integral,'  expected=', expected
!   print *,
!
!   integral=Integrate(0.0d0,10.0d0,Integrand3,integrandFunction&
!        &,integrationWorkspace,toleranceRelative=1.0d-6, reset=integrationReset)
!   print *, "integrate f(x)=1/√x       from x=0…10"
!   expected = 2.0d0*sqrt(10.0d0)
!   rel_err = abs((integral-expected)/expected)
!   print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
!   print *,
!
!   ! Test 2D integrations.
!   integral=Integrate(0.0d0,2.0d0*m_pi,Integrand4,integrandFunction&
!        &,integrationWorkspace,toleranceRelative=1.0d-6, reset=integrationReset)
!   print *, "integrate f(x,y)=y·cos(x) from x=0…2π and y=0…x"
!   expected = 2.0d0*m_pi
!   rel_err = abs((integral-expected)/expected)
!   print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
!   print *,
!
!
!   integrationReset = .true.
!   call cpu_time(cpu_time_start)
!   do i = 1, 100000
!     integral=Integrate(1d-20,-1.0d0,test_Pdm, integrandFunction,   &
!                        integrationWorkspace, toleranceRelative=1.0d-7, &
!                        toInfinity=.true., reset=integrationReset)
!   end do
!   call cpu_time(cpu_time_finish)
!   print *, "integrate f(x)=tanh(x)/(x/30)/(1+x/30)²       from x=0…∞"
!   expected = 97.946378405461
!   rel_err = abs((integral-expected)/expected)
!   print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
!   print *, '  time: ',cpu_time_finish-cpu_time_start, '(100,000 repetitions)'
!   print *,
!
!   integrationReset = .true.
!   call cpu_time(cpu_time_start)
!   do i = 1, 100000
!     integral=Integrate(0.0d0,100.0d0,test_Pdm, integrandFunction,  &
!                        integrationWorkspace, toleranceRelative=1.0d-7, &
!                        hasSingularities=.true.,toInfinity=.false.,     &
!                        reset=integrationReset)
!   end do
!   call cpu_time(cpu_time_finish)
!   print *, "integrate f(x)=tanh(x)/(x/30)/(1+x/30)²       from x=0…100"
!   expected = 97.946378405461
!   rel_err = abs((integral-expected)/expected)
!   print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
!   print *, '  time: ',cpu_time_finish-cpu_time_start, '(100,000 repetitions)'
!   print *,
!
!   integrationReset = .true.
!   call cpu_time(cpu_time_start)
!   do i = 1, 100000
!     integral=Integrate(1d-10,-1.0d0,test_Pstars, integrandFunction,   &
!                        integrationWorkspace, toleranceRelative=1.0d-7, &
!                        toInfinity=.true., reset=integrationReset)
!   end do
!   call cpu_time(cpu_time_finish)
!   print *, "integrate f(x)=tanh(x)/cosh²(x/3.)       from x=0…∞"
!   expected = 2.34839248149319
!   rel_err = abs((integral-expected)/expected)
!   print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
!   print *, '  time: ',cpu_time_finish-cpu_time_start, '(100,000 repetitions)'
!   print *,
!
!
!   integrationReset = .true.
!   call cpu_time(cpu_time_start)
!   do i = 1, 100000
!     integral=Integrate(1d-10,100.0d0,test_Pstars, integrandFunction,  &
!                        integrationWorkspace, toleranceRelative=1.0d-7, &
!                        hasSingularities=.true.,toInfinity=.false.,     &
!                        reset=integrationReset)
!   end do
!   call cpu_time(cpu_time_finish)
!   print *, "integrate f(x)=tanh(x)/cosh²(x/3.)       from x=0…100"
!   expected = 2.34839248149319
!   rel_err = abs((integral-expected)/expected)
!   print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
!   print *, '  time: ',cpu_time_finish-cpu_time_start, '(100,000 repetitions)'
!   print *,
!
!
!   integral=Integrate(1d-10,100.0d0,test_Pdm_alt, integrandFunction,  &
!                      integrationWorkspace, toleranceRelative=1.0d-7, &
!                      toInfinity=.true., reset=integrationReset)
!   print *, "integrate f(x)=exp(-x)/(x/30)/(1+x/30)²       from x=0…∞"
!   expected = 97.946378405461
!   rel_err = abs((integral-expected)/expected)
!   print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
!   print *,
!
!
!   integral=Integrate(1d-10,-1.0d0,test_Pstars_alt, integrandFunction,   &
!                       integrationWorkspace, toleranceRelative=1.0d-7, &
!                       toInfinity=.true., reset=integrationReset)
!   print *, "integrate f(x)=exp(-x)/cosh²(x/3.)       from x=0…∞"
!   expected = 2.34839248149319
!   rel_err = abs((integral-expected)/expected)
!   print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
!   print *,
!


end program Test_Integration
