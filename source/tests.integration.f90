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
  type   (fgsl_function             ), save, private :: integrandFunction
  type   (fgsl_integration_workspace), save, private :: integrationWorkspace

contains
  double precision function Integrand1(x)
    !% Integral for unit testing.
    implicit none
    double precision, intent(in   ) :: x

    Integrand1=x
    return
  end function Integrand1

  double precision function Integrand2(x)
    !% Integral for unit testing.
    implicit none
    double precision, intent(in   ) :: x

    Integrand2=sin(x)
    return
  end function Integrand2

  double precision function Integrand3(x)
    !% Integral for unit testing.
    implicit none
    double precision, intent(in   ) :: x

    Integrand3=1.0d0/sqrt(x)
    return
  end function Integrand3

  double precision function Integrand4(x)
    !% Integral for unit testing.
    use Integration
    implicit none
    double precision, intent(in   ) :: x

    Integrand4=cos(x)*Integrate(0.0d0,x,Integrand1,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6)
    return
  end function Integrand4

  double precision function test_Pstars(x)
    !% Integral for unit testing.
    use Integration
    implicit none
    double precision, intent(in   ) :: x
    double precision, parameter :: xmax = 600
    double precision, parameter :: hstars_to_hd = 3.0
    if (x<xmax) then
      test_Pstars= tanh(x)/cosh(x/hstars_to_hd)**2
    else
      test_Pstars= 0d0
    endif
    return
  end function test_Pstars

  double precision function test_Pdm(x)
    !% Integral for unit testing.
    use Integration
    implicit none
    double precision, intent(in   ) :: x
    double precision, parameter :: xmax = 10000
    double precision, parameter :: hstars_to_hd = 30.0

    if (x<xmax) then
      test_Pdm = tanh(x)/(x/hstars_to_hd)/(1+(x/hstars_to_hd))**2
    else
      test_Pdm = 0d0
    endif
    return
  end function test_Pdm

  double precision function test_Pstars_alt(x)
    !% Integral for unit testing.
    use Integration
    implicit none
    double precision, intent(in   ) :: x
    double precision, parameter :: xmax = 600
    double precision, parameter :: hstars_to_hd = 3.0
    if (x<xmax) then
      test_Pstars_alt= exp(-x)/cosh(x/hstars_to_hd)**2
    else
      test_Pstars_alt= 0d0
    endif
    return
  end function test_Pstars_alt

  double precision function test_Pdm_alt(x)
    !% Integral for unit testing.
    use Integration
    implicit none
    double precision, intent(in   ) :: x
    double precision, parameter :: xmax = 10000
    double precision, parameter :: hstars_to_hd = 30.0

    if (x<xmax) then
      test_Pdm_alt = exp(-x)/(x/hstars_to_hd)/(1+(x/hstars_to_hd))**2
    else
      test_Pdm_alt = 0d0
    endif
    return
  end function test_Pdm_alt

end module
program Test_Integration
  !% Tests that numerical integration routines work.
  use Integration
  use Test_Integration_Functions
  implicit none
  type            (fgsl_function             ) :: integrandFunction
  type            (fgsl_integration_workspace) :: integrationWorkspace
  logical                                      :: integrationReset
  double precision                             :: integral, expected, rel_err
  integer i
  double precision :: cpu_time_start, cpu_time_finish

  integrationReset=.true.
!   ! Test simple integrations.
  integral=Integrate(0.0d0, 1.0d0, Integrand1, integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6, reset=integrationReset)
  print *, "integrate f(x)=x          from x=0……1"
  expected = 0.5d0
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *,

!   integrationReset=.false.
  integral=Integrate(0.0d0,2.0d0*m_pi,Integrand2,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6, reset=integrationReset)
  print *, "integrate f(x)=sin(x)     from x=0…2π"
  expected = 0.0d0
  print *, '  result=',integral,'  expected=', expected
  print *,

  integral=Integrate(0.0d0,10.0d0,Integrand3,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6, reset=integrationReset)
  print *, "integrate f(x)=1/√x       from x=0…10"
  expected = 2.0d0*sqrt(10.0d0)
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *,

  ! Test 2D integrations.
  integral=Integrate(0.0d0,2.0d0*m_pi,Integrand4,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6, reset=integrationReset)
  print *, "integrate f(x,y)=y·cos(x) from x=0…2π and y=0…x"
  expected = 2.0d0*m_pi
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *,


  integrationReset = .true.
  call cpu_time(cpu_time_start)
  do i = 1, 100000
    integral=Integrate(1d-20,-1.0d0,test_Pdm, integrandFunction,   &
                       integrationWorkspace, toleranceRelative=1.0d-7, &
                       toInfinity=.true., reset=integrationReset)
  end do
  call cpu_time(cpu_time_finish)
  print *, "integrate f(x)=tanh(x)/(x/30)/(1+x/30)²       from x=0…∞"
  expected = 97.946378405461
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *, '  time: ',cpu_time_finish-cpu_time_start, '(100,000 repetitions)'
  print *,

  integrationReset = .true.
  call cpu_time(cpu_time_start)
  do i = 1, 100000
    integral=Integrate(0.0d0,100.0d0,test_Pdm, integrandFunction,  &
                       integrationWorkspace, toleranceRelative=1.0d-7, &
                       hasSingularities=.true.,toInfinity=.false.,     &
                       reset=integrationReset)
  end do
  call cpu_time(cpu_time_finish)
  print *, "integrate f(x)=tanh(x)/(x/30)/(1+x/30)²       from x=0…100"
  expected = 97.946378405461
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *, '  time: ',cpu_time_finish-cpu_time_start, '(100,000 repetitions)'
  print *,

  integrationReset = .true.
  call cpu_time(cpu_time_start)
  do i = 1, 100000
    integral=Integrate(1d-10,-1.0d0,test_Pstars, integrandFunction,   &
                       integrationWorkspace, toleranceRelative=1.0d-7, &
                       toInfinity=.true., reset=integrationReset)
  end do
  call cpu_time(cpu_time_finish)
  print *, "integrate f(x)=tanh(x)/cosh²(x/3.)       from x=0…∞"
  expected = 2.34839248149319
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *, '  time: ',cpu_time_finish-cpu_time_start, '(100,000 repetitions)'
  print *,


  integrationReset = .true.
  call cpu_time(cpu_time_start)
  do i = 1, 100000
    integral=Integrate(1d-10,100.0d0,test_Pstars, integrandFunction,  &
                       integrationWorkspace, toleranceRelative=1.0d-7, &
                       hasSingularities=.true.,toInfinity=.false.,     &
                       reset=integrationReset)
  end do
  call cpu_time(cpu_time_finish)
  print *, "integrate f(x)=tanh(x)/cosh²(x/3.)       from x=0…100"
  expected = 2.34839248149319
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *, '  time: ',cpu_time_finish-cpu_time_start, '(100,000 repetitions)'
  print *,


  integral=Integrate(1d-10,100.0d0,test_Pdm_alt, integrandFunction,  &
                     integrationWorkspace, toleranceRelative=1.0d-7, &
                     toInfinity=.true., reset=integrationReset)
  print *, "integrate f(x)=exp(-x)/(x/30)/(1+x/30)²       from x=0…∞"
  expected = 97.946378405461
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *,


  integral=Integrate(1d-10,-1.0d0,test_Pstars_alt, integrandFunction,   &
                      integrationWorkspace, toleranceRelative=1.0d-7, &
                      toInfinity=.true., reset=integrationReset)
  print *, "integrate f(x)=exp(-x)/cosh²(x/3.)       from x=0…∞"
  expected = 2.34839248149319
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *,



end program Test_Integration
