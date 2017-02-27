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
  logical                            , save, private :: integrationReset    =.true.

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

  double precision function fakePressure(x)
    !% Integral for unit testing.
    use Integration
    implicit none
    double precision, intent(in   ) :: x

    fakePressure=exp(-x)
    return
  end function fakePressure

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

  ! Begin unit tests.
!   call Unit_Tests_Begin_Group("Numerical integration")

  ! Test simple integrations.
  integral=Integrate(0.0d0, 1.0d0, Integrand1, integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6)
  print *, "integrate f(x)=x          from x=0……1"
  expected = 0.5d0
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *,

  integral=Integrate(0.0d0,2.0d0*m_pi,Integrand2,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6)
  print *, "integrate f(x)=sin(x)     from x=0…2π"
  expected = 0.0d0
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected
  print *,

  integrationReset=.true.
  integral=Integrate(0.0d0,10.0d0,Integrand3,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6)
  print *, "integrate f(x)=1/√x       from x=0…10"
  expected = 2.0d0*sqrt(10.0d0)
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *,

  ! Test 2D integrations.
  integrationReset=.true.
  integral=Integrate(0.0d0,2.0d0*m_pi,Integrand4,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6)
  print *, "integrate f(x,y)=y·cos(x) from x=0…2π and y=0…x"
  expected = 2.0d0*m_pi
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *,

  integrationReset=.true.
  integral=Integrate(0.0d0,-1.0d0,fakePressure, integrandFunction&
       &,integrationWorkspace, toInfinity=.true.,toleranceRelative=1.0d-10)
  print *, "integrate f(x)=exp(-x**2)       from x=0…oo"
  expected = 1d0
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *,

    integrationReset=.true.
  integral=Integrate(0.0d0,1000.0d0,fakePressure, integrandFunction&
       &,integrationWorkspace, toInfinity=.false.,toleranceRelative=1.0d-10)
  print *, "integrate f(x)=exp(-x**2)       from x=0…100"
  expected = 1d0
  rel_err = abs((integral-expected)/expected)
  print *, '  result=',integral,'  expected=', expected, 'rel_err =', rel_err
  print *,

end program Test_Integration
