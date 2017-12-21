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
! Contains a module which implements calculations of Bessel functions.
! (currently this is a wrapper to use the bessel functions from FGSL)
module Bessel_Functions
  use FGSL
  implicit none
  private
  public :: K0, K1, K2, I0, I1, I2

contains

  double precision function K0(argument)
    ! Computes the $K_0$ Bessel function.
    implicit none
    double precision, intent(in) :: argument

    K0=FGSL_SF_Bessel_Kc0(argument)
    return
  end function K0

  double precision function K1(argument)
    ! Computes the $K_1$ Bessel function.
    implicit none
    double precision, intent(in) :: argument

    K1=FGSL_SF_Bessel_Kc1(argument)
    return
  end function K1

  double precision function K2(argument)
    ! Computes the $K_1$ Bessel function.
    implicit none
    double precision, intent(in) :: argument

    K2=FGSL_SF_Bessel_Kcn(2,argument)
    return
  end function K2

  double precision function I0(argument)
    ! Computes the $I_0$ Bessel function.
    implicit none
    double precision, intent(in) :: argument

    I0=FGSL_SF_Bessel_Ic0(argument)
    return
  end function I0

  double precision function I1(argument)
    ! Computes the $I_1$ Bessel function.
    implicit none
    double precision, intent(in) :: argument

    I1=FGSL_SF_Bessel_Ic1(argument)
    return
  end function I1

  double precision function I2(argument)
    ! Computes the $I_1$ Bessel function.
    implicit none
    double precision, intent(in) :: argument

    I2=FGSL_SF_Bessel_Icn(2,argument)
    return
  end function I2

end module Bessel_Functions
