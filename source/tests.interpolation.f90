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
program testInterpolation
  ! Small program for testing the interpolation routines
  use interpolation
  implicit none
  double precision, dimension(30) :: x_ref,y_ref
  double precision, dimension(20) :: x_20,y_20
  double precision, dimension(40) :: x_40, y_40
  integer :: i
  character(len=100), parameter :: method ='akima'
  integer, parameter :: u= 17

  ! Produces a uniform grid from 0 to 100 (for each variable)
  x_ref = [ (100*dble(i-1)/dble(size(x_ref)-1), i=1,size(x_ref),1) ]
  x_20 = [ (100*dble(i-1)/dble(size(x_20)-1), i=1,size(x_20),1) ]
  x_40 = [ (100*dble(i-1)/dble(size(x_40)-1), i=1,size(x_40),1) ]

  ! Computes the reference values (something complicated...)
  y_ref = x_ref**2*sin(x_ref/10)/100+cos(x_ref/4.)*3

  call interpolate(x_ref, y_ref, x_40, y_40)
  call interpolate(x_ref, y_ref, x_20, y_20)

  open(unit=u, file='/tmp/interp_test_ref.dat')
  do i=1,size(x_ref)
      write(u, '(2(F9.2))') x_ref(i), y_ref(i)
  end do
  close(u)

  open(unit=u, file='/tmp/interp_test_40.dat')
  do i=1,size(x_40)
      write(u, '(2(F9.2))') x_40(i), y_40(i)
  end do
  close(u)

  open(unit=u, file='/tmp/interp_test_20.dat')
  do i=1,size(x_20)
      write(u, '(2(F9.2))') x_20(i), y_20(i)
  end do
  close(u)


  ! Produces a uniform grid from 0 to 100 (for each variable)
  x_ref = [ (100*dble(i-1)/dble(size(x_ref)-1), i=size(x_ref),1,-1) ]
  x_20 = [ (100*dble(i-1)/dble(size(x_20)-1), i=size(x_20),1,-1) ]
  x_40 = [ (100*dble(i-1)/dble(size(x_40)-1), i=size(x_40),1,-1) ]

  ! Computes the reference values (something complicated...)
  y_ref = x_ref**2*sin(x_ref/10)/100+cos(x_ref/4.)*3

  call interpolate(x_ref, y_ref, x_40, y_40)
  call interpolate(x_ref, y_ref, x_20, y_20)

  open(unit=u, file='/tmp/interp_test_ref_r.dat')
  do i=1,size(x_ref)
      write(u, '(2(F9.2))') x_ref(i), y_ref(i)
  end do
  close(u)

  open(unit=u, file='/tmp/interp_test_40_r.dat')
  do i=1,size(x_40)
      write(u, '(2(F9.2))') x_40(i), y_40(i)
  end do
  close(u)

  open(unit=u, file='/tmp/interp_test_20_r.dat')
  do i=1,size(x_20)
      write(u, '(2(F9.2))') x_20(i), y_20(i)
  end do
  close(u)

end program testInterpolation
