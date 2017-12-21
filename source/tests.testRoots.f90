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
program testRoot
  use root_finder
  implicit none
  double precision :: x, a3, a2, a1, a0, guess

  a3 = 1d0
  a2 = -1d0
  a1 = -2d0
  a0 = -1d0
  guess = 20

  x = CubicRootClose(a3, a2, a1, a0, guess)
  print *, '(a3,a2,a1,a0)=',a3, a2, a1, a0
  print *, 'Found x=',x
  print *, 'Value at solution:', a3*x**3+a2*x**2+a1*x+a0

  a3 = 99.7d0
  a2 = 0.78787878d0
  a1 = 45d1
  a0 = -17d2
  guess = 10000

  x = CubicRootClose(a3, a2, a1, a0, guess)
  print *, '(a3,a2,a1,a0)=',a3, a2, a1, a0
  print *, 'Found x=',x
  print *, 'Value at solution:', a3*x**3+a2*x**2+a1*x+a0

end program testRoot
