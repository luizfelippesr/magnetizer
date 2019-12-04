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
module tools
  ! A set of general purpose tools
  implicit none

contains

  function linspace(z_min, z_max, n) result(array)
    double precision, intent(in) :: z_min,z_max
    double precision, allocatable, dimension(:) :: array
    double precision :: step
    integer :: n,i
    step = (z_max-z_min)/dble(n-1)
    allocate(array(n))
    array = (/( (i-1)*step+z_min, i=1 , n)/)
  end function linspace

end module tools
