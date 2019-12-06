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
module FRB
  use random
  use surface_density, only: molecular_gas_surface_density
  use tools, only: linspace
  implicit none
  private

  public get_FRB_position
contains
  function get_FRB_position(h_FRB, r_disk, Mgas_disk, Mstars_disk, rmax) result(FRB_pos)
    integer, parameter :: n = 161
    double precision, dimension(n) :: R, Sigma_m, vertical_density
    double precision, intent(in) :: h_FRB, r_disk, Mgas_disk, Mstars_disk, rmax
    double precision :: R_FRB, x_FRB, y_FRB, z_FRB, tmp
    double precision, dimension(3) :: FRB_pos

    ! Evaluates the molecular surface density in the interval R=0..rmax
    R = linspace(0d0,rmax,n)
    Sigma_m = molecular_gas_surface_density(R, r_disk, Mgas_disk, Mstars_disk)
    ! Gets a random radius with probability proportional to the surface density
    R_FRB = draw_from_pdf(R, Sigma_m)

    ! Knowing R, gets x
    call random_number(x_FRB) ! uniform random between 0..1
    x_FRB = 2*x_FRB*R_FRB - R_FRB ! uniform random between -R_FRB and R_FRB
    ! Knowing R and x, gets y
    y_FRB = random_sign()*sqrt(R_FRB**2 - x_FRB**2) ! random y compatible with x and R
    ! This procedure leads to the correct distribution of x^2+y^2, but may
    ! introduce bias, as y and x are chosen differently.
    ! To avoid this, they are randomly swapped below
    if (random_sign()>0d0) then
      tmp = y_FRB
      y_FRB = x_FRB
      x_FRB = tmp
    endif

    ! Evaluates the vertical density in the interval R=-rmax,rmax
    R = linspace(-4.5*h_FRB, 4.5*h_FRB, n)
    vertical_density =  exp(-abs(R)/h_FRB)
    z_FRB = draw_from_pdf(R, vertical_density)

    FRB_pos(1) = x_FRB; FRB_pos(2) = y_FRB; FRB_pos(3) = z_FRB
  end function get_FRB_position
end module FRB
