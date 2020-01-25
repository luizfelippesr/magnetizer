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

  public get_FRB_position, get_FRB_LoS_z_position
contains
  function get_FRB_position(h_FRB, r_disk, Mgas_disk, Mstars_disk, rmax) result(FRB_pos)
    ! Computes a FRB position, under the assumption that it occurs randomly
    ! with a probability proportional to the molecular density
    !
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

  function get_FRB_LoS_z_position(x, y, z, h_FRB, r_disk, Mgas_disk, &
                                  Mstars_disk, rmax) result(z_FRB)
    ! Given a line of sight, computes a FRB position, under the assumption that
    ! it occurs randomly with a probability proportional to the molecular density
    !
    double precision, intent(in) :: y
    double precision, dimension(:), intent(in) :: x, z
    double precision, intent(in) :: h_FRB, r_disk, Mgas_disk, Mstars_disk, rmax
    double precision, dimension(size(x)) :: R, Sigma_m, vertical_density
    double precision, dimension(size(x)) :: molecular_density
    double precision :: z_FRB
    double precision, parameter :: MINIMAL_MOLECULAR_DENSITY = 1d-6
    double precision, parameter :: z_FRB_INVALID_MARKER = -1d8

    ! Computes the values of the cylindrical radius, R, for the line of sight
    R = sqrt(x**2 + y**2)
    ! The probability is proportional to the molecular density
    Sigma_m = molecular_gas_surface_density(R, r_disk, Mgas_disk, Mstars_disk)
    vertical_density =  exp(-abs(z)/h_FRB)
    molecular_density = Sigma_m*vertical_density ! indicative molecular density

    ! Traps the case where there is negligible molecular gas on the sightline
    if (maxval(molecular_density) < MINIMAL_MOLECULAR_DENSITY) then
      z_FRB = z_FRB_INVALID_MARKER
      return
    endif

    ! Finds the z position of the FRB, assuming the prob. is prop. to density
    z_FRB = draw_from_pdf(z, molecular_density) ! NB draw_from_pdf renormalizes
    ! Also, x_FRB can be found *later* since it is on the same LoS
    return
  end function get_FRB_LoS_z_position
end module FRB
