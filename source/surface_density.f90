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
module surface_density
  ! Functions related to computing surface densities and molecular fractions
  use math_constants
  use global_input_parameters
  implicit none

contains

  function exp_surface_density(rs, r, M) result(Sigma)
    ! Construncts exponential surface density profile
    ! Input:  r -> radii array
    !         rs -> scale radius of the disk (same units as r)
    !         M -> Mass
    ! Output: Array containing Sigma (units: [M]/[r^2])
    double precision, dimension(:), intent(in) :: r
    double precision, intent(in) :: rs, M
    double precision, dimension(size(r)) :: Sigma

    Sigma = M/2./pi/rs**2 * exp(-r/rs)
    return
  end function exp_surface_density


  function molecular_gas_surface_density(r, r_disk, Mgas_disk, Mstars_disk) result(Sigma_m)
    ! Construncts surface density profile of the molecular gas
    ! Input:  r -> radii array
    !         rs -> scale radius of the disk (same units as r)
    !         M -> Mass
    ! Output: Array containing Sigma_m (units: [M]/[r^2])
    use input_constants
    double precision, dimension(:), intent(in) :: r
    double precision :: r_disk, Mgas_disk, Mstars_disk, rs, rs_g
    double precision, dimension(size(r)) :: Rm, Sigma_m, Sigma_star, Sigma_g

    ! Prepares constants
    rs = constDiskScaleToHalfMassRatio * r_disk ! kpc
    rs_g = p_gasScaleRadiusToStellarScaleRadius_ratio * rs ! kpc

    ! Computes (total) gas and stars surface densities
    Sigma_g = exp_surface_density(rs_g, abs(r), Mgas_disk) ! Msun/kpc^2
    Sigma_star = exp_surface_density(rs, abs(r), Mstars_disk) ! Msun/kpc^2

    ! Computes R_mol
    Rm = molecular_to_diffuse_ratio(r_disk, Sigma_g, Sigma_star)

    ! Computes diffuse gas surface density
    Sigma_m = Sigma_g * Rm/(Rm+1d0)
    return
  end function molecular_gas_surface_density


  function molecular_to_diffuse_ratio(rdisk, Sigma_g, Sigma_star) result(Rmol)
    ! Computes the ratio of molecular gas to diffuse gas as a function of
    ! radius using the empirical relation found by Blitz & Rosolowsky (2006)
    ! Input: Sigma_g -> Surface density profile of gas (Msun/kpc^2)
    !        Sigma_stars -> Surface density profile of stars (Msun/kpc^2)
    !        rdisk -> half mass radius of the disk, in kpc
    ! Output: Array $f_\text{mol}$
    double precision, dimension(:), intent(in) :: Sigma_g, Sigma_star
    double precision, intent(in) :: rdisk
    double precision, dimension(size(Sigma_g)) :: Rmol, Pnot

    Pnot = midplane_pressure_Elmegreen(rdisk, Sigma_g, Sigma_star)
    Rmol = (Pnot/p_Rmol_P0)**(p_Rmol_alpha)
    return
  end function molecular_to_diffuse_ratio


  function midplane_pressure_Elmegreen(rdisk, Sigma_g, Sigma_star, v_gas) result(P)
    ! Computes the pressure following the prescription of Elmegreen (1989)
    ! Input: Sigma_g -> Surface density profile of gas (Msun/kpc^2)
    !        Sigma_stars -> Surface density profile of stars (Msun/kpc^2)
    !        rdisk -> half mass radius of the disk, in kpc
    !        v_gas -> optional: velocity dispersion of gas in km/s.
    !                 Default: 10km/s
    ! Output: array containing the pressure in Gaussian units
    use input_constants
    use fgsl
    double precision, intent(in) :: rdisk
    double precision, optional, intent(in) :: v_gas
    double precision, dimension(:), intent(in) :: Sigma_g, Sigma_star
    double precision, dimension(size(Sigma_star)) :: P, v_star_SI
    double precision, dimension(size(Sigma_star)) :: Sigma_star_SI, Sigma_g_SI
    double precision :: h_star_SI, v_gas_SI
    double precision, parameter :: rs_to_r50=constDiskScaleToHalfMassRatio

    ! Computes the stellar scale-height
    h_star_SI = p_stellarHeightToRadiusScale*rs_to_r50*rdisk*kpc_SI
    ! Sets the velocity dispersion of the gas (if absent, uses sound speed)
    if (present(v_gas)) then
      v_gas_SI = v_gas * 1d3
    else
      v_gas_SI = p_ISM_sound_speed_km_s * 1d3
    endif
    ! Adjusts units of stellar and gas surface densities
    Sigma_g_SI = Sigma_g * Msun_SI / kpc_SI**2
    Sigma_star_SI = Sigma_star * Msun_SI / kpc_SI**2

    ! Computes the velocity dispersion of the stars using:
    ! $\v_\star = \sqrt{ \pi G h_\star \Sigma }$
    v_star_SI = sqrt( pi * G_SI * h_star_SI * Sigma_star_SI )
    ! Do not allow (v_gas_SI/v_star_SI > 1)
    where (v_star_SI < v_gas_SI)
      v_star_SI = v_gas_SI
    endwhere
    ! Finishes calculation
    P = pi/2d0 * G_SI * Sigma_g_SI  &
       * (Sigma_g_SI + (v_gas_SI/v_star_SI)*Sigma_star_SI) &
       * convertPressureSItoGaussian
    return
  end function midplane_pressure_Elmegreen

end module surface_density
