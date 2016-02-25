! Contains subroutines which compute the pressure and density at the midplane
! and disk scaleheight assuming (magneto)hydrostatic equilibrium and accounting
! for the various pressure contributions.
module pressureEquilibrium
  use math_constants
  implicit none

contains
  function exp_surface_density(r, M, rdisk) result Sigma
    ! Construncts exponential surface density profile
    ! Input: r -> radii array
    !        rdisk -> half mass radius of the disk (same units as r)
    !        M -> Mass
    ! Output: Array containing Sigma (units: [M]/[r^2])
    use input_constants
    double precision, dimension(:), intent(in) :: r
    double precision, intent(in) :: rdisk, M
    double precision, dimension(size(r)) :: Sigma
    double precision :: rs

    rs = constDiskScaleToHalfMassRatio*rdisk
    Sigma = M/2./pi/rs**2 * exp(-r/rs)
    return
  end function exp_surface_density

  function molecular_fraction(Sigma_g, Sigma_star)
    ! Computes the fraction of molecular gas as function of radius
    ! using the empirical relation found by Blitz & Rosolowsky (2004,2006)
    ! Input: Sigma_g -> Surface density profile of gas (Msun/kpc^2)
    !        Sigma_stars -> Surface density profile of stars (Msun/kpc^2)
    !        beta_star -> optional, stellar scale-height to scale-length
    !                     ratio. Default: beta_star = p_beta_star

    !        M -> Mass
    ! Output: Array containing Sigma (units: [M]/[r^2])

    double precision, intent(in) :: rdisk, Mgas, Mstars
    double precision, dimension(size(r)) :: P
  end molecular_fraction


  function midplane_pressure_Elmegreen(Sigma_g, Sigma_star, rdisk, v_gas) result(P)
    ! Computes the pressure following the prescription of Elmegreen (1989)
    ! Input: Sigma_g -> Surface density profile of gas (Msun/kpc^2)
    !        Sigma_stars -> Surface density profile of stars (Msun/kpc^2)
    !        rdisk -> half mass radius of the disk, in kpc
    !        v_gas -> optional: velocity dispersion of gas in km/s.
    !                 Default: 10km/s
    ! Output: array containing the pressure in Gaussian units
    use global_input_parameters
    double precision, intent(in) :: rdisk
    double precision, optional, intent(in) :: v_gas
    double precision, dimension(:), intent(in) :: Sigma_g, Sigma_star
    double precision, dimension(size(Sigma_star)) :: P, v_star_SI
    double precision, dimension(size(Sigma_star)) :: Sigma_star_SI, Sigma_g_SI
    double precision, parameter :: rs_to_r50=constDiskScaleToHalfMassRatio
    double precision, parameter :: convertPressureSItoGaussian=10
    double precision, parameter :: kpc_SI = FGSL_CONST_MKSA_PARSEC*1d3
    double precision, parameter :: G_SI=FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT
    double precision, parameter :: Msun_SI = FGSL_CONST_MKSA_SOLAR_MASS
    double precision, parameter :: stellarHeightToRadiusScale = 1d0/7.3

    ! Some unit adjustments
    h_star_SI = stellarHeightToRadiusScale*rs_to_r50*rdisk*kpc_SI
    if (present(v_gas)) then
      v_gas_SI = v_gas * 1d3
    else:
      v_gas_SI = 10d3
    endif
    Sigma_g_SI = Sigma_g * Msun_SI / kpc_SI**2
    Sigma_star_SI = Sigma_star * Msun_SI / kpc_SI**2

    ! \v_stars = \sqrt( \pi G h \Sigma )
    v_star_SI = sqrt( pi * G_SI * h_star_SI * Sigma_star_SI )
    ! Ceiling
    where (v_star_SI < v_gas_SI)
      v_star_SI = v_gas_SI
    endwhere

    P = pi/2d0 * G_SI * Sigma_g_SI  &
       * (Sigma_g_SI + (v_gas_SI/v_star_SI)*Sigma_star_SI) &
       * convertPressureSItoGaussian
    return
  end function midplane_pressure_Elmegreen


end module pressureEquilibrium
