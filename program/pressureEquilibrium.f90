! Contains subroutines which compute the pressure and density at the midplane
! and disk scaleheight assuming (magneto)hydrostatic equilibrium and accounting
! for the various pressure contributions.
module pressureEquilibrium
  use math_constants
  implicit none

contains
  subroutine solves_hytrostatic_equilibrium(rdisk, M_g, M_star, r, B, rho_d, h_d)
    ! Computes the mid-plane density and scale height under the assumption of
    ! hydrostatic equilibrium
    ! Input:  rdisk -> half mass radius of the disk (kpc)
    !         M_g -> total gas mass (Msun)
    !         M_star -> total stellar mass (Msun)
    !         r -> radii array (kpc)
    !         B -> magnetic field array (microgauss)
    ! Output: rho_d -> density profile of the gas at the midplane (g/cm^3)
    !         h_d   -> scaleheight profile of the gas at the midplane (kpc)
    use input_constants
    use global_input_parameters
    use root_finder
    use fgsl
    double precision, intent(in) :: rdisk, M_g, M_star
    double precision, dimension(:), intent(in) :: r, B
    double precision, dimension(size(r)), intent(out) :: rho_d, h_d
    double precision, dimension(size(r)) :: a0, a1, a2, a3
    double precision, dimension(size(r)) :: Sigma_g_nonSI, Sigma_star_nonSI
    double precision, dimension(size(r)) :: Sigma_d, Sigma_star, Sigma_d_nonSI
    double precision, dimension(size(r)) :: B2_4pi, Rm
    double precision :: rs, A, bm, bs, v0
    double precision, parameter :: G=FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT
    double precision, parameter :: Msun_SI = FGSL_CONST_MKSA_SOLAR_MASS
    double precision, parameter :: kpc_SI = FGSL_CONST_MKSA_PARSEC*1d3
    double precision, parameter :: h_guess = 0.2d0*kpc_SI
    double precision, parameter :: density_SI_to_cgs = 1d-3
    integer :: i
    ! Prepares constants
    rs = constDiskScaleToHalfMassRatio*rdisk
    A = (p_ISM_csi + p_ISM_kappa/3d0 + 1d0/p_ISM_gamma)
    ! Nicknames
    bm = p_molecularHeightToRadiusScale
    bs = p_stellarHeightToRadiusScale
    v0 = p_ISM_csi
    ! Another shorthand: B^2/(4\pi)
    B2_4pi = B**2/4d0/pi
    ! Computes initial surface density profiles
    Sigma_g_nonSI = exp_surface_density(rs, r, M_g)
    Sigma_star_nonSI = exp_surface_density(rs, r, M_star)
    ! Computes R_mol
    Rm = molecular_to_diffuse_ratio(rdisk, Sigma_g_nonSI, Sigma_star_nonSI)
    ! Adjusts units of stellar surface density (to SI!)
    Sigma_star = (Sigma_star_nonSI*Msun_SI/kpc_SI/kpc_SI)
    ! Computes diffuse gas surface density (to SI!)
    Sigma_d_nonSI = Sigma_g_nonSI / (Rm+1d0)
    ! Adjusts units of diffuse gas surface density (to SI!)
    Sigma_d = (Sigma_d_nonSI*Msun_SI/kpc_SI/kpc_SI)

    ! Computes the coefficients of the cubic equation
    a0 = A/2d0 * Sigma_d * v0**2 * rs**2 * bm * bs

    a1 = (B2_4pi - pi/2d0*G*Sigma_d**2)*rs**2*bm*bs &
        + A/2d0 * Sigma_d * v0**2 * rs * (bm+bs)

    a2 = B2_4pi*(bm+bs)*rs + A/2d0*Sigma_d*v0**2 &
        - pi*G*Sigma_d*(Sigma_star*bm + Sigma_d*(bs*Rm + 0.5*bm + 0.5*bs))*rs

    a3 = B2_4pi - pi*G*Sigma_d*(Sigma_star + Sigma_d*(Rm + 0.5))

    do i=1,size(r)
      h_d(i) = CubicRootClose(a3(i), a2(i), a1(i), a0(i), h_guess)
      rho_d(i) = Sigma_d(i)/h_d(i)/2d0 * density_SI_to_cgs
      h_d(i) = h_d(i)/kpc_SI
    end do
    return
  end subroutine

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


  function molecular_to_diffuse_ratio(rdisk, Sigma_g, Sigma_star) result(Rmol)
    ! Computes the ratio of molecular gas to diffuse gas as a function of
    ! radius using the empirical relation found by Blitz & Rosolowsky (2006)
    ! Input: Sigma_g -> Surface density profile of gas (Msun/kpc^2)
    !        Sigma_stars -> Surface density profile of stars (Msun/kpc^2)
    !        rdisk -> half mass radius of the disk, in kpc
    ! Output: Array $f_\text{mol}$
    use global_input_parameters
    double precision, dimension(:), intent(in) :: Sigma_g, Sigma_star
    double precision, intent(in) :: rdisk
    double precision, dimension(size(Sigma_g)) :: fmol, Rmol, Pnot

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
    use global_input_parameters
    use input_constants
    use fgsl
    double precision, intent(in) :: rdisk
    double precision, optional, intent(in) :: v_gas
    double precision, dimension(:), intent(in) :: Sigma_g, Sigma_star
    double precision, dimension(size(Sigma_star)) :: P, v_star_SI
    double precision, dimension(size(Sigma_star)) :: Sigma_star_SI, Sigma_g_SI
    double precision :: h_star_SI, v_gas_SI
    double precision, parameter :: rs_to_r50=constDiskScaleToHalfMassRatio
    double precision, parameter :: convertPressureSItoGaussian=10
    double precision, parameter :: kpc_SI = FGSL_CONST_MKSA_PARSEC*1d3
    double precision, parameter :: G_SI=FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT
    double precision, parameter :: Msun_SI = FGSL_CONST_MKSA_SOLAR_MASS

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


  function computes_midplane_ISM_pressure_using_scaleheight(   &
                                    rdisk, Sigma_g, Sigma_star, h_d) result(P)
    ! Computes the pressure in the midplane using Sigma_g, Sigma_star and h_d
    ! (This helper function can be used for checking up the solution of
    !  other functions in this module)
    ! Input: Sigma_g -> Surface density profile of (total) gas (Msun/kpc^2)
    !        Sigma_stars -> Surface density profile of stars (Msun/kpc^2)
    !        h_d -> Scale-height profile of the diffuse gas (kpc)
    !        rdisk -> half mass radius of the disk (kpc)
    ! Output: array containing the pressure in Gaussian units
    use global_input_parameters
    use input_constants
    use fgsl
    double precision, intent(in) :: rdisk
    double precision, dimension(:), intent(in) :: Sigma_g, Sigma_star, h_d
    double precision, dimension(size(Sigma_star)) :: weight_d, weight_star
    double precision, dimension(size(Sigma_star)) :: P, R_m
    double precision, dimension(size(Sigma_star)) :: Sigma_d_SI, Sigma_star_SI
    double precision, parameter :: rs_to_r50=constDiskScaleToHalfMassRatio
    double precision, parameter :: convertPressureSItoGaussian=10
    double precision, parameter :: kpc_SI = FGSL_CONST_MKSA_PARSEC*1d3
    double precision, parameter :: G_SI=FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT
    double precision, parameter :: Msun_SI = FGSL_CONST_MKSA_SOLAR_MASS
    double precision :: h_star, h_m

    ! Computes missing scale heights (in kpc)
    h_star = rdisk * rs_to_r50 * p_stellarHeightToRadiusScale
    h_m = rdisk * rs_to_r50 * p_molecularHeightToRadiusScale
    ! Computes R_mol
    R_m = molecular_to_diffuse_ratio(rdisk, Sigma_g, Sigma_star)
    ! Adjusts units of stellar surface density
    Sigma_star_SI = (Sigma_star*Msun_SI/kpc_SI/kpc_SI)
    ! Computes diffuse gas surface density
    Sigma_d_SI = (Sigma_g*Msun_SI/kpc_SI/kpc_SI) / (R_m+1d0)
    ! Computes weighting associated with stars
    weight_d = 2d0*h_d/(h_star+h_d)
    ! Computes weighting associated with diffuse gas
    weight_d = R_m*2d0*h_d/(h_m+h_d) + 1d0
    ! Finishes calculation
    P = pi/2d0 * G_SI * Sigma_d_SI * (Sigma_d_SI*weight_d         &
                    + Sigma_star_SI*weight_star) * convertPressureSItoGaussian
    return
  end function computes_midplane_ISM_pressure_using_scaleheight

  function computes_midplane_ISM_pressure_from_B_and_rho(B, rho) result(P)
    ! Computes the ISM pressure, knowing the magnetic field and density
    ! (This helper function can be used for checking up the solution of
    !  other functions in this module)
    ! Input:  B   -> magnetic field (in microgauss)
    !         rho -> density (in g/cm^3)
    ! Output: P -> midplane pressure (in erg/cm^3)
    use global_input_parameters
    double precision, dimension(:), intent(in) :: B, rho
    double precision, dimension(size(rho)) :: P

    P = (B*1d6)**2/4d0/pi + (p_ISM_csi + p_ISM_kappa/3d0 + 1d0/p_ISM_gamma) &
                            * p_ISM_sound_speed_km_s**2 * rho
    return
  end function computes_midplane_ISM_pressure_from_B_and_rho

end module pressureEquilibrium
