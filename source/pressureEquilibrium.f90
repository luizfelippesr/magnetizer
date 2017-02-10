! Contains subroutines which compute the pressure and density at the midplane
! and disk scaleheight assuming (magneto)hydrostatic equilibrium and accounting
! for the various pressure contributions.
module pressureEquilibrium
  use math_constants
  use fgsl
  use, intrinsic :: iso_c_binding

  implicit none
  private

  public solve_hydrostatic_equilibrium_cubic
  public solve_hydrostatic_equilibrium_numerical
  public computes_midplane_ISM_pressure_using_scaleheight
  public computes_midplane_ISM_pressure_from_B_and_rho
  public computes_midplane_ISM_pressure_using_rotation

  double precision, parameter :: G_SI=FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT
  double precision, parameter :: Msun_SI = FGSL_CONST_MKSA_SOLAR_MASS
  double precision, parameter :: kpc_SI = FGSL_CONST_MKSA_PARSEC*1d3
  double precision, parameter :: density_SI_to_cgs = 1d-3
  double precision, parameter :: convertPressureSItoGaussian=10

contains
  subroutine solve_hydrostatic_equilibrium_numerical(rdisk, M_g, M_star, r, B, &
                                                      Om, G, &
                                                      rho_d, h_d, &
                                                      Rm_out, &
                                                      Sigma_star_out, &
                                                      Sigma_d_out)
    ! Computes the mid-plane density and scale height under the assumption of
    ! hydrostatic equilibrium
    ! Input:  rdisk -> half mass radius of the disk (kpc)
    !         M_g -> total gas mass (Msun)
    !         M_star -> total stellar mass (Msun)
    !         r -> radii array (kpc)
    !         B -> magnetic field array (microgauss)
    ! Output: rho_d -> density profile of the gas at the midplane (g/cm^3)
    !         h_d   -> scaleheight profile of the gas at the midplane (kpc)
    ! Optional outputs (all as function of r)
    !         Sigma_star_out -> The surface density profile of stars (n-array)
    !         Sigma_d_out -> The surface density profile of diffuse gas (n-array)
    !         Rm_out -> Molecular to atomic gas ratio profile (n-array)
    !
    use input_constants
    use global_input_parameters
    use root_finder, only: FindRoot
    use messages, only: error_message
    double precision, intent(in) :: rdisk, M_g, M_star
    double precision, dimension(:), intent(in) :: r, B, Om, G
    double precision, dimension(size(r)), intent(out) :: rho_d, h_d
    double precision, dimension(size(r)), intent(out), optional ::    &
                                          Rm_out, Sigma_d_out, Sigma_star_out
    double precision, dimension(size(r)) :: Sigma_d, Sigma_star, Sigma_g
    double precision, dimension(size(r)) :: Rm
    double precision :: rs, rs_g
    !
    real(fgsl_double) :: minimum_h, maximum_h
    real(fgsl_double) :: pressure_equation_min, pressure_equation_max
    real(fgsl_double), target, dimension(7) :: parameters_array
    integer :: i, i_rreg
    type(c_ptr) :: params_ptr

    ! Prepares constants
    rs = constDiskScaleToHalfMassRatio * rdisk ! kpc
    rs_g = p_gasScaleRadiusToStellarScaleRadius_ratio * rs ! kpc

    ! Computes (total) gas and stars surface densities
    Sigma_g = exp_surface_density(rs_g, abs(r), M_g) ! Msun/kpc^2
    Sigma_star = exp_surface_density(rs, abs(r), M_star) ! Msun/kpc^2

    ! Computes R_mol
    Rm = molecular_to_diffuse_ratio(rdisk, Sigma_g, Sigma_star)

    ! Computes diffuse gas surface density
    Sigma_d = Sigma_g / (Rm+1d0)

    ! Prepares a pointer to parameters_array
    params_ptr = c_loc(parameters_array)

    ! Sets up a guessed initial search interval for r=0
    minimum_h = 1d-3*rdisk
    maximum_h = 1d0*rdisk
    parameters_array(1) = rdisk
    i_rreg = 2

    do i=1, size(r)
      ! If truncation is requested, store the index and skip
      if (p_truncates_within_rreg .and. r(i) < p_rreg_to_rdisk*rdisk) then
        ! This avoid being bitten by the not-very physical behaviour at r=0
        i_rreg = i
        cycle
      endif

      ! A bit of FGSL interfacing gymnastics
      parameters_array(2) = Sigma_d(i)
      parameters_array(3) = Sigma_star(i)
      parameters_array(4) = Rm(i)
      if (p_enable_P2) then
        parameters_array(5) = Om(i)
        parameters_array(6) = G(i)
      else
        parameters_array(5) = 0d0
        parameters_array(6) = 0d0
      endif
      parameters_array(7) = B(i)

      ! Initially tests whether the interval contains a root
      ! (this assumes that is an odd number of root in the interval,
      ! hopefully only one root!)
      pressure_equation_min = pressure_equation(minimum_h, params_ptr)
      pressure_equation_max = pressure_equation(maximum_h, params_ptr)
      if (pressure_equation_min/pressure_equation_max > 0d0) then
        minimum_h = 1d-5*rdisk
        maximum_h = 1d0*rdisk
        h_d(i) = 0d0
        call error_message('solve_hydrostatic_equilibrium_numerical',         &
                           'Initial guess does NOT include a change of sign.' &
                           // 'Unable to find the root (i.e. h and rho).',    &
                           code='P')
        cycle
      endif

      ! Finds the root
      h_d(i) = FindRoot(pressure_equation, params_ptr, &
                        [minimum_h, maximum_h])

      ! Uses previous h value to set up search interval for the next one
      ! (making things much faster!)
      maximum_h = h_d(i)*1.35d0
      minimum_h = h_d(i)*0.65d0
    end do

    ! If truncation was requested, uses the same value as rreg for r<rreg
    if (i_rreg>1) then
      h_d(1:i_rreg) = h_d(i_rreg+1)
    endif

    ! Outputs optional quantities
    if (present(Rm_out)) Rm_out=Rm
    if (present(Sigma_star_out)) Sigma_star_out=Sigma_star
    if (present(Sigma_d_out)) Sigma_d_out=Sigma_d

    ! Computes density
    where(h_d/=0d0)
      rho_d = Sigma_d*Msun_SI/kpc_SI/kpc_SI /h_d/2d0/kpc_SI  * density_SI_to_cgs
    elsewhere
      rho_d = 0d0
    endwhere
  end subroutine solve_hydrostatic_equilibrium_numerical


  function pressure_equation(h, params) bind(c)
    ! Hydrostatic pressure equation for interfacing with the
    ! root finder GSL routine
    real(c_double), value :: h
    type(c_ptr), value :: params
    real(c_double) :: pressure_equation
    real(fgsl_double), dimension(1) :: Sigma_d, Sigma_star, Rm
    real(fgsl_double), dimension(1) :: B, Om, G, h_d, rho_cgs
    real(fgsl_double) :: r_disk
    real(fgsl_double), dimension(1) :: P1, P2, Pgas
    real(fgsl_double), pointer, dimension(:) :: p

    call c_f_pointer(params, p, [7])

    ! Converts argument and array of parameters into small "0-arrays"
    ! (this allows using the previous routines...)
    r_disk     = p(1)
    Sigma_d    = [p(2)]
    Sigma_star = [p(3)]
    Rm         = [p(4)]
    Om         = [p(5)]
    G          = [p(6)]
    B          = [p(7)]
    h_d        = [h]

    ! Computes the midplane pressure, from gravity (without radial part)
    P1 = computes_midplane_ISM_pressure_using_scaleheight(  &
                                        r_disk, Sigma_d, Sigma_star, Rm, h_d)
    ! Computes the midplane pressure, from gravity (the radial part)
    P2 = computes_midplane_ISM_pressure_using_rotation(Sigma_d, Om, G, h_d)

    ! Computes the midplane pressure, from the density
    rho_cgs = Sigma_d*Msun_SI/kpc_SI/kpc_SI /(2d0*h_d*kpc_SI) * density_SI_to_cgs
    Pgas = computes_midplane_ISM_pressure_from_B_and_rho(B, rho_cgs)

    pressure_equation = P1(1) + P2(1) - Pgas(1)

  end function pressure_equation


  subroutine solve_hydrostatic_equilibrium_cubic(rdisk, M_g, M_star, r, B,  &
                                          rho_d, h_d, &
                                          Sigma_star_out, Sigma_d_out, Rm_out,&
                                          all_roots)
    ! Computes the mid-plane density and scale height under the assumption of
    ! hydrostatic equilibrium
    ! Input:  rdisk -> half mass radius of the disk (kpc)
    !         M_g -> total gas mass (Msun)
    !         M_star -> total stellar mass (Msun)
    !         r -> radii array (kpc)
    !         B -> magnetic field array (microgauss)
    ! Output: rho_d -> density profile of the gas at the midplane (g/cm^3)
    !         h_d   -> scaleheight profile of the gas at the midplane (kpc)
    ! Optional outputs (all as function of r)
    !         Sigma_star_out -> The surface density profile of stars (n-array)
    !         Sigma_d_out -> The surface density profile of diffuse gas (n-array)
    !         Rm_out -> Molecular to atomic gas ratio profile (n-array)
    !         all_roots -> outputs all roots from the solve (nx3-array)
    !
    use input_constants
    use global_input_parameters
    use root_finder
    double precision, intent(in) :: rdisk, M_g, M_star
    double precision, dimension(:), intent(in) :: r, B
    double precision, dimension(size(r)), intent(out) :: rho_d, h_d
    double precision, dimension(size(r),3), intent(out),optional :: all_roots
    double precision, dimension(size(r)), intent(out), optional ::    &
                                          Rm_out, Sigma_d_out, Sigma_star_out
    double precision, dimension(size(r)) :: a0, a1, a2, a3
    double precision, dimension(size(r)) :: Sigma_g_nonSI, Sigma_star_nonSI
    double precision, dimension(size(r)) :: Sigma_d, Sigma_star, Sigma_d_nonSI
    double precision, dimension(size(r)) :: B2_4pi, Rm
    double precision :: rs, rs_g_nonSI, A, bm, bs, v0, rs_nonSI
    double precision, parameter :: h_guess = 0.2d0*kpc_SI
    double precision, dimension(3) :: roots
    integer :: i
    ! Prepares constants
    rs_nonSI = constDiskScaleToHalfMassRatio * rdisk ! kpc
    rs_g_nonSI = p_gasScaleRadiusToStellarScaleRadius_ratio * rs_nonSI ! kpc
    rs = rs_nonSI * kpc_SI ! m
    A = (p_ISM_kappa/3d0*(1d0+2d0*p_ISM_xi) + 1d0/p_ISM_gamma)
    
    ! Nicknames
    bm = p_molecularHeightToRadiusScale
    bs = p_stellarHeightToRadiusScale
    v0 = p_ISM_sound_speed_km_s * 1d3 ! m/s

    ! Another shorthand: (B_Gauss)^2/(4\pi)
    B2_4pi = (B*1d-6)**2/4d0/pi * 0.1d0 ! 0.1 J/m^3 = 1 erg/cm^3

    ! Computes initial surface density profiles
    if (rs_g_nonSI==0 .or. rs_nonSI==0) then
      stop 'solve_hydrostatic_equilibrium: Fatal error. Disk of negligible size'
    endif
    
    Sigma_g_nonSI = exp_surface_density(rs_g_nonSI, abs(r), M_g)
    Sigma_star_nonSI = exp_surface_density(rs_nonSI, abs(r), M_star)

    ! Computes R_mol
    Rm = molecular_to_diffuse_ratio(rdisk, Sigma_g_nonSI, Sigma_star_nonSI)

    ! Computes diffuse gas surface density
    Sigma_d_nonSI = Sigma_g_nonSI / (Rm+1d0)
    ! Adjusts units of stellar surface density (to SI!)
    Sigma_star = (Sigma_star_nonSI*Msun_SI/kpc_SI/kpc_SI)
    ! Adjusts units of diffuse gas surface density (to SI!)
    Sigma_d = (Sigma_d_nonSI*Msun_SI/kpc_SI/kpc_SI)

    ! Computes the coefficients of the cubic equation
    a3 = pi*G_SI*Sigma_d*(Sigma_star + Sigma_d*(Rm + 0.5)) - B2_4pi

    a2 =  -B2_4pi*(bm+bs)*rs - A/2d0*Sigma_d*v0**2 &
        + pi*G_SI*Sigma_d*(Sigma_star*bm + Sigma_d*(bs*Rm + 0.5*bm + 0.5*bs))*rs

    a1 = - (B2_4pi - pi/2d0*G_SI*Sigma_d**2)*rs**2*bm*bs &
        - A/2d0 * Sigma_d * v0**2 * rs * (bm+bs)

    a0 = -A/2d0 * Sigma_d * v0**2 * rs**2 * bm * bs

    ! Solves the cubic equation
    ! (In the case of multiple solutions, considers only those with the
    ! correct order of magnitude, i.e. close to h_guess=200pc. This is a quick
    ! and dirty way of avoiding the unphysical solutions without having to
    ! substitute back in the equation.)
    do i=1,size(r)

      if (Sigma_d_nonSI(i) < 0.01) then
        rho_d(i) = 1d-28
        h_d(i) = Sigma_d(i)/rho_d(i)/2d0 / density_SI_to_cgs
        print *, 'r', r(i), h_d(i), i, size(r)
      else
        h_d(i) = CubicRootClose(a3(i), a2(i), a1(i), a0(i), h_guess, roots)
      endif

      if (h_d(i)>1d-10) then
        rho_d(i) = Sigma_d(i)/h_d(i)/2d0 * density_SI_to_cgs
      else
        rho_d(i) = 1d-100
      end if
      h_d(i) = h_d(i)/kpc_SI
      ! Outputs all roots if required
      if (present(all_roots)) all_roots(i,:) = roots(:)/kpc_SI
    end do
    ! Outputs optional quantities
    if (present(Rm_out)) Rm_out=Rm
    if (present(Sigma_star_out)) Sigma_star_out=Sigma_star_nonSI
    if (present(Sigma_d_out)) Sigma_d_out=Sigma_d_nonSI

    return
  end subroutine solve_hydrostatic_equilibrium_cubic

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
                                rdisk, Sigma_d, Sigma_star, R_m, h_d) result(P)
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
    double precision, dimension(:), intent(in) :: Sigma_d, Sigma_star
    double precision, dimension(:), intent(in) :: h_d, R_m
    double precision, dimension(size(Sigma_star)) :: weight_d, weight_star
    double precision, dimension(size(Sigma_star)) :: P
    double precision, dimension(size(Sigma_star)) :: Sigma_d_SI, Sigma_star_SI
    double precision, parameter :: rs_to_r50=constDiskScaleToHalfMassRatio
    double precision :: h_star, h_m

    ! Computes missing scale heights (in kpc)
    h_star = rdisk * rs_to_r50 * p_stellarHeightToRadiusScale
    h_m = rdisk * rs_to_r50 * p_molecularHeightToRadiusScale

    ! Adjusts units of surface densities
    Sigma_star_SI = (Sigma_star*Msun_SI/kpc_SI/kpc_SI)
    Sigma_d_SI = (Sigma_d*Msun_SI/kpc_SI/kpc_SI)

    ! Computes weighting associated with stars
    weight_star = 2d0*h_d/(h_star+h_d)

    ! Computes weighting associated with diffuse gas
    weight_d = R_m*2d0*h_d/(h_m+h_d) + 1d0

    ! Finishes calculation
    P = pi/2d0 * G_SI * Sigma_d_SI * (Sigma_d_SI*weight_d   &
                    + Sigma_star_SI*weight_star) * convertPressureSItoGaussian

    return
  end function computes_midplane_ISM_pressure_using_scaleheight


  function computes_midplane_ISM_pressure_from_B_and_rho(B, rho) result(P)
    ! Computes the ISM pressure, knowing the magnetic field and density
    ! (This helper function can be used for checking up the solution of
    !  other functions in this module)
    ! Input:  B   -> magnetic field (in microgauss)
    !         rho -> density of diffuse gas (in g/cm^3)
    ! Output: P -> midplane pressure (in erg/cm^3)
    use global_input_parameters
    double precision, dimension(:), intent(in) :: B, rho
    double precision, dimension(size(rho)) :: P

    P = (B*1d-6)**2/4d0/pi + (p_ISM_kappa/3d0*(1.+2.*p_ISM_xi) + 1d0/p_ISM_gamma) &
                             * (p_ISM_sound_speed_km_s*1d5)**2 * rho
    return

  end function computes_midplane_ISM_pressure_from_B_and_rho


  function computes_midplane_ISM_pressure_using_rotation(Sigma_d, Om, S, h_d) result(P)
    double precision, dimension(:), intent(in) :: Sigma_d, Om, S, h_d
    double precision, dimension(size(h_d)) :: P
    double precision, dimension(size(h_d)) :: Sigma_d_SI, Om_SI, S_SI, h_d_SI
    double precision, parameter :: km_SI = 1d3
    ! Computes the non-local contribution to the ISM pressure at
    ! the midplane from the rotation curve
    !
    ! Input: Sigma_d -> Surface density profile of diffuse gas (Msun/kpc^2)
    !        Om -> Angular velocity in km/s/kpc
    !        S -> Shear rate velocity in km/s/kpc
    !        h_d -> Scale-height profile of the diffuse gas (kpc)
    ! Output: array containing the pressure in Gaussian units

    Sigma_d_SI = (Sigma_d*Msun_SI/kpc_SI/kpc_SI)
    Om_SI = Om * km_SI/kpc_SI
    S_SI  = S  * km_SI/kpc_SI
    h_d_SI = h_d * kpc_SI

    P = - Om_SI*(Om_SI + S_SI)*h_d_SI*Sigma_d_SI*convertPressureSItoGaussian
  end function

end module pressureEquilibrium
