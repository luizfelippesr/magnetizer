! Contains subroutines which compute the pressure and density at the midplane
! and disk scaleheight assuming (magneto)hydrostatic equilibrium and accounting
! for the various pressure contributions.
module pressureEquilibrium
  use math_constants
  use fgsl
  use global_input_parameters
  use input_params, only: r_bulge , r_halo, r_disk, Mhalo, Mstars_disk, &
                          Mstars_bulge, Mgas_disk, nfw_cs1
  use grid, only: r_kpc, nxghost
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public solve_hydrostatic_equilibrium_cubic
  public solve_hydrostatic_equilibrium_numerical
  public computes_midplane_ISM_pressure_P1
  public computes_midplane_ISM_pressure_from_B_and_rho
  public computes_midplane_ISM_pressure_P2

  double precision, parameter :: G_SI=FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT
  double precision, parameter :: Msun_SI = FGSL_CONST_MKSA_SOLAR_MASS
  double precision, parameter :: kpc_SI = FGSL_CONST_MKSA_PARSEC*1d3
  double precision, parameter :: density_SI_to_cgs = 1d-3
  double precision, parameter :: convertPressureSItoGaussian=10


  ! Global variables used in the integration
  double precision :: i_R_to_Rsdm, i_hd_to_RsDM
  double precision :: i_R_to_Rb,   i_hd_to_Rb
  double precision :: i_hd_to_hm, i_hd_to_hstar
  double precision, dimension(:), pointer :: pOm_h, pG_h
  !double precision, dimension(:), pointer :: pOm_h_e, pG_h_e

  type            (fgsl_function             ) :: integrandFunction
  type            (fgsl_integration_workspace) :: integrationWorkspace
  logical                                      :: integrationReset = .true.

contains
  subroutine solve_hydrostatic_equilibrium_numerical( r, B, &
                                                      Om, G, &
                                                      Om_h, G_h, &
                                                      Om_b, G_b, &
                                                      rho_d, h_d, &
                                                      Rm_out, &
                                                      Sigma_star_out, &
                                                      Sigma_d_out, &
                                                      nghost)
    ! Computes the mid-plane density and scale height under the assumption of
    ! hydrostatic equilibrium
    ! Input:  r -> radii array (kpc)
    !         B -> magnetic field array (microgauss)
    ! Output: rho_d -> density profile of the gas at the midplane (g/cm^3)
    !         h_d   -> scaleheight profile of the gas at the midplane (kpc)
    ! Optional outputs (all as function of r)
    !         Sigma_star_out -> The surface density profile of stars (n-array)
    !         Sigma_d_out -> The surface density profile of diffuse gas (n-array)
    !         Rm_out -> Molecular to atomic gas ratio profile (n-array)
    !
    use input_constants
    use root_finder, only: FindRoot
    use messages, only: error_message
    double precision, dimension(:), intent(in) :: r, B, Om, G
    double precision, dimension(:), target, intent(in) :: Om_h, G_h
    double precision, dimension(:), target, intent(in) :: Om_b, G_b
    integer, intent(in), optional :: nghost
    double precision, dimension(size(r)), intent(out) :: rho_d, h_d
    double precision, dimension(size(r)), intent(out), optional ::    &
                                          Rm_out, Sigma_d_out, Sigma_star_out
    double precision, dimension(size(r)) :: Sigma_d, Sigma_star, Sigma_g
    double precision, dimension(size(r)) :: Rm
    double precision :: rs, rs_g
    !
    real(fgsl_double) :: minimum_h, maximum_h, guess_h, tmp
    real(fgsl_double) :: pressure_equation_min, pressure_equation_max
    real(fgsl_double), target, dimension(11) :: parameters_array, parameters_array_alt
    integer :: i, i_rreg, nghost_actual, error_count
    type(c_ptr) :: params_ptr, params_ptr_alt

    if (present(nghost)) then
      nghost_actual = nghost
    else
      nghost_actual = 0
    endif

    error_count = 0

    ! Prepares constants
    rs = constDiskScaleToHalfMassRatio * r_disk ! kpc
    rs_g = p_gasScaleRadiusToStellarScaleRadius_ratio * rs ! kpc

    ! Computes (total) gas and stars surface densities
    Sigma_g = exp_surface_density(rs_g, abs(r), Mgas_disk) ! Msun/kpc^2
    Sigma_star = exp_surface_density(rs, abs(r), Mstars_disk) ! Msun/kpc^2

    ! Computes R_mol
    Rm = molecular_to_diffuse_ratio(r_disk, Sigma_g, Sigma_star)

    ! Computes diffuse gas surface density
    Sigma_d = Sigma_g / (Rm+1d0)

    ! Prepares a pointer to parameters_array
    params_ptr = c_loc(parameters_array)
    params_ptr_alt = c_loc(parameters_array_alt)

    ! Sets up a guessed initial search interval for r=0
    minimum_h = 1d-3*r_disk
    maximum_h = 2d0*r_disk

    ! Sets global pointers to acess Omega and Shear profiles
    !pOm_h => Om_h
    !pOm_h_e => Om_h_e
    !pG_h => G_h
    !pG_h_e => G_h_e

    ! Trapping possible errors
    i_rreg = -1
    h_d = -17d0

    ! Loops through the radii, skipping ghost zone and centre if nghost
    ! was specified
    do i=nghost_actual+2, size(r)
      if (p_truncates_within_rreg .and. r(i) < p_rreg_to_rdisk*r_disk) then
        ! If truncation is requested, stores the index and skips
        i_rreg = i
        cycle
      endif

      ! A bit of FGSL interfacing gymnastics
      parameters_array(1) = r(i)
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
      parameters_array(7) = Om_h(i)
      parameters_array(8) = G_h(i)
      parameters_array(9) = Om_b(i)
      parameters_array(10) = G_b(i)
      parameters_array(11) = B(i)

      parameters_array_alt = parameters_array
      parameters_array_alt(5) = 0d0
      parameters_array_alt(6) = 0d0

      ! Initially tests whether the interval contains a root
      ! (this assumes that is an odd number of roots in the interval,
      ! hopefully only one root!)
      pressure_equation_min = pressure_equation(minimum_h, params_ptr)
      pressure_equation_max = pressure_equation(maximum_h, params_ptr)

      if (sign(1d0,pressure_equation_min) /= sign(1d0,pressure_equation_max)) then
        ! If there is a change in sign, finds the root
        h_d(i) = FindRoot(pressure_equation, params_ptr, &
                          [minimum_h, maximum_h])
      else
        ! Otherwise, tries again with a larger interval (i.e. a second chance)
        minimum_h = minimum_h/2d0
        maximum_h = maximum_h*2d0

        pressure_equation_min = pressure_equation(minimum_h, params_ptr)
        pressure_equation_max = pressure_equation(maximum_h, params_ptr)

        if (sign(1d0,pressure_equation_min) /= sign(1d0,pressure_equation_max)) then
          h_d(i) = FindRoot(pressure_equation, params_ptr, &
                        [minimum_h, maximum_h])
        else
          ! If the second chance fails, reports the error
          if (error_count<3) then
            ! If it only happened a few times, apply a rough patch
            call error_message('solve_hydrostatic_equilibrium_numerical',        &
                              'Initial guess does NOT include a change of sign.'&
                              // ' Unable to find the root (i.e. h and rho).'   &
                              // ' Will continue using previous h value.', code='P')
            h_d(i) = h_d(i-1)
            error_count = error_count + 1
          else
            ! Otherwise signals serious error!
            h_d = -1
            rho_d = 0d0
            return
          endif
        endif
      endif
      rho_d(i) = density_to_scaleheight(h_d(i), Sigma_d(i))

      ! emergency fix --------------
      !if (i> nghost_actual+2 .and. h_d(i)<0) then
      !  h_d(i) = guess_h
      !  stop
      !endif
      ! end emergency fix --------------

      ! Uses previous h value to set up search interval for the next one.
      ! The guess assumes an exponentially increasing scaleheight
      if (i/=size(r)) then
!         if (h_d(i)>0) then
            guess_h = h_d(i)*exp( (r(i+1)-r(i))/rs )
!         else
!             ! If previous calculation was an error, apply a patch
!             guess_h = h_d(i-1)*exp( (r(i+1)-r(i-1))/rs )
!         endif
        maximum_h = guess_h*1.2d0 ! assumes the root will be within \pm 20%
        minimum_h = guess_h*0.8d0 ! assumes the root will be within \pm 20%
      end if

      if (rho_d(i) < p_minimum_density) then
        ! Keeps density larger than the minimum
        rho_d(i) = p_minimum_density
        ! Ajusts scaleheight accordingly
        h_d(i) = density_to_scaleheight(rho_d(i), Sigma_d(i))
      endif
    end do

    if (nghost_actual>0) then
      do i=1, nghost_actual
        h_d(i) = h_d(2*nghost_actual+2-i)
      end do

      ! Extrapolates linearly to get the middlepoint (to avoid r=0 problems)
      h_d(nghost_actual+1) = h_d(nghost_actual+2) + (r(nghost_actual+1) -r(nghost_actual+2))&
          *(h_d(nghost_actual+2)-h_d(nghost_actual+3))/(r(nghost_actual+2)-r(nghost_actual+3))

      do i=1, nghost_actual+1
        rho_d(i) = density_to_scaleheight(h_d(i), Sigma_d(i))
      enddo
    endif

    ! If truncation was requested, uses the same value as rreg for r<rreg
    if (i_rreg>1) then
      h_d(1:i_rreg) = h_d(i_rreg+1)
      rho_d(1:i_rreg) = rho_d(i_rreg+1)
    endif

    ! Outputs optional quantities
    if (present(Rm_out)) Rm_out=Rm
    if (present(Sigma_star_out)) Sigma_star_out=Sigma_star
    if (present(Sigma_d_out)) Sigma_d_out=Sigma_d

  end subroutine solve_hydrostatic_equilibrium_numerical

  elemental function density_to_scaleheight(h_d, Sigma_d) result(rho_d)
    ! Converts a scale-height in kpc into a density in g/cm^3
    ! (or vice-versa)
    double precision, intent(in) :: h_d, Sigma_d
    double precision :: rho_d
    rho_d = Sigma_d*Msun_SI/kpc_SI/kpc_SI /h_d/2d0/kpc_SI  * density_SI_to_cgs
  end function density_to_scaleheight


  function pressure_equation(h, params) bind(c)
    ! Hydrostatic pressure equation for interfacing with the
    ! root finder GSL routine
    real(c_double), value :: h
    type(c_ptr), value :: params
    real(c_double) :: pressure_equation
    real(fgsl_double) :: r, Sigma_d, Sigma_star, R_m
    real(fgsl_double) :: B, Om, G, rho_cgs
    real(fgsl_double) :: Om_h, G_h, Om_b, G_b
    real(fgsl_double) :: P1, P2, Pgas
    real(fgsl_double), pointer, dimension(:) :: p

    call c_f_pointer(params, p, [11])

    ! Converts argument and array of parameters into small "0-arrays"
    ! (this allows using the previous routines...)
    r          = p(1)
    Sigma_d    = p(2)
    Sigma_star = p(3)
    R_m        = p(4)
    Om         = p(5)
    G          = p(6)
    Om_h       = p(7)
    G_h        = p(8)
    Om_b       = p(9)
    G_b        = p(10)
    B          = p(11)

    ! Computes the midplane pressure, from gravity (without radial part)
    P1 = computes_midplane_ISM_pressure_P1(r, Sigma_d, Sigma_star, R_m, h, &
                                            Om_h, G_h, Om_b, G_b)
    ! Computes the midplane pressure, from gravity (the radial part)
    P2 = computes_midplane_ISM_pressure_P2(Sigma_d, Om, G, h)
    if (P2>0) P2 = 0

    ! Computes the midplane pressure, from the density
    rho_cgs = density_to_scaleheight(h, Sigma_d)
    Pgas = computes_midplane_ISM_pressure_from_B_and_rho(B, rho_cgs)

    pressure_equation = P1 + P2 - Pgas

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


  function computes_midplane_ISM_pressure_P1(r, Sigma_d, Sigma_star, R_m, &
                                             h_d, Om_h, G_h, Om_b, G_b,   &
                                             Pd_out, Pm_out, Pstars_out,  &
                                             Pbulge_out, Pdm_out) result(P)
    ! Computes the pressure in the midplane using Sigma_g, Sigma_star and h_d
    ! Input: Sigma_g -> Surface density profile of (total) gas (Msun/kpc^2)
    !        Sigma_stars -> Surface density profile of stars (Msun/kpc^2)
    !        h_d -> Scale-height profile of the diffuse gas (kpc)
    !        rdisk -> half mass radius of the disk (kpc)
    ! Output: array containing the pressure in Gaussian units
    use input_constants
    use Integration
    double precision, intent(in) :: r, Sigma_d, Sigma_star, R_m, h_d
    double precision, intent(in) :: Om_h, G_h, Om_b, G_b
    double precision,intent(out),optional :: Pd_out,Pm_out,Pstars_out,Pbulge_out,Pdm_out
    double precision :: P, Pm, Pstars, Pdm, Pbulge, Pd
    double precision :: Sigma_d_SI, Sigma_star_SI
    double precision, parameter :: rs_to_r50=constDiskScaleToHalfMassRatio
    double precision :: h_star, h_m, rs_DM, rho_dm_SI, rho_b_SI, rb
    double precision :: I_m, I_stars, I_dm, I_bulge, preFactor
    double precision, parameter :: I_d = 0.5
    double precision, parameter :: rel_TOL = 1d-8
    double precision, parameter :: small = 1d-12
    double precision, parameter :: rb_to_r50 = (sqrt(2.0d0)-1.0d0)
    double precision, parameter :: Mstars_bulge_min = 1d3 ! Msun
    double precision, parameter :: km_SI = 1d3
    integer :: i, i_start, j

    ! Computes missing scale heights (in kpc)
    h_star = r_disk * rs_to_r50 * p_stellarHeightToRadiusScale
    h_m = r_disk * rs_to_r50 * p_molecularHeightToRadiusScale
    rs_DM = r_halo * nfw_cs1

    ! Computes DM and bulge reference densities
    rho_dm_SI = Mhalo*Msun_SI*(rs_DM*kpc_SI)**(-3)
    rho_dm_SI = rho_dm_SI * (log(1d0+1/nfw_cs1) - 1d0/(1d0+1))

    ! Adjusts units of surface densities
    Sigma_star_SI = (Sigma_star*Msun_SI/kpc_SI/kpc_SI)
    Sigma_d_SI = (Sigma_d*Msun_SI/kpc_SI/kpc_SI)

    ! Sets global variables (used in integrations)
    i_hd_to_hm = h_d/h_m
    i_hd_to_hstar = h_d/h_star

    preFactor = pi * G_SI * Sigma_d_SI * convertPressureSItoGaussian

    Pd = preFactor * Sigma_d_SI * 1d0 * I_d

    if (.not.p_sech2_profile) then
      ! If a exp(-z/h) profile is assumed, uses the analytic simple solution
      I_stars = h_d/(h_star+h_d)
      I_m = h_d/(h_m+h_d)
    else
      ! If a sech^2(z/h) is chosen, do the numerical integration
      I_m = i_hd_to_hm *Integrate(small, 0d0, integrand_m,   &
                                  integrandFunction,         &
                                  integrationWorkspace,      &
                                  toleranceRelative=rel_TOL, &
                                  toInfinity=.true.,         &
                                  reset=integrationReset)

      I_stars = i_hd_to_hstar * Integrate(small, 0d0, integrand_stars, &
                                          integrandFunction,           &
                                          integrationWorkspace,        &
                                          toleranceRelative=rel_TOL,   &
                                          toInfinity=.true.,           &
                                          reset=integrationReset)
    endif
    Pm = preFactor * Sigma_d_SI*R_m *  I_m
    Pstars = preFactor * Sigma_star_SI * I_stars

    ! Galaxy stellar bulge
    if (p_use_Pbulge .and. Mstars_bulge>Mstars_bulge_min) then
      if (p_test_bulge_old) then
        ! Prepares global variables
        rb = r_bulge * rb_to_r50
        rho_b_SI = Mstars_bulge*Msun_SI/(2d0*pi)*(rb*kpc_SI)**(-3)
        i_hd_to_rb = h_d/rb
        i_R_to_Rb = r/rb
        ! Numerically integrates and computes result
        I_bulge = Integrate(small, 7d0, integrand_bulge, integrandFunction,  &
                            integrationWorkspace, toleranceRelative=rel_TOL, &
                            toInfinity=.false., reset=integrationReset)
        Pbulge = preFactor * 2d0 * rho_b_SI * (h_d*kpc_SI) * I_bulge
      else
        Pbulge =  Om_b*(1.5d0*Om_b+G_b) * (km_SI/kpc_SI)**2
        Pbulge = Pbulge * Sigma_d_SI * h_d * kpc_SI * convertPressureSItoGaussian
      endif
    else
      Pbulge = 0d0
    endif

    ! Halo dark matter halo
    if (p_use_Pdm) then
      if (p_test_DM_old) then
        ! Numerically integrates and computes result
        ! NB this does NOT account for halo contraction!
        ! Prepares global variables
        i_hd_to_Rsdm = h_d/rs_DM
        i_R_to_Rsdm = r/rs_DM

        I_dm = Integrate(0d0, 100d0, integrand_dm, integrandFunction,     &
                        integrationWorkspace, toleranceRelative=rel_TOL, &
                        toInfinity=.true., reset=integrationReset)
        Pdm = preFactor * 2d0 * rho_dm_SI * (h_d*kpc_SI) * I_dm

      else
        ! NB this does account for halo contraction but involves
        ! major approximations
        Pdm =  Om_h*(1.5d0*Om_h+G_h) * (km_SI/kpc_SI)**2
        Pdm = Pdm * Sigma_d_SI * h_d * kpc_SI * convertPressureSItoGaussian
      endif
    endif
    ! Finishes calculation
    P = Pd + Pm + Pstars + Pbulge + Pdm

    ! Returns extra optional outputs
    if (present(Pd_out)) Pd_out = Pd
    if (present(Pm_out)) Pm_out = Pm
    if (present(Pstars_out)) Pstars_out = Pstars
    if (present(Pbulge_out)) Pbulge_out = Pbulge
    if (present(Pdm_out)) Pdm_out = Pdm

    return
  end function computes_midplane_ISM_pressure_P1

  elemental function computes_midplane_ISM_pressure_from_B_and_rho(B, rho) result(P)
    ! Computes the ISM pressure, knowing the magnetic field and density
    ! Input:  B   -> magnetic field (in microgauss)
    !         rho -> density of diffuse gas (in g/cm^3)
    ! Output: P -> midplane pressure (in erg/cm^3)
    double precision, intent(in) :: B, rho
    double precision :: P

    P = (B*1d-6)**2/4d0/pi + (p_ISM_kappa/3d0*(1.+2.*p_ISM_xi) + 1d0/p_ISM_gamma) &
                             * (p_ISM_sound_speed_km_s*1d5)**2 * rho
    return

  end function computes_midplane_ISM_pressure_from_B_and_rho

  elemental function computes_midplane_ISM_pressure_P2(Sigma_d, &
                                                          Om, S, h_d) result(P)
    double precision, intent(in) :: Sigma_d, Om, S, h_d
    double precision :: P, Sigma_d_SI, Om_SI, S_SI, h_d_SI
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

    P = Om *(Om + S) * (km_SI/kpc_SI)**2
    P = P * Sigma_d_SI * h_d * kpc_SI * convertPressureSItoGaussian

  end function computes_midplane_ISM_pressure_P2

  pure function integrand_stars(x)
    ! Integrand in stars case
    ! Assumes the gas follows a sech^2 vertical profile
    ! Relies on global variable i_hd_to_hstar
    double precision, intent(in) :: x
    double precision :: integrand_stars, sech, y
    y = x*i_hd_to_hstar
    if (.not.p_sech2_profile) then
      integrand_stars = exp(-x)*exp(-y)
    else
      if (abs(y)<500) then
        sech = 1d0/cosh(y)
      else
        sech = 0d0
      endif
      integrand_stars = tanh(x)*sech*sech
    endif
  end function integrand_stars

  pure function integrand_m(x)
    ! Integrand in molecular gas case
    ! Assumes the gas follows a sech^2 vertical profile
    ! Relies on global variable i_hd_to_hm
    double precision, intent(in) :: x
    double precision :: integrand_m, sech, y
    y = x*i_hd_to_hm
    if (.not.p_sech2_profile) then
      integrand_m = exp(-x)*exp(-y)
    else
      if (abs(y)<500) then
        sech = 1d0/cosh(y)
      else
        sech = 0d0
      endif
      integrand_m = tanh(x)*sech*sech
    endif
  end function integrand_m

  pure function integrand_bulge(x)
    ! Integrand in bulge case
    ! Assumes the bulge follows a Hernquist profile
    ! Assumes the diffuse gas follows a sech^2 vertical profile
    ! Relies on global variable i_hd_to_hm
    double precision, intent(in) :: x
    double precision :: y, integrand_bulge
    y = sqrt(i_R_to_Rb**2 + (i_hd_to_Rb*x)**2)
    if (.not.p_sech2_profile) then
      integrand_bulge = exp(-x)/y/(1d0+y)/(1d0+y)/(1d0+y)
    else
      integrand_bulge = tanh(x)/y/(1d0+y)/(1d0+y)/(1d0+y)
    endif
  end function integrand_bulge

  pure function integrand_dm(x)
    ! Integrand in the dark matter halo case
    ! Assumes a Navarro-Frenk-White DM profile
    ! Assumes the diffuse gas follows a sech^2 vertical profile
    ! Relies on global variables i_hd_to_Rsdm and i_R_to_Rsdm
    double precision, intent(in) :: x
    double precision :: y, integrand_dm
    y = sqrt(i_R_to_Rsdm**2 + (i_hd_to_RsDM*x)**2)
    if (.not.p_sech2_profile) then
      integrand_dm = exp(-x)/y/(1d0+y)/(1d0+y)
    else
      integrand_dm = tanh(x)/y/(1d0+y)/(1d0+y)
    endif
  end function integrand_dm

end module pressureEquilibrium
