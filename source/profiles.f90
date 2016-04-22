! Contains the subroutines which compute initial profiles used in the dynamo calculations
module profiles
  use grid
  implicit none

  double precision, dimension(nx) :: h, h_kpc
  double precision, dimension(nx) :: Om, G, Om_kmskpc, G_kmskpc
  double precision, dimension(nx) :: Uz, Uz_kms
  double precision, dimension(nx) :: Ur, dUrdr, d2Urdr2, Ur_kms
  double precision, dimension(nx) :: n, n_cm3
  double precision, dimension(nx) :: l, l_kpc, dldr
  double precision, dimension(nx) :: v, v_kms, dvdr
  double precision, dimension(nx) :: etat, etat_cm2s, etat_kmskpc
  double precision, dimension(nx) :: tau, tau_Gyr, tau_s
  double precision :: tau_sol, tau_sol_Gyr, tau_sol_s
  double precision, dimension(nx) :: Beq, Beq_mkG
  double precision, dimension(nx) :: alp_k, alp_k_kms
  double precision, dimension(nx), private :: Om_d, Om_b, Om_h
  double precision, dimension(nx), private :: G_d, G_b, G_h
  double precision, private :: rreg
  double precision :: Uphi_halfmass_kms  = -1 ! Negative value when unitialized

contains
  logical function construct_profiles(B)
    ! Computes the radial of variation of all relevant physical constants
    ! Input: B -> magnetic field profile (if abscent, zero is used)
    ! Other input variables are obtained from the public global variables
    ! defined in the input_parameters module.
    ! The output is written in public global variables
    ! Returns True if succesful. False otherwise.
    use global_input_parameters
    use calc_params
    use rotationCurves
    use outflow
    use pressureEquilibrium
    use input_constants
    use messages
    double precision, dimension(nx), intent(in), optional :: B
    double precision, dimension(nx) :: B_actual
    double precision, dimension(nx) :: rho_cgs
    double precision, dimension(nx) :: Sigma_d, Sigma_star, Pgrav, Pgas, Rm
    double precision, dimension(nx,3) :: all_roots
    double precision, parameter :: P_TOL=1e-10
    double precision :: r_disk_min
    integer :: i_halfmass

    if (present(B)) then
      B_actual = B
    else
      B_actual = 0.0
    endif

    construct_profiles = .true.

    ! Sets the 'reference radius' to the disk half-mass radius
    ! (actually, the maximum half mass radius over history)
    r_sol = r_disk/r_max_kpc
    ! Sets the minimum radius to be followed (for the disk rotation curve)
    r_disk_min = r_max_kpc*rmin_over_rmax

    ! ROTATION CURVE
    ! Computes the profile associated with each component
    call disk_rotation_curve(r_kpc, r_disk, r_disk_min, v_disk, Om_d, G_d)
    call bulge_rotation_curve(r_kpc, r_bulge, v_bulge, Om_b, G_b)
    call halo_rotation_curve(r_kpc, r_halo, v_halo, nfw_cs1, Om_h, G_h)

    ! Combines the various components
    Om_kmskpc = sqrt( Om_d**2 + Om_b**2 + Om_h**2 )
    G_kmskpc = (Om_d*G_d + Om_b*G_b + Om_h*G_h)/Om_kmskpc

    ! Regularises
    rreg = r_kpc(minloc(abs(r_kpc - p_rreg_to_rdisk*r_disk),1))
    call regularize(abs(r_kpc), rreg, Om_kmskpc, G_kmskpc)

    ! If required, avoid bumps in the shear
    if (.not.p_allow_positive_shears) then
      where (G_kmskpc>0d0)
        G_kmskpc = 0d0
      end where
    endif

    ! Adjusts units to code units (set in units module)
    Om = Om_kmskpc/h0_km*h0_kpc*t0_s/t0
    G  = G_kmskpc /h0_km*h0_kpc*t0_s/t0

    ! EXTRA (for debugging/diagnostic): computes quantities at r_disk
    ! Finds the position in the grid closest to r_disk (disk half-mass radius)
    i_halfmass = minloc(abs(r_kpc - r_disk), 1)
    ! Computes the rotation velocity at the disk half-mass radius
    Uphi_halfmass_kms = Om_kmskpc(i_halfmass)*r_kpc(i_halfmass)

    ! RADIAL VELOCITY PROFILE (not currently being used?)
    Ur = Ur_sol
    dUrdr = 0.d0
    d2Urdr2 = 0.d0

    Ur_kms = Ur*h0_km/h0/t0_s*t0

    ! TURBULENT VELOCITY PROFILE
    v_kms = p_ISM_sound_speed_km_s * p_ISM_kappa
    v = v_kms / h0_km * h0 * t0_s / t0

    ! Traps case of negligible disks (produced by a recent major merger)
    if (r_disk < r_disk_min .or. Mgas_disk< Mgas_disk_min) then
      ! Sets the profiles to zero and returns
      n = 0.0
      l = 0.0
      h = 0.0
      Uz = 0.0
      etat = 0.0
      alp_k = 0.0
      return
    endif


    ! Solves for density and height
    if (.not.p_check_hydro_solution) then
      call solves_hytrostatic_equilibrium(r_disk, Mgas_disk, Mstars_disk, &
                          abs(r_kpc), B_actual, rho_cgs, h_kpc, Rm_out=Rm)
      if (any(h_kpc<0)) then
        call message('construct_profiles: Error. Negative scaleheight detected.')
        construct_profiles = .false.
      endif
      if (any(h_kpc>1d3)) then
        call message('construct_profiles: Error. Huge scaleheight detected.')
        construct_profiles = .false.
      endif
    else
      ! If required, checks the solution.
      ! (This is a slow step. Should be used for debugging only.)
      call solves_hytrostatic_equilibrium(r_disk, Mgas_disk, Mstars_disk, &
                                      abs(r_kpc), B_actual, rho_cgs, h_kpc, &
                                           Sigma_star, Sigma_d, Rm, all_roots)
      ! Computes the midplane pressure, from gravity
      Pgrav = computes_midplane_ISM_pressure_using_scaleheight(  &
                                        r_disk, Sigma_d, Sigma_star, Rm, h_kpc)
      ! Computes the midplane pressure, from the density
      Pgas = computes_midplane_ISM_pressure_from_B_and_rho(B_actual, rho_cgs)

      if (any(abs(Pgrav-Pgas)/Pgas>P_TOL)) then
        call message('construct_profiles: Warning. Invalid solution for the hydrostatic equilibrium.')
        construct_profiles = .false.
      end if
    end if

    ! NUMBER DENSITY PROFILE
    ! At some point, we have to improve this accounting for the metallicity of
    ! of the galaxy (i.e. multiplying by the mean molecular weight)
    n_cm3 = rho_cgs/Hmass
    n = n_cm3 / n0_cm3 * n0

    ! TURBULENT SCALE PROFILE
    if (p_use_fixed_turbulent_to_scaleheight_ratio) then
      l_kpc = h_kpc * p_turbulent_to_scaleheight_ratio

      l = l_kpc*h0/h0_kpc
      dldr = 0.d0 !TODO make this derivative work!?
                ! It is not clear that this would be a good idea...
    else
      l_kpc = p_ISM_turbulent_length
      if (.not.p_limit_turbulent_scale) then
        l = l_kpc*h0/h0_kpc
        dldr = 0.d0
      else
        ! Enforces l<=h
        where (l_kpc>h_kpc)
          l_kpc = h_kpc
        endwhere
        l = l_kpc*h0/h0_kpc
        dldr = 0.d0 !TODO make this derivative work!?
                ! It is not clear that this would be a good idea...
      endif
    endif

    ! EQUIPARTITION MAGNETIC FIELD STRENGTH PROFILE
    Beq = sqrt(4d0*pi*n) * v  !Formula for equiparition field strength
    Beq_mkG = Beq * B0_mkG / B0

    ! SCALE HEIGHT PROFILE
    h = h_kpc*h0/h0_kpc

    ! VERTICAL VELOCITY PROFILE
    Uz_kms = outflow_speed(r_kpc, rho_cgs, h_kpc, v_kms, &
                                                    r_disk, v_disk, SFR, Rm)
    Uz = Uz_kms/h0_km*h0*t0_s/t0

    ! TURBULENT DIFFUSIVITY PROFILE
    etat = 1.d0/3*l*v  !Formula for etat from mixing length theory
    etat_cm2s = etat*h0_cm**2/h0**2/t0_s*t0
    etat_kmskpc = etat*h0_km*h0_kpc/h0**2/t0_s*t0

    ! TURBULENT CORRELATION TIME PROFILE
    tau=          ctau*l/v  !Formula for tau from mixing length theory
    tau_Gyr=      ctau*tau*t0_Gyr/t0
    tau_s=        ctau*tau*t0_s/t0
    tau_sol=      ctau*l_sol/v_sol
    tau_sol_Gyr = ctau*tau_sol*t0_Gyr/t0
    tau_sol_s=    ctau*tau_sol*t0_s/t0

    ! KINETIC ALPHA PROFILE
    if (.not.Krause) then
      alp_k = C_alp  !No variation of alpha
    else
      where (h>1d-70) ! Just to avoid numerical problems
        alp_k = C_alp*l**2/h*Om
      elsewhere
        alp_k = alpceil*v ! Same ceiling as below
      endwhere
    endif
    if (Alp_ceiling) then
      where (alp_k > alpceil*v)
        alp_k = alpceil*v
      endwhere
    endif
    alp_k_kms = alp_k*h0_km/h0/t0_s*t0

    ! Overrides previous definitions if these options were selected in the
    ! global parameters file
    if (.not.Turb_dif) then
        etat=0.d0
      endif

    if (.not.Advect) then
      om=0.d0
    endif

  end function construct_profiles


  subroutine regularize(rx, r_xi, Omega, Shear)
    ! Regularizes the angular velocity and shear, making the centre of the
    ! galaxy behave as a rotating solid body (with constant angular velocity
    ! Omega = Omega(r_xi)).
    !
    ! Input: rx -> the radii where the rotation curve will be computed
    ! Inout: Omega -> angular velocity profile
    !        Shear  -> shear profile
    ! Info:  \Omega(r) = exp(-r_\xi/r) [\tilde\Omega(r) -\Omega(r_\xi)]
    !                    + \Omega(r_\xi)
    !
    double precision, dimension(:), intent(in)  :: rx
    double precision, intent(in)  :: r_xi
    double precision, dimension(size(rx)), intent(inout) :: Omega, Shear
    double precision, dimension(size(rx)) :: exp_minus_rxi_over_r
!     double precision, dimension(size(rx)) :: exp_things
    double precision, dimension(size(rx)) :: rxi_over_r_2
    double precision :: Omega_xi
    double precision, parameter :: small_factor=1e-10
    double precision, parameter :: s=1d0
    ! Finds the index of r=r_xi and sets Omega_xi
    Omega_xi = 1.5d0*Omega( minloc(abs(rx-r_xi),1) )

    ! Be cautious about very small radii
    ! (to avoid NaNs still keeping Omega -> Omega_xi for r->0)
    where (rx > small_factor*r_xi)
      rxi_over_r_2 = (r_xi/rx)**2
    elsewhere
      exp_minus_rxi_over_r = 0.0
      Omega = 0.0
      Shear = 0.0
      rxi_over_r_2 = 0.0
    endwhere
    ! Explicitly avoids underflows too
    where (rxi_over_r_2 < 1e8)
      exp_minus_rxi_over_r = exp( -rxi_over_r_2 )
    elsewhere
      exp_minus_rxi_over_r =0
    endwhere

!     exp_things = exp(s*(rx-r_xi)**2)
!     exp_things = exp(-r_xi/2d0/rx)
!     Omega = (Omega*exp_things+Omega_xi)/(1d0+exp_things)

    ! Regularises Shear
    Shear = exp_minus_rxi_over_r*(2.0*rxi_over_r_2*(Omega-Omega_xi)+Shear)
    ! Regularises Omega
    Omega = exp_minus_rxi_over_r*(Omega-Omega_xi)+Omega_xi

  end subroutine regularize
end module profiles
