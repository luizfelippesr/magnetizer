! Contains the subroutines which compute initial profiles used in the dynamo calculations
module profiles
  use grid
  implicit none
  double precision, dimension(:), allocatable :: h, n
  double precision, dimension(:), allocatable :: l, dldr, l_kpc
  double precision, dimension(:), allocatable :: v, dvdr
  double precision, dimension(:), allocatable :: etat, tau, Beq, alp_k
  double precision, dimension(:), allocatable :: Uz, Ur, dUrdr
  double precision, dimension(:), allocatable :: Om, G
  double precision :: delta_r
  private :: prepare_profiles_module_public_variables
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
    use messages, only: error_message
    use deriv, only: xder
    use grid, only: x
    double precision :: Uphi_halfmass_kms  = -1 ! Negative value when unitialized
    double precision, dimension(nx), intent(in), optional :: B
    double precision, dimension(nx) :: d2Urdr2, Ur_kms, Uz_kms
    double precision, dimension(nx) :: Om_kmskpc, G_kmskpc
    double precision, dimension(nx) :: etat_cm2s, etat_kmskpc
    double precision, dimension(nx) :: h_kpc, n_cm3, v_kms, alp_k_kms
    double precision, dimension(nx) :: Om_d, Om_b, Om_h, Beq_mkG
    double precision, dimension(nx) :: G_d, G_b, G_h
    double precision, dimension(nx) :: B_actual, tau_Gyr, tau_s
    double precision, dimension(nx) :: rho_cgs
    double precision, dimension(nx) :: Sigma_d, Sigma_star, Pgas, Rm, Pgrav
    double precision, dimension(nx,3) :: all_roots
    double precision, parameter :: P_TOL=1e-10
    double precision :: rreg
    double precision :: r_disk_min
    integer :: i_halfmass

    call prepare_profiles_module_public_variables()

    if (present(B)) then
      B_actual = B
    else
      B_actual = 0.0
    endif

    construct_profiles = .true.

    ! Sets the minimum radius to be followed (for the disk rotation curve)
    r_disk_min = r_max_kpc*rmin_over_rmax
    ! Adjust units of delta_r (which is used to get the seed field)
    delta_r = p_delta_r_kpc/r_max_kpc

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
    G_kmskpc = xder(Om_kmskpc)*x

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
    Ur = 0d0
    dUrdr = 0.d0
    d2Urdr2 = 0.d0

    Ur_kms = Ur*h0_km/h0/t0_s*t0

    ! TURBULENT VELOCITY PROFILE
    v_kms = p_ISM_sound_speed_km_s * p_ISM_kappa
    v = v_kms / h0_km * h0 * t0_s / t0
    dvdr = 0.0

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
    if (.not.p_use_legacy_cubic_solver) then
      ! Computes h and rho solving for hydrostatic equilibrium
      call solve_hydrostatic_equilibrium_numerical(abs(r_kpc), B_actual, &
                                                   Om_kmskpc, G_kmskpc, &
                                                   rho_cgs, h_kpc, Rm, &
                                                   Sigma_d_out=Sigma_d)
    else
      if (.not.p_check_hydro_solution) then
        call solve_hydrostatic_equilibrium_cubic(r_disk, Mgas_disk, Mstars_disk, &
                            abs(r_kpc), B_actual, rho_cgs, h_kpc, Rm_out=Rm)
      else
        ! If required, checks the solution.
        ! (This is a slow step. Should be used for debugging only.)
        call solve_hydrostatic_equilibrium_cubic(r_disk, Mgas_disk, Mstars_disk, &
                                        abs(r_kpc), B_actual, rho_cgs, h_kpc, &
                                            Sigma_star, Sigma_d, Rm, all_roots)
        ! Computes the midplane pressure, from gravity
        Pgrav = computes_midplane_ISM_pressure_using_scaleheight(  &
                                          r_disk, Sigma_d, Sigma_star, Rm, h_kpc)
        ! Computes the midplane pressure, from the density
        Pgas = computes_midplane_ISM_pressure_from_B_and_rho(B_actual, rho_cgs)

        if (any(abs(Pgrav-Pgas)/Pgas>P_TOL)) then
          call error_message('construct_profiles', &
                            'Invalid solution for the hydrostatic equilibrium.',&
                            code='P')
          construct_profiles = .false.
        end if
      end if
    endif

    ! Traps possible aberrant values for the scale height
    if (any(h_kpc<0)) then
      call error_message('construct_profiles', 'Negative scaleheight detected.',&
                         code='H')
      construct_profiles = .false.
    endif

    ! Stricter(ish) trap: h/r<1.5 at a half mass radius
    if (h_kpc(i_halfmass)>r_disk) then
      call error_message('construct_profiles','Huge scaleheight detected.', &
                         code='h')
!       construct_profiles = .false.
    endif

    ! NUMBER DENSITY PROFILE
    ! At some point, we may want to improve this accounting for the metallicity of
    ! of the galaxy (e.g. using mean molecular weight of the diffuse phase)
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
    double precision, parameter :: small_factor=1e-15
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
    where (rxi_over_r_2 < 1e12)
      exp_minus_rxi_over_r = exp( -rxi_over_r_2 )
    elsewhere
      exp_minus_rxi_over_r =0
    endwhere

    ! Regularises Shear
    Shear = exp_minus_rxi_over_r*(2.0*rxi_over_r_2*(Omega-Omega_xi)+Shear)
    ! Regularises Omega
    Omega = exp_minus_rxi_over_r*(Omega-Omega_xi)+Omega_xi

  end subroutine regularize

  subroutine prepare_profiles_module_public_variables()
    ! Helper to check whether the module variables have the correct shape
    call check_allocate(h)
    call check_allocate(n)
    call check_allocate(l)
    call check_allocate(dldr)
    call check_allocate(l_kpc)
    call check_allocate(v)
    call check_allocate(dvdr)
    call check_allocate(etat)
    call check_allocate(tau)
    call check_allocate(Beq)
    call check_allocate(alp_k)
    call check_allocate(Uz)
    call check_allocate(Ur)
    call check_allocate(dUrdr)
    call check_allocate(Om)
    call check_allocate(G)
  end subroutine prepare_profiles_module_public_variables
end module profiles
