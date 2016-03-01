! Contains the subroutines which compute initial profiles used in the dynamo calculations
module profiles
  use global_input_parameters
  use calc_params
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
  integer :: ialp_k
  double precision, dimension(nx) :: alp_k, alp_k_kms
  double precision, dimension(nx), private :: Om_d, Om_b, Om_h
  double precision, dimension(nx), private :: G_d, G_b, G_h, B2
  double precision, dimension(nx), private :: midplanePressure
  double precision :: rreg
  double precision, parameter :: RREG_TO_RDISK = 0.15
  double precision :: Uphi_halfmass_kms  = -1 ! Negative value when unitialized

contains
  subroutine construct_profiles(B)
    use outflow
    use pressureEquilibrium
    use input_constants
    double precision, dimension(nx), intent(in), optional :: B
    double precision, dimension(nx) :: B_actual
    double precision, dimension(nx) :: B_aux_muG, rho_cgs
    double precision :: rho_ref
    double precision, parameter :: INIT_RHO_TOL = 1d-6
    integer, parameter :: I_REF = 4
    logical :: initial_actual
    integer :: i_halfmass, i

    if (present(B)) then
      B_actual = B
    else
      B_actual = 0.0
    endif

    ! Sets the 'reference radius' to the disk half-mass radius
    r_sol = r_disk/r_max_kpc

    ! ROTATION CURVE
    ! Computes the profile associated with each component
    call disk_rotation_curve(r_kpc, r_disk, v_disk, Om_d, G_d)
    call bulge_rotation_curve(r_kpc, r_bulge, v_bulge, Om_b, G_b)
    call halo_rotation_curve(r_kpc, r_halo, v_halo, nfw_cs1, Om_h, G_h)

    ! Combines the various components
    Om_kmskpc = sqrt( Om_d**2 + Om_b**2 + Om_h**2 )
    G_kmskpc = (Om_d*G_d + Om_b*G_b + Om_h*G_h)/Om_kmskpc

    ! Regularises
    rreg = r_kpc(minloc(abs(r_kpc - RREG_TO_RDISK*r_disk),1))
    call regularize(r_kpc, rreg, Om_kmskpc, G_kmskpc)

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

    ! TURBULENT SCALE PROFILE
    if (.not.Var_l) then
      l = l_sol
      dldr = 0.d0
    else
      l = l_sol * exp((r-r_sol)/r_l)
      dldr = l / r_l
    endif
    l_kpc = l * h0_kpc / h0

    ! TURBULENT VELOCITY PROFILE
    v_kms = p_ISM_sound_speed_km_s * p_ISM_kappa
    v = v_kms / h0_km * h0 * t0_s / t0

    ! Solves for density and height
    call solves_hytrostatic_equilibrium(r_disk, Mgas_disk, Mstars_disk, r, &
                                        B_actual, rho_cgs, h_kpc)
    ! NUMBER DENSITY PROFILE
    n_cm3 = rho_cgs/Hmass
    n = n_cm3 / n0_cm3 * n0

    ! EQUIPARTITION MAGNETIC FIELD STRENGTH PROFILE
    Beq = sqrt(4d0*pi*n) * v  !Formula for equiparition field strength
    Beq_mkG = Beq * B0_mkG / B0

    ! SCALE HEIGHT PROFILE
    h = h_kpc*h0/h0_kpc

    ! VERTICAL VELOCITY PROFILE
    Uz_kms = outflow_speed(r_kpc, rho_cgs, h_kpc, v_kms, r_disk, SFR, Mgas_disk, Mstars_disk)
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
      alp_k = C_alp*l**2/h*Om  !Decreasing with radius
    endif
    if (Alp_ceiling) then
      do ialp_k=1,nx
        if (alp_k(ialp_k) > alpceil*v(ialp_k)) then
          alp_k(ialp_k) = alpceil*v(ialp_k)
        endif
      enddo
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

  end subroutine construct_profiles

  subroutine disk_rotation_curve(rx, r_disk, v_disk, Omega, Shear)
    ! Computes the rotation curve associated with an exponential disk
    ! Input:  rx -> the radii where the rotation curve will be computed
    !         r_disk ->  the half mass radius of the disk
    !         v_disk -> the circular velocity at r_disk
    ! Output: Omega -> angular velocity profile
    !         Shear  -> shear profile
    ! Info:   V^2 \propto y^2 [I0(y)K0(y)-I1(y)K1(y)] where y=rx/r_s
    !         r_s --> the scale radius
    !         For the shear, see  http://is.gd/MnDKyS
    ! Ref:    Binney & Tremaine or  Mo, Bosch & White

    use Bessel_Functions
    use input_constants
    implicit none
    double precision, intent(in) :: r_disk, v_disk
    double precision, dimension(:), intent(in)  :: rx
    double precision, dimension(size(rx)) :: A
    double precision, dimension(size(rx)),intent(out) :: Omega, Shear
    double precision, parameter :: rmin_over_rmax=0.1
    double precision, parameter :: TOO_SMALL=2e-7 ! chosen empirically
    double precision, parameter :: rs_to_r50=constDiskScaleToHalfMassRatio
    double precision, dimension(size(rx)) :: y
    double precision :: rs, y50
    integer :: i

    ! Scale radius
    rs = r_disk*constDiskScaleToHalfMassRatio
    ! (Half-mass radius) / (scale radius)
    y50 = 1d0/constDiskScaleToHalfMassRatio

    ! Traps disks of negligible size
    if (r_disk < r_max_kpc*rmin_over_rmax) then
      Shear = 0
      Omega = 0
      return
    end if

    ! Sets y=r/rs, trapping r~0 (which breaks the evaluation
    ! of the Bessel functions)
    where (rx > TOO_SMALL)
      ! The absolute value is taken to deal with the ghost zone
      y = abs(rx) / rs
    elsewhere
      y = TOO_SMALL / rs
    endwhere

    do i=1,size(rx)
      A(i) = (  I0(y(i)) * K0(y(i))            &
              - I1(y(i)) * K1(y(i)) )          &
              ! Check this!
            / ( I0(y50) * K0(y50)  &
              - I1(y50) * K1(y50))
      A(i) = sqrt(A(i))
    end do
    Omega = A*v_disk/r_disk

    do i=1,size(rx)
      Shear(i) =    I1(y(i)) * K0(y(i))                       &
                  - I0(y(i)) * K1(y(i))                       &
                  -0.5d0 * K1(y(i)) *( I0(y(i)) + I2(y(i)) )  &
                  +0.5d0 * I1(y(i)) *( K0(y(i)) + K2(y(i)) )
    end do

    Shear = Shear/A/2d0*v_disk/r_disk
    return
  end subroutine disk_rotation_curve

  subroutine bulge_rotation_curve(rx, r_bulge, v_bulge, Omega, Shear)
    ! Computes the rotation curve associated with an Hernquist profile
    ! Input:  rx -> the radii where the rotation curve will be computed
    !         r_bulge ->  the half mass radius of the spheroid
    !         v_bulge -> the circular velocity at r_bulge
    ! Output: Omega -> angular velocity profile
    !         Shear  -> shear profile
    ! Info:   V^2 \propto y/(y+1)^2  where y=r/a
    ! Ref:    http://adsabs.harvard.edu/abs/1990ApJ...356..359H

    implicit none
    double precision, intent(in) :: r_bulge, v_bulge
    double precision, dimension(:), intent(in)  :: rx
    double precision, dimension(size(rx)) :: v, dvdr
    double precision, dimension(size(rx)),intent(out) :: Omega, Shear
    double precision, parameter :: a_to_r50 = 1.0d0/(sqrt(2.0d0)-1.0d0)
    double precision, dimension(size(rx)) :: y
    double precision :: A
    integer :: i

    y = abs(rx)/ (r_bulge/a_to_r50)

    A = v_bulge * (1d0 + a_to_r50) / sqrt(a_to_r50)

    ! v =  (y/a_to_r50) * (a_to_r50 +1d0)**2 * (y +1d0)**(-2)
    ! v = v_bulge * sqrt(v)
    v = A * sqrt(y)/(y+1d0)

    Omega = v/abs(rx)

    dvdr = A/2d0/((y+1d0)*sqrt(y)) - A*sqrt(y)/(1d0+y)**2


    Shear = dvdr - Omega

  end subroutine bulge_rotation_curve

  subroutine halo_rotation_curve(rx, r_halo, v_halo, nfw_cs1, Omega, Shear)
    ! Computes the rotation curve associated with an NFW halo
    ! Warning: This ignores effects of adiabatic contraction!
    !          It need to be accounted for later
    ! Input: rx -> the radii where the rotation curve will be computed
    !        r_halo -> the virial radius of the halo
    !        v_halo -> the circular velocity at r_halo
    !        nfw_cs1 -> 1/c_s -> inverse of the NFW concentration parameter
    !                        i.e. (NFW scale radius)/(virial radius)
    ! Output: Omega -> angular velocity profile
    !         Shear  -> shear profile
    ! Info:  V^2 \propto {ln[(cs1+y)/cs1] - y/(cs1+y)} /
    !                    {ln[(cs1+1)/cs1] - 1/(cs1+1)} / y
    ! Ref: NFW profile

    implicit none
    double precision, intent(in) :: r_halo, v_halo, nfw_cs1
    double precision, dimension(:), intent(in)  :: rx
    double precision, dimension(size(rx)), intent(out) :: Omega, Shear
    double precision, dimension(size(rx)) :: v, dv
    double precision, dimension(size(rx)) :: y
    double precision, dimension(size(rx)) :: B
    double precision :: A
    integer :: i

    y = abs(rx) / r_halo

    A = 1.0 / ( log((nfw_cs1+1d0)/nfw_cs1) - 1d0/(nfw_cs1+1d0) )
    A = sqrt(A)

    B = (log((nfw_cs1+y)/nfw_cs1) - y/(nfw_cs1+y)) / y
    v = A * v_halo * sqrt(B)

    dv = 1.0/(nfw_cs1+y)**2-(log((nfw_cs1+y)/nfw_cs1)-y/(nfw_cs1+y))/y**2/2.0
    dv = dv * A * v_halo

    Omega = v/rx
    Shear = dv - Omega

  end subroutine halo_rotation_curve

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
    double precision, dimension(size(rx)) :: rxi_over_r_2
    double precision :: Omega_xi
    double precision, parameter :: small_factor=1e-10
    integer i
    ! Finds the index of r=r_xi and sets Omega_xi
    Omega_xi = Omega( minloc(abs(rx-r_xi),1) )

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

    ! Regularises Shear
    Shear = exp_minus_rxi_over_r*(2.0*rxi_over_r_2*(Omega-Omega_xi)+Shear)
    ! Regularises Omega
    Omega = exp_minus_rxi_over_r*(Omega-Omega_xi)+Omega_xi

  end subroutine regularize
end module profiles
