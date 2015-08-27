! Contains the subroutines which compute initial profiles used in the dynamo calculations
module profiles
  use global_input_parameters
  use calc_params
  use grid
!
  implicit none
!
  double precision, dimension(nx) :: h, h_kpc
  double precision, dimension(nx) :: om, G, om_kmskpc, G_kmskpc
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
!
  contains
    subroutine construct_profiles
!     SCALE HEIGHT PROFILE
      if (Flaring) then
        h= h_sol*dexp((r-r_sol)/r_h)
      else
        h= h_sol
      endif
      h_kpc=h*h0_kpc/h0

      ! ROTATION CURVE
      if (Om_Brandt) then
        om=om0/((1.d0+(r/r_om)**2))**(1.d0/2)  !Brandt angular velocity profile
        if (Shear) then
          G=-om0*(r/r_om)**2/(1.d0+(r/r_om)**2)**(3.d0/2)
        else
          G=0.
        endif
!       else if (three_components_rotation_curve)
!         call compute_rotation_curve(M_bulge,r_bulge,M_disk,r_disk,Mhalo,cs, om, G)
      endif
      om_kmskpc=om*h0_km/h0_kpc/t0_s*t0
      G_kmskpc=  G*h0_km/h0_kpc/t0_s*t0
!
!     VERTICAL VELOCITY PROFILE
      if (.not.Var_Uz) then
        Uz= Uz_sol  !No variation of Uz
      else
        Uz= Uz_sol*dexp(-(r-r_sol)/r_Uz)  !Decreasing with radius according to exponential
      endif
      Uz_kms=Uz*h0_km/h0/t0_s*t0
!
!     RADIAL VELOCITY PROFILE
      Ur=Ur_sol
      dUrdr=0.d0
      d2Urdr2=0.d0
      Ur_kms= Ur*h0_km/h0/t0_s*t0
!
!     NUMBER DENSITY PROFILE
      if (.not.Var_n) then
        n=n_sol
      else
        n=n_sol*exp(-(r-r_sol)/r_n)
      endif
      n_cm3=n*n0_cm3/n0
!
!     TURBULENT SCALE PROFILE
      if (.not.Var_l) then
        l=l_sol
        dldr=0.d0
      else
        l=l_sol*exp((r-r_sol)/r_l)
        dldr=l/r_l
      endif
      l_kpc=l*h0_kpc/h0
!
!     RMS TURBULENT VELOCITY PROFILE
      if (.not.Var_v) then
        v=v_sol
        dvdr=0.d0
      else
        v=v_sol*exp(-(r-r_sol)/r_v)
        dvdr=v/r_v
      endif
      v_kms=v*h0_km/h0/t0_s*t0
!
!     TURBULENT DIFFUSIVITY PROFILE
      etat=1.d0/3*l*v  !Formula for etat from mixing length theory
      etat_cm2s=etat*h0_cm**2/h0**2/t0_s*t0
      etat_kmskpc=etat*h0_km*h0_kpc/h0**2/t0_s*t0
!
!     TURBULENT CORRELATION TIME PROFILE
      tau=        ctau*l/v  !Formula for tau from mixing length theory
      tau_Gyr=    ctau*tau*t0_Gyr/t0
      tau_s=      ctau*tau*t0_s/t0
      tau_sol=    ctau*l_sol/v_sol
      tau_sol_Gyr=ctau*tau_sol*t0_Gyr/t0
      tau_sol_s=  ctau*tau_sol*t0_s/t0
!
!     EQUIPARTITION MAGNETIC FIELD STRENGTH PROFILE
      Beq=dsqrt(4*pi*n)*v  !Formula for equiparition field strength
      Beq_mkG=Beq*B0_mkG/B0
!
!     KINETIC ALPHA PROFILE
      if (.not.Krause) then
        alp_k= C_alp  !No variation of alpha
      else
        alp_k= C_alp*l**2/h*om  !Decreasing with radius
      endif
      if (Alp_ceiling) then
        do ialp_k=1,nx
          if (alp_k(ialp_k)>alpceil*v(ialp_k)) then
            alp_k(ialp_k)=alpceil*v(ialp_k)
          endif
        enddo
      endif
      alp_k_kms=alp_k*h0_km/h0/t0_s*t0
    end subroutine construct_profiles

!     subroutine compute_rotation_curve(M_bulge,r_bulge,M_disk,r_disk,Mhalo,cs,om,G)
!       implicit none
!       double precision, intent(in) :: M_bulge,r_bulge,M_disk,r_disk,Mhalo,cs
!       double precision, dimension(nx), intent(out) :: om, G
!
!     end subroutine compute_rotation_curve

end module profiles
