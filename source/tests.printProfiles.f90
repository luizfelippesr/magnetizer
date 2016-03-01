program testPressure
  use fgsl
  use pressureEquilibrium
!   use outflow
  implicit none

  integer, parameter :: n = 10
  double precision, dimension(n) :: r
  double precision, dimension(n) :: rho_d, h_d, B, Sigma_d, Sigma_star, Rm
  double precision, dimension(n) :: Pgas, Pgrav
  double precision :: rdisk, M_star, M_g, cs, SFR
  double precision, parameter :: G=FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT
  double precision, parameter :: Msun_SI = FGSL_CONST_MKSA_SOLAR_MASS
  double precision, parameter :: kpc_SI = FGSL_CONST_MKSA_PARSEC*1d3
  integer :: i

  ! Andromeda
  M_g = 7.7d8 ! Msun  2010A&A...511A..89C
  M_star = 1.4d11 ! Msun  2010A&A...511A..89C
  cs = 10d0 !km/s
  SFR = 1d0 ! solar mass per year
  ! Half mass radius
  rdisk = 1.711 ! kpc   2003AJ....125..525J
  !  450 arcsec * 0.784 Mpc                http://ned.ipac.caltech.edu/cgi-bin/nDistance?name=MESSIER+031

  ! MW
  M_g = 8d9 ! Msun  ?
  M_star = 2.85d10 ! Msun  2010A&A...511A..89C
  cs = 10d0 !km/s
  SFR = 1.5d0 ! solar mass per year
  rdisk = 5.035 ! kpc
  !  450 arcsec * 0.784 Mpc                http://ned.ipac.caltech.edu/cgi-bin/nDistance?name=MESSIER+031

  do i=1,size(r)
    r(i) = 0.16*i*rdisk
  end do
  B = 1d0

!
!   print *, 'rho (cgs)', pi*G/2d0 *(Msun_SI/kpc_SI/kpc_SI)**2 / 1d4**2 /1.9333333333333 *1e-3
!   print *, 'Pgas (cgs)', pi*G/2d0 *(Msun_SI/kpc_SI/kpc_SI)**2 / 1d4**2  *1e-3 * 1d6**2
!   print *, 'Pgrav (cgs)', pi*G/2d0 *(Msun_SI/kpc_SI/kpc_SI)**2 * 10
!
!   print *, 'h (kpc)', 1d4**2*1.9333333333333/pi/G/(Msun_SI/kpc_SI/kpc_SI) /kpc_SI

  call solves_hytrostatic_equilibrium(rdisk, M_g, M_star, r, B, rho_d, h_d, &
                                      Sigma_star, Sigma_d, Rm)
!   print *, 'rho_comp=', rho_d(3)
!   print *, 'h_d=', h_d(3)

  Pgrav = computes_midplane_ISM_pressure_using_scaleheight(  &
                                          rdisk, Sigma_d, Sigma_star, Rm, h_d)
  Pgas = computes_midplane_ISM_pressure_from_B_and_rho(B, rho_d)
!   print *, 'Pgrav', Pgrav(3)
!   print *, 'Pgas', Pgas(3)
!   print *, Pgrav(3)/Pgas(3)
!   stop
  print *, "#r / kpc     ",  &
           "#rho / g/cm^-3",  &
           "#h/pc      ",  &
           "#1/(1+Rm)  ",  &
           "#P_grav / (erg/cm^3)",  &
           "#P_gas / (erg/cm^3)",  &
           "#Sigma_d / (Msun/kpc^2)",  &
           "#Sigma_s / (Msun/kpc^2)"
  do i=1,size(r)
    print *, r(i), rho_d(i), h_d(i)*1d3, 1d0/(1d0+Rm(i)), Pgrav(i), Pgas(i), Sigma_d(i), Sigma_star(i)
  end do

end program testPressure
