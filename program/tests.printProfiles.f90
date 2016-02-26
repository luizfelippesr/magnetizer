program testPressure
  use pressureEquilibrium
!   use outflow
  implicit none

  integer, parameter :: n = 25
  double precision, dimension(n) :: r
  double precision, dimension(n) :: rho_d, h_d, B, Sigma_d, Sigma_star, Rm
  double precision, dimension(n) :: Pgas, Pgrav
  double precision :: rdisk, M_star, M_g, cs, SFR
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
  B = 0.0

  call solves_hytrostatic_equilibrium(rdisk, M_g, M_star, r, B, rho_d, h_d, &
                                      Sigma_star, Sigma_d, Rm)

  Pgrav = computes_midplane_ISM_pressure_using_scaleheight(   &
                                          rdisk, Sigma_d, Sigma_star, Rm, h_d)
  Pgas = computes_midplane_ISM_pressure_from_B_and_rho(B, rho_d)

  print *, "#r/kpc     ",  &
           "#(rho/m_H)/cm^-3 ",  &
           "#h/pc      ",  &
           "#1/(1+Rm)  ",  &
           "#P_grav/(erg/cm^2)",  &
           "#P_gas/(erg/cm^2)",  &
           "#Sigma_d/(Msun/pc^2)",  &
           "#Sigma_s/(Msun/pc^2)"
  do i=1,size(r)
    print *, r(i), rho_d(i)/1.67372d-24, h_d(i)*1d3, 1d0/(1d0+Rm(i)), Pgrav(i), Pgas(i), Sigma_d(i)/1d6, Sigma_star(i)/1d6
  end do

end program testPressure
