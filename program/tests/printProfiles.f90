program testPressure
  use pressureEquilibrium
  use outflow
  implicit none

  integer, parameter :: n = 25
  double precision, dimension(n) :: r
  double precision, dimension(n) :: P, rho, h, v
  double precision :: rdisk, Mstars, Mgas, cs, SFR
  integer :: i

  ! Andromeda
  ! Half mass radius
  rdisk = 1.711 ! kpc   2003AJ....125..525J
  !  450 arcsec * 0.784 Mpc                http://ned.ipac.caltech.edu/cgi-bin/nDistance?name=MESSIER+031
  Mgas = 7.7d8 ! Msun  2010A&A...511A..89C
  Mstars = 1.4d11 ! Msun  2010A&A...511A..89C
  cs = 10d0 !km/s
  SFR = 1d0 ! solar mass per year

  do i=1,size(r)
    r(i) = 0.2*i*rdisk
  end do

  call set_density_procedures('simple','simple')

  P = midplane_pressure(r, rdisk, Mgas, Mstars, sigma_scaling=.true.)
  rho = midplane_density(r, P, 0d0*r, cs, 5d0/3d0, 1d0)
  h = scaleheight(r, rdisk, Mgas, rho)
  v = outflow_speed(r, rho, h, cs*r/r, rdisk, SFR, Mgas, Mstars, 'superbubble_simple')

  print *, '#r(i)  ', '  P(i)', "  n'(i)", '    n','    h','    v'
  do i=1,size(r)
    print *, r(i), P(i), P(i)/(cs*1d5)**2/1.67372d-24, rho(i)/1.67372d-24, h(i)*1d3, v(i)
  end do

end program testPressure
