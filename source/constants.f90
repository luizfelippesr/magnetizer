! Contains math constants and unit conversions

module math_constants

  implicit none

  double precision, parameter :: pi=     3.14156295358d0
  double precision, parameter :: mp_g=   1.67262158d-24  !proton mass in grams
  double precision, parameter :: cm_kpc= 3.08567758d21  !number of cm in 1 kpc
  double precision, parameter :: km_kpc= 3.08567758d16  !number of km in 1 kpc
  double precision, parameter :: cm_km=  1.d5  !number of cm in 1 km
  double precision, parameter :: s_Gyr=  1.d9*365.25d0*24*3600  !number of s in 1 Gyr
  double precision, parameter :: mkG_G=  1.d6  !number of mkG in 1 G

end module math_constants


module units  !Specifies code units
  use math_constants

  implicit none

  ! DIMENSIONLESS UNITS
  double precision, parameter :: etat0 = 1.0d0
  double precision, parameter :: h0 = 1.0d0
  double precision, parameter :: B0 = 1.0d0
  double precision, parameter :: t0 = h0**2/etat0
  double precision, parameter :: n0 = B0**2/h0**2*t0**2

  ! DIMENSIONAL UNITS
  ! Unit of etat, corresp to a typical value, in km*kpc/s
  double precision, parameter :: etat0_kmskpc = 1.d0/3*0.1d0*10.d0
  ! Unit of etat, corresp to a typical value, in cm^2/s
  double precision, parameter :: etat0_cm2s = etat0_kmskpc*cm_km*cm_kpc
  ! Unit of length, corresp to a typical half-disk thickness, in kpc
  double precision, parameter :: h0_kpc = 1d0
  ! Unit of length, corresp to a typical half-disk thickness, in km
  double precision, parameter :: h0_km=h0_kpc*km_kpc
  ! Unit of length, corresp to a typical half-disk thickness, in cm
  double precision, parameter :: h0_cm=h0_km*cm_km
  ! Unit of magnetic field, in mkG
  double precision, parameter :: B0_mkG = 1.0d0
  ! Unit of time, corresp to a typical vert diff time, in s*kpc/km
  double precision, parameter :: t0_kpcskm = h0_kpc**2/etat0_kmskpc
  ! Unit of time in Gyr
  double precision, parameter :: t0_Gyr = t0_kpcskm/s_Gyr*km_kpc
  ! Unit of time in s
  double precision, parameter :: t0_s = t0_kpcskm*km_kpc
  ! Unit of number density, in cm^-3
  double precision, parameter :: n0_cm3 = (B0_mkG/mkG_G)**2/h0_cm**2*t0_s**2/mp_g

end module units


module input_constants  ! Contains physical, model-dependent, constants

  use math_constants
  implicit none

  ! TURBULENCE
  !Ratio of tau (eddy turnover time) to correlation time of the turbulence
  double precision, parameter :: ctau = 1.d0
  ! ADVECTION NO-z
  !No-z correction term for terms involving Uz*Br, Uz*Bp
  double precision, parameter :: C_U = 0.25d0
  !No-z correction term for terms involving Uz*alp_m
  double precision, parameter :: C_a = 1.d0!0.3d0
  ! DIFFUSIVE MAGNETIC HELICITY FLUX
  !No-z correction term for terms involving diffusive flux
  double precision, parameter :: C_d = -pi**2
  ! RESISTIVITY
  !1.e-5 !Inverse magnetic Reynolds number
  double precision, parameter :: Rm_inv = 0.d0
   !Relevant only if module Alp_ceiling=1; Limits maximum value of alpha to ceiling*v_kms
  double precision, parameter :: alpceil = 1d0 !0.5d0
  ! RANDOM MAGNETIC FIELD
  !Strength of the rms random magnetic field brms in units of Beq
  double precision, parameter :: fmag = 0.5d0
  ! Exponential disk properties
  ! Ratio between the scale radius and the half mass radius
  double precision, parameter :: constDiskScaleToHalfMassRatio = 1d0/1.678346990d0
  ! Hydrogren mass
  double precision, parameter :: Hmass = 1.67372d-24 !g

end module input_constants
