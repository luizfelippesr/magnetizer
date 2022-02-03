##### This code computes synchrotron emission for a perfectly disc galaxy with a constant value for magnetic field ####

import numpy as np
import matplotlib.pyplot as plt
import math as ma

c_cms         = 2.99792458e10 #cm s^-1 units
me_g          = 9.1094e-28 # weight in gram. cgs units
e_cgs         = 4.8032e-10 # cm^3/2 g^1/2 s^-1 units
Jy_to_cgs     = 1e-23  #1 Jy = 10^{-23} erg s^-1 cm^-2 Hz^-1
muG_to_G      = 1.e-6
GeV_to_erg    = 0.00160218
kpc_to_cm     = 3.086e+21 # kparsec in cm
Mpc_to_cm     = 3.086e+24 # megaparsec in cm
deg_to_arcsec = 3600e0
deg_to_rad    = ma.pi/180

##########Various constants appearing in the equation for emissivity #######
def emissivity_syn_constant(s):
  const  =  e_cgs**3.0/(me_g*c_cms**2.0) * ( 3.0*e_cgs/( 4*ma.pi * me_g**3.0 * c_cms**5.0) )**((s-1.0)/2.0) 
  const  =  const * 3.0**0.5 / ( 4*ma.pi * (s+1.0) ) * ma.gamma( (3.0*s-1.0)/12.0 ) * ma.gamma( (3.0*s+19.0)/12.0 ) 
  const  = const / Jy_to_cgs    #converting to Janski

  return const


################## Ke calculation using equipartition between magnetic field and cosmic rays #########
def Ke_equipart_cgs( Bsq_muG2,s) :
  ratio_np   = 0.01
  E1_GeV     = 5
  E2_GeV     = 100
  Ke_cgs     =  ratio_np * (s-2.0) * (muG_to_G**2.0 * Bsq_muG2/8/ma.pi)/( (GeV_to_erg**(2-s)) * (E1_GeV**(2.0-s) - E2_GeV**(2.0-s)) )
           #1.0e-12 -> Magnetic field from micro-gauss * Gauss.  
           #(0.00160218)**(2.0-s) -> Energy in GeV to erg
           # Ke/Bsq ~ 3.35*10**-18
  return Ke_cgs


##########Convert absolute mangitude into apparanet magnitude #######
def ext_factor_luminosity(z1):  
 omega_v    = 0.7 
 omega_m    = 1e0-omega_v
 h          = 0.7
 zfinal     = z1 
 cbyH0_Mpch = 2997.92  
 if (z1 <= 0.002):
   zfinal = 0.002

 var1 = 0.0 ;
 z2   = 0.0
 dz2  = 0.0003

 while True: 
   z2   = z2 + dz2 ;
   var2 = 1.0 / ( omega_v + omega_m * (1.0+z2)**3.0 )**0.5
   var1 = var1 + 2.0*dz2 * var2 
   z2   = z2 + dz2 
   if (z2>zfinal):
       break
 dL_Mpc     = (1.0+z1) * (cbyH0_Mpch/h) * var1 # var1 is luminosity distance. Varified by online calculator
 ext_factor  =  (1.0+z1) / ( 4*ma.pi * (dL_Mpc * Mpc_to_cm)**2.0 )  
         #The last (1.0+z1) factor accounts for K correction shinking of observed frequency band. 
 return ext_factor

#####################################################################
############ Provide Input Parameters Here ##########################
#####################################################################
#rad_kpc = 7 #radius of galaxy in kpc
# From Fletcher+11, Table 2 I take r=4.8 kpc
rad_kpc     = 4.8 #radius of galaxy in kpc
rad_min_kpc = 0.0
h_kpc       = 0.5  #scale height of galaxy in kpc
#B_muG   = 15 # average magnetic field in micro gauss
B_muG       = 30e0 #LC: Fletcher+11 p. 2411, 1st column, mentions 20 muG. But Beck+19 Table 3 gives 16-17 in the relevant radial range
#z_ay  = np.arange(0.0016, 3.0, 0.01) #The synchrotron flux will be calculated for all these galaxies
#LC: start at z=0
z_ay        = np.arange(0.0, 3.0, 0.01) #The synchrotron flux will be calculated for all these galaxies
#wavelength_m   = 0.21  #wavelength in meter
wavelength_m   = 0.06  #wavelength in meter
s              = 3e0 #spectral index of cosmic ray spectrum



#constants used for conversion between units
sp1b2      = (s+1.0)/2.0
sm1b2      = (s-1.0)/2.0
Bsq_muG2   = B_muG*B_muG
vol_cm3    = 2.0*h_kpc*ma.pi*(rad_kpc*2.0 -rad_min_kpc**2.0) *kpc_to_cm**3.0  #volume of galaxy in cm^3 units 

#Computing the synchrotron intensity produced by the specified galaxy at various redshifts
Jy_Iz_ay = z_ay * 0.0
for (iz, z) in enumerate(z_ay):
  wavelength_cm  = wavelength_m *(1e0+z) * 100.0 #redshifted wavelength in cm
  nu_s           = c_cms/wavelength_cm #wavelength converted to frequency in Hz
  Ke_cm3erg2     = Ke_equipart_cgs( Bsq_muG2,s) #Ke from equipartition
  #Jy_Iz_ay[iz]  = vol_cm3 * Ke_cm3erg2 * emissivity_syn_constant(s) *  (muG_to_G * B_muG/3.0)**sp1b2  * nu_s**(-sm1b2)
  #LC: I removed the factor of 3 in the denominator as could not see where it comes from
  Jy_Iz_ay[iz]  = vol_cm3 * Ke_cm3erg2 * emissivity_syn_constant(s) *  (muG_to_G * B_muG)**sp1b2  * nu_s**(-sm1b2)
  #Iz_ay[iz]  = Iz_ay[iz] * (1.0-1.0/ma.e**4.0)/4 #this factor to account for z-dependence of B i.e. int_0^h exp(-4z/h) dz
  #LC: neglect the above factor factor to test z=0 case?
  #Iz_ay[iz]  = Iz_ay[iz] * 2.0/3.0 #b^2 perpenticular to LoS is taken as 2/3rd of b^2
  Jy_Iz_ay[iz]  = Jy_Iz_ay[iz] * (2.0/3.0)**(sp1b2/2) #b^2 perpenticular to LoS is taken as 2/3rd of b^2
  #LC: neglect geometry, assume face-on, and B_LOS~=0, so that B_Perp ~= B
  Jy_Iz_ay[iz]  = Jy_Iz_ay[iz] * ext_factor_luminosity(z) #converts total flux to observed luminosity

print()
print('I(z=0) [Jy] = ',Jy_Iz_ay[0])
print('ext_factor_luminosity(z=0)',ext_factor_luminosity(0))
print('-------------------')



#Distance to M51 is known which can be used directly to compute observed luminosity from total flux
Mpc_D = 7.6	#used by Fletcher+11 for distance to M51
cm_D  = Mpc_to_cm * Mpc_D
ext_factor_nocosm = 1e0/4/ma.pi/cm_D**2
print('factor neglecting cosmology',ext_factor_nocosm)
Jy_Iz_ay[0]  = Jy_Iz_ay[0]*ext_factor_nocosm/ext_factor_luminosity(0)	#take out cosmology for now


#To find what is the flux per beam
arcsec_Dbeam = 15e0
kpc_Dbeam = arcsec_Dbeam/deg_to_arcsec*deg_to_rad*cm_D/kpc_to_cm  #See caption of Table 2 of Fletcher+11
#Beam diameter: Fletcher+11, p. 2405, Eq. 3 mentions 0.6kpc, while Fletcher+11 p.2411 1st column mentions 0.15 kpc, 15 arcsec~=0.55kpc
print('Effective beam diameter [kpc] = ',kpc_Dbeam)
kpc_Rbeam = kpc_Dbeam/2
Nbeam = (rad_kpc/kpc_Rbeam)**2


print('I(z=0) per beam [Jy/beam] = ',Jy_Iz_ay[0]/Nbeam)
print('Fletcher Table 2: total lambda 6cm radio intensity, mean of arm and interarm regions between r={0} and {1} kpc ~= 8e-4 Jy/beam'.format(rad_min_kpc, rad_kpc))
print()

# #Plotting z vs I
# plt.scatter(z_ay, Jy_Iz_ay)
# plt.yscale('log')
# plt.xlabel('z')
# plt.ylabel('I(21cm) in Jy')
# plt.show()
