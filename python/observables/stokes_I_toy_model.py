##### This code computes synchrotron emission for a perfectly disc galaxy with a constant value for magnetic field ####

import numpy as np
import matplotlib.pyplot as plt
from math import gamma as gamma
import math

##########Various constants appearing in the equation for emissivity #######
def emissivity_syn_constant(s):
  c_cms  =  2.99792458e10 #cm s^-1 units
  me_g   =  9.1094e-28 # weight in gram. cgs units
  e_cgs  =  4.8032e-10 # cm^3/2 g^1/2 s^-1 units
  const  =  e_cgs**3.0/(me_g*c_cms**2.0) * ( 3.0*e_cgs/( 12.5663706144 * me_g**3.0 * c_cms**5.0) )**((s-1.0)/2.0) 
  const  =  const * 3.0**0.5 / ( 12.5663706144 * (s+1.0) ) * gamma( (3.0*s-1.0)/12.0 ) * gamma( (3.0*s+19.0)/12.0 ) 
  const  = const / 1.e-23    #converting to Janski

  return const


################## Ke calculation using equipartition between magnetic field and cosmic rays #########
def Ke_equipart_cgs( Bsq_muG2,s) :
  ratio_np   = 0.01
  E1_GeV     = 5
  E2_GeV     = 100
  muG_to_G   = 1.e-6
  Ke_cgs     =  ratio_np * (s-2.0) * (muG_to_G**2.0 * Bsq_muG2/25.1327412287)/( ((0.00160218)**(2-s)) * (E1_GeV**(2.0-s) - E2_GeV**(2.0-s)) )
           #1.0e-12 -> Magnetic field from micro-gauss * Gauss.  
           #(0.00160218)**(2.0-s) -> Energy in GeV to erg
           # Ke/Bsq ~ 3.35*10**-18
  return Ke_cgs


##########Convert absolute mangitude into apparanet magnitude #######
def ext_factor_luminosity(z1):  
 omega_v    = 0.7 
 omega_m    = 1-omega_v
 h          = 0.7
 Mpc_to_cm  = 3.086e+24 # megaparsec in cm
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
 ext_factor  =  (1.0+z1) / ( 4*3.1415 * (dL_Mpc * Mpc_to_cm)**2.0 )  
         #The last (1.0+z1) factor accounts for K correction shinking of observed frequency band. 
 return ext_factor

#####################################################################
############ Provide Input Parameters Here ##########################
#####################################################################
rad_kpc = 7 #radius of galaxy in kpc
h_kpc   = 0.5  #scale height of galaxy in kpc
B_muG   = 15 # average magnetic field in micro gauss
z_ay  = np.arange(0.0016, 3.0, 0.01) #The synchrotron flux will be calculated for all these galaxies
wavelength_m   = 0.21  #wavelength in meter
s       = 3 #spectral index of cosmic ray spectrum



#constants used for conversion between units
kpc_to_cm   =  3.086e+21 # kparsec in cm
sp1b2      = (s+1.0)/2.0
sm1b2      = (s-1.0)/2.0
Bsq_muG2   = B_muG*B_muG
muG_to_G   = 1.e-6
c_cms      = 2.99792458e10 #speed of light in cm s^-1 units
vol_cm3    = 2.0*h_kpc*3.1415*rad_kpc*rad_kpc*kpc_to_cm**3.0  #volume of galaxy in cm^3 units 

#Computing the synchrotron intensity produced by the specified galaxy at various redshifts
Iz_ay = z_ay * 0.0
for (iz, z) in enumerate(z_ay):
  wavelength_cm  = wavelength_m *(1+z) * 100.0 #redshifted wavelength in cm
  nu_s           = c_cms/wavelength_cm #wavelength converted to frequency in Hz
  Ke_cm3erg2     = Ke_equipart_cgs( Bsq_muG2,s) #Ke from equipartition
  Iz_ay[iz]  = vol_cm3 * Ke_cm3erg2 * emissivity_syn_constant(s) *  (muG_to_G * B_muG/3.0)**sp1b2  * nu_s**(-sm1b2)
  Iz_ay[iz]  = Iz_ay[iz] * (1.0-1.0/math.e**4.0)/4 #this factor to account for z-dependence of B i.e. int_0^h exp(-4z/h) dz
  Iz_ay[iz]  = Iz_ay[iz] * 2/3.0 #b^2 perpenticular to LoS is taken as 2/3rd of b^2 
  Iz_ay[iz]  = Iz_ay[iz] * ext_factor_luminosity( z) #converts total flux to observed luminosity

#Plotting z vs I
plt.scatter(z_ay, Iz_ay)
plt.yscale('log')
plt.xlabel('z')
plt.ylabel('I(21cm) in Jy')
plt.show()
