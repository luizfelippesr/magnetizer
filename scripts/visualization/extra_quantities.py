""" 
Contains routines that compute auxiliary quantities.
This shall be updated later to include more complex properties.
"""
from numpy import pi, arctan, sqrt
import numpy as np
    
def compute_extra_quantity(qname, f, select_gal=slice(None,None,None),
                           select_r=slice(None,None,None),
                           select_z=slice(None,None,None),
                           return_units=False):
    """ 
    Computes extra Magnetizer quantites.
    
    Input: qname -> code-name of the quantity 
           f -> object containing the quantities keys (e.g. hdf5file['Output'])
           select_gal -> optional, index or mask for the galaxy indices
           select_r -> optional, index or mask for the radius indices
           select_z -> optional, index or mask for the redshift indices
    Returns: the computed extra quantity
    """
    # Shorthands
    ig, ir, iz = select_gal, select_r, select_z
    unit = None

    if qname == 'V':
        quantity = f['Omega'][ig,ir,iz]*f['r'][ig,ir,iz]
        unit = r'km/s'

    elif qname == 'D':
        S = f['Shear'][ig,ir,iz]
        eta_t= f['etat'][ig,ir,iz]
        alpha = f['alp'][ig,ir,iz]
        h = f['h'][ig,ir,iz]

        quantity = alpha * S * (h/1e3)**3 / (eta_t**2)

    elif qname in (r'D_{{\rm crit}}','Dc','Dcrit'):
        Ru = compute_extra_quantity('R_u', f, ig,ir,iz)
        Cu = 0.25
        quantity = - (pi/2.)**5 * (1. + 4./pi**2 *Cu * Ru)**2

    elif qname == 'R_u':
        quantity = f['Uz'][ig,ir,iz] * (h/1e3) / f['etat'][ig,ir,iz]

    elif qname == '|Bp|':
        quantity = np.abs(f['Bp'][ig,ir,iz])
        unit = r'microgauss'

    elif qname == '|Br|':
        quantity = np.abs(f['Bp'][ig,ir,iz])
        unit = r'microgauss'

    elif qname == '|Bz|':
        quantity = f['Bzmod'][ig,ir,iz]
        unit = r'microgauss'

    elif qname == 'p':
        quantity = arctan(f['Br'][ig,ir,iz]/f['Bp'][ig,ir,iz])*180/pi

    elif qname == 'q':
        quantity = -f['Shear'][ig,ir,iz]/f['Omega'][ig,ir,iz],

    elif qname == 'l/h':
        quantity = f['l'][ig,ir,iz]/f['h'][ig,ir,iz]

    elif qname == 'h/r':
        quantity = f['h'][ig,ir,iz]/f['r'][ig,ir,iz]/1e3

    elif qname == 'q/(l/h)^2':
        quantity = -f['Shear'][ig,ir,iz]/f['Omega'][ig,ir,iz]/(f['l'][ig,ir,iz]/f['h'][ig,ir,iz])**2

    elif qname == r'\tau\Omega':
        quantity = f['tau'][ig,ir,iz]*f['Omega'][ig,ir,iz]

    elif qname == r'Btot':
        quantity = sqrt(f['Bp'][ig,ir,iz]**2 +
                    f['Br'][ig,ir,iz]**2 +
                    f['Bzmod'][ig,ir,iz]**2)
        unit = r'microgauss'

    elif qname == r'Bfloor':
        h = f['h'][ig,ir,iz]/1e3
        l = f['l'][ig,ir,iz]/1e3
        Delta_r = 2.
        fmag = 0.5
        r = f['r'][ig,ir,iz]
        Ncells= np.abs(3.*r*Delta_r*h/l**3)
        Beq = f['Beq'][ig,ir,iz]
        brms= fmag*Beq
        quantity = np.exp(-Delta_r/2./r)*brms/Ncells**(0.5)*l/Delta_r/3.

    elif qname == r'growth':
        if iz==0:
            quantity = 0.0 * f['r'][ig,ir,iz]
            unit = r'Gyr^-1'
        else:
            Btot1 = sqrt(f['Bp'][ig,ir,iz]**2 +
                         f['Br'][ig,ir,iz]**2 +
                         f['Bzmod'][ig,ir,iz]**2)
            Btot0 = sqrt(f['Bp'][ig,ir,iz-1]**2 +
                         f['Br'][ig,ir,iz-1]**2 +
                         f['Bzmod'][ig,ir,iz-1]**2)

            delta_logB = np.log(Btot1) - np.log(Btot0)
            delta_t = f['t'][iz] - f['t'][iz-1]

            quantity = delta_logB/delta_t
        unit = r'Gyr^-1'
    elif qname == r'Bmax':
        if len(ig)==1:
            quantity = sqrt(f['Bp'][ig,:,iz]**2 +
                            f['Br'][ig,:,iz]**2 +
                            f['Bzmod'][ig,:,iz]**2).max()
        else:
            quantity = sqrt(f['Bp'][ig,:,iz]**2 +
                            f['Br'][ig,:,iz]**2 +
                            f['Bzmod'][ig,:,iz]**2).max(axis=1)

        unit = r'microgauss'
    elif qname == r'rmax':
        print 'asdasdasdsasd'
        btot = sqrt(f['Bp'][ig,:,iz]**2 +
                    f['Br'][ig,:,iz]**2 +
                    f['Bzmod'][ig,:,iz]**2)
        if len(ig)==1:
            quantity = f['r'][ig,np.argmax(btot),iz]
            print 'TEST'
        else:
            quantity = np.empty(btot.shape[0])
            ok = np.argmax(btot,axis=1)
            print btot.shape[0]
            for i in range(btot.shape[0]):
                quantity[i] = f['r'][i,ok[i],iz]

                print 'test', f['Bp'][i,ok[i],iz]
        unit = r'kpc'
    else:
        raise ValueError, qname + ' is unknown.'

    if not return_units:
        return quantity
    else:
        return quantity, unit

