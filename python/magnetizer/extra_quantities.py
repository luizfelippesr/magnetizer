"""
Contains routines that compute auxiliary quantities.
This shall be updated later to include more complex properties.
"""
from numpy import pi, arctan, sqrt
import numpy as np
import re
import astropy.units as u

units_dict = {
              'Gyr' : u.Gyr,
              'Mpc^-3' : u.Mpc**-3,
              'Msun' : u.Msun,
              'Msun/yr' : u.Msun/u.yr,
              'cm^-3' : u.cm**-3,
              'erg cm^-3' : u.erg*u.cm**-3,
              'km/s' : u.km/u.s,
              'km/s/kpc': u.km/u.s/u.kpc,
              'kpc' : u.kpc,
              'kpc km/s' : u.kpc*u.km/u.s,
              'microgauss' : u.microgauss,
              'pc' : u.pc,
              's' : u.s
             }


def __compute_extra_quantity_max(qname, f, ig=slice(None), iz=slice(None)):

        if qname in f:
            quantity = f[qname][ig,:,iz]
            unit = ''
        else:
            quantity, unit = compute_extra_quantity(qname, f,
                                                    ig, slice(None), iz,
                                                    return_units=True)
        if ig==slice(None):
            quantity = quantity.max(axis=1)
        else:
            quantity = quantity.max(axis=0)

        return quantity, unit


def compute_extra_quantity(qname, f, select_gal=slice(None),
                           select_r=slice(None), select_z=slice(None),
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
    elif qname[:2] == 'V_':
        quantity = f['Omega_'+qname[2]][ig,ir,iz]*f['r'][ig,ir,iz]
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
        quantity = f['Uz'][ig,ir,iz] * ( f['h'][ig,ir,iz]/1e3) / f['etat'][ig,ir,iz]

    elif qname == '|Bp|':
        quantity = np.abs(f['Bp'][ig,ir,iz])
        unit = units_dict['microgauss']

    elif qname == '|Br|':
        quantity = np.abs(f['Bp'][ig,ir,iz])
        unit = units_dict['microgauss']

    elif qname == '|Bz|':
        quantity = f['Bzmod'][ig,ir,iz]
        unit = units_dict['microgauss']

    elif qname == 'p':
        quantity = arctan(f['Br'][ig,ir,iz]/f['Bp'][ig,ir,iz])*180/pi

    elif qname == 'q':
        quantity = -f['Shear'][ig,ir,iz]/f['Omega'][ig,ir,iz],

    elif qname == 'l/h':
        quantity = f['l'][ig,ir,iz]/f['h'][ig,ir,iz]

    elif qname == 'h/r':
        quantity = f['h'][ig,ir,iz]/f['r'][ig,ir,iz]/1e3

    elif qname == r'\tau\Omega':
        quantity = f['tau'][ig,ir,iz]*f['Omega'][ig,ir,iz]

    elif qname == r'Btot':
        quantity = sqrt(f['Bp'][ig,ir,iz]**2 +
                    f['Br'][ig,ir,iz]**2 +
                    f['Bzmod'][ig,ir,iz]**2)
        unit = units_dict['microgauss']

    elif qname == r'B_Beq':
        quantity = sqrt(f['Bp'][ig,ir,iz]**2 +
                    f['Br'][ig,ir,iz]**2 +
                    f['Bzmod'][ig,ir,iz]**2)
        quantity /= f['Beq'][ig,ir,iz]
        unit = units_dict['microgauss']

    elif qname == r'D_Dc':
        D = compute_extra_quantity('D', f, ig,ir,iz)
        Dc = compute_extra_quantity('Dcrit', f, ig,ir,iz)
        quantity = D/Dc

    elif qname == r'D_Dc_max':
        D_Dc = compute_extra_quantity('D_Dc', f, ig,slice(None),iz)

        if select_gal==slice(None,None,None):
            quantity = D_Dc.max(axis=1)
        else:
            quantity = D_Dc.max(axis=0)

    elif qname in ('active_dynamo', 'active'):
        D_Dc = compute_extra_quantity('D_Dc_max', f, ig,slice(None),iz)

        quantity = D_Dc > 1

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
            unit = 1./u.Gyr
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
        unit = 1./u.Gyr


    elif qname == r'growth_max':
        quantity, unit = __compute_extra_quantity_max('growth', f, ig, iz)


    elif qname == 'pmax':
        quantity, unit = __compute_extra_quantity_max('p', f, ig, iz)

    elif qname == r'Bmax':
        quantity, unit = __compute_extra_quantity_max('Btot', f, ig, iz)


    elif qname == r'bmax':
        quantity, unit = __compute_extra_quantity_max('Beq', f, ig, iz)
        quantity /= 2.
        unit = units_dict['microgauss']


    elif qname == r'rmax':
        btot = compute_extra_quantity('Btot', f, ig,slice(None),iz)
        r = f['r'][ig,:,iz]

        if select_gal==slice(None):
            ok = np.argmax(btot,axis=1)
        else:
            ok = np.argmax(btot,axis=0)

        quantity = r[ok]

        unit = units_dict['kpc']

    else:
        raise ValueError, qname + ' is unknown.'

    if not return_units:
        return quantity
    else:
        return quantity, unit

