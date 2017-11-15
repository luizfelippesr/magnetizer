"""
Contains routines that compute auxiliary quantities.
This shall be updated later to include more complex properties.
"""
from numpy import pi, arctan, sqrt
import numpy as np
import re
import astropy.units as u


def compute_extra_quantity(qname, mag_run, gal_id=None, z=None):
    """
    Computes extra Magnetizer quantites.

    This function needs to be called with either gal_id or z different
    from None

    Parameters
    ----------
    qname : string
        code-name of the quantity
        
    mag_run : MagnetizerRun
        
    gal_id : integer
        index of the required galaxy. Having this different from None
        signals that the redshift evolution of a given galaxy should be
        computed

    z : float
        target redshift . Having this different from None signals that
        the quantity needs to be computed at this redshift only for a
        range of galaxies.


    Returns
    -------
    astropy quantity
        if gal_id!=None, this will have shape (ngrid,nz) or (nz)
        if z!=None, this will have shape (ngals,ngrid) or (ngals)
    """

    # Makes sure the correct method is used
    if (gal_id is not None) and (z is None) :
        get = lambda quantity: mag_run.get_galaxy(quantity, gal_id)
    elif (gal_id is None) and (z is not None):
        get = lambda quantity, redshift=z: mag_run.get(quantity, redshift)
    else:
        raise ValueError, 'Must choose either gal_id or z'

    if qname == 'V':
        quantity = get('Omega')*get('r')

    elif qname == 'D':
        S = get('Shear')
        eta_t= get('etat')
        alpha = get('alp')
        h = get('h')

        quantity = alpha * S * h**3 / (eta_t**2)
        quantity = quantity.cgs

    elif qname in (r'D_{{\rm crit}}','Dc','Dcrit'):
        Ru = get('R_u')
        Cu = 0.25 # hard-coded at input_constants module
        quantity = - (pi/2.)**5 * (1. + 4./pi**2 *Cu * Ru)**2
        quantity = quantity.cgs

    elif qname == 'R_u':
        quantity = get('Uz') * get('h') / get('etat')

    elif qname == '|Bp|':
        quantity = np.abs(get('Bp'))

    elif qname == '|Br|':
        quantity = np.abs(get('Bp'))

    elif qname == '|Bz|':
        quantity = get('Bzmod')

    elif qname == 'p':
        quantity = arctan(get('Br')/get('Bp'))
        quantity = quantity.to(u.degree)

    elif qname == 'q':
        quantity = -get('Shear')/get('Omega')

    elif qname == 'l/h':
        quantity = get('l')/get('h')

    elif qname == 'h/r':
        quantity = get('h')/get('r')

    elif qname == r'\tau\Omega':
        quantity = get('tau')*get('Omega')

    elif qname == r'Btot':
        quantity = sqrt(get('Bp')**2 +
                    get('Br')**2 +
                    get('Bzmod')**2)

    elif qname == r'B_Beq':
        quantity = sqrt(get('Bp')**2 +
                    get('Br')**2 +
                    get('Bzmod')**2)
        quantity /= get('Beq')

    elif qname == r'D_Dc':
        quantity = get('D')/get('Dc')

    elif qname == r'Bfloor':
        r = get('r')
        h = get('h')
        Beq = get('Beq')
        #l = get('l') # This is not what is actually used in the code!
        l = mag_run.parameters.ISM_and_disk['P_ISM_TURBULENT_LENGTH']*u.pc
        kappa = mag_run.parameters.dynamo['P_FLOOR_KAPPA']

        Delta_r = l * kappa
        fmag = 0.5
        Ncells= np.abs(3.*r*Delta_r*h/l**3)
        brms= fmag*Beq

        quantity = np.exp(-Delta_r/2./r)*brms/Ncells**(0.5)*l/Delta_r/3.
        quantity = quantity.to(u.microgauss)

    elif qname == r'growth':
        if z is not None:
            # Redshift selection
            iz = mag_run._closest_redshift_idx(z)
            if iz==0:
                return 0.0 * get('r').base / u.Gyr

            prev_z = mag_run.redshifts[iz-1]

            Btot1 = sqrt(get('Bp')**2 +
                         get('Br')**2 +
                         get('Bzmod')**2)
            Btot0 = sqrt(get('Bp',redshift=prev_z)**2 +
                         get('Br',redshift=prev_z)**2 +
                         get('Bzmod',redshift=prev_z)**2)

            delta_lnB = np.log(Btot1/Btot0)
            delta_t = mag_run.times[iz] - mag_run.times[iz-1]

            quantity = delta_lnB/delta_t
        else:
            # Profile selection
            Btot = sqrt(get('Bp')**2 +
                         get('Br')**2 +
                         get('Bzmod')**2)

            quantity = np.empty_like(Btot.base)
            quantity[:,0] = 0.0*u.Gyr

            delta_lnB = np.log(Btot[:,1:]/Btot[:,:-1])
            delta_t = mag_run.times[1:] - mag_run.times[:-1]

            for prof in delta_lnB:
                prof /=delta_t

            quantity[:,1:] = delta_lnB

    elif qname == r'growth_max':
        if z is not None:
            # Redshift selection
            iz = mag_run._closest_redshift_idx(z)

            if iz==0:
                return np.nan * get('Mstars_disk').base / u.Gyr

            prev_z = mag_run.redshifts[iz-1]

            Bmax1 = get('Bmax')
            Bmax0 = get('Bmax',redshift=prev_z)

            delta_lnB = np.log(Bmax1/Bmax0)
            delta_t = mag_run.times[iz] - mag_run.times[iz-1]

            quantity = delta_lnB/delta_t
        else:
            # Profile selection
            Bmax = get('Bmax')

            quantity = np.empty_like(Bmax.base)/u.Gyr

            quantity[0] = 0.0

            delta_lnB = np.log(Bmax[1:]/Bmax[:-1])
            delta_t = mag_run.times[1:] - mag_run.times[:-1]

            quantity[1:] = delta_lnB/delta_t

    elif qname == 'pmax':
        quantity = get('p')
        quantity = __get_profile_max(quantity, gal_id=gal_id, z=z)

    elif qname == r'Bmax':
        quantity = get('Btot')
        quantity = __get_profile_max(quantity, gal_id=gal_id, z=z)

    elif qname == r'bmax':
        quantity = get('Beq')
        quantity = __get_profile_max(quantity, gal_id=gal_id, z=z)
        quantity *= mag_run.parameters.dynamo['FMAG']

    elif qname == r'max_r':
        raise NotImplementedError

    else:
        raise ValueError, qname + ' is unknown.'

    return quantity


def __get_profile_max(quantity, gal_id=None, z=None):

    if z is not None:
        quantity = quantity.max(axis=1)
    else:
        quantity = quantity.max(axis=0)
    return quantity
