""" 
Contains routines that compute auxiliary quantities.
This shall be updated later to include more complex properties.
"""

from numpy import pi, arctan, sqrt
import numpy as np
    
    
def compute_extra_quantity(qname, f, select_gal=np.s_[:], select_r=np.s_[:],
                           select_z=np.s_[:]):
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
    
    if qname == 'V':
        return f['Omega'][ig,ir,iz]*f['r'][ig,ir,iz]

    elif qname == 'D':
        S = f['Shear'][ig,ir,iz]
        eta_t= f['etat'][ig,ir,iz]
        alpha = f['alp'][ig,ir,iz]
        h = f['h'][ig,ir,iz]

        return alpha * S * (h/1e3)**3 / (eta_t**2)

    elif qname in (r'D_{{\rm crit}}','Dc','Dcrit'):
        Ru = compute_extra_quantity('R_u', f, ig,ir,iz)
        Cu = 0.25
        return - (pi/2.)**5 * (1. + 4./pi**2 *Cu * Ru)**2

    elif qname == 'R_u':
        return f['Uz'][ig,ir,iz] * (h/1e3) / f['etat'][ig,ir,iz]

    elif qname == 'p':
        return arctan(f['Br'][ig,ir,iz]/f['Bp'][ig,ir,iz])*180/pi

    elif qname == 'q':
        return -f['Shear'][ig,ir,iz]/f['Omega'][ig,ir,iz],

    elif qname == 'l/h':
        return f['l'][ig,ir,iz]/f['h'][ig,ir,iz]

    elif qname == 'h/r':
        return f['h'][ig,ir,iz]/f['r'][ig,ir,iz]/1e3

    elif qname == 'q/(l/h)^2':
        return -f['Shear'][ig,ir,iz]/f['Omega'][ig,ir,iz]/(f['l'][ig,ir,iz]/f['h'][ig,ir,iz])**2

    elif qname == r'\tau\Omega':
        return f['tau'][ig,ir,iz]*f['Omega'][ig,ir,iz]

    elif qname == r'Btot':
        return sqrt(f['Bp'][ig,ir,iz]**2 + f['Br'][ig,ir,iz]**2 + f['Bzmod'][ig,ir,iz]**2)

    else:
        raise ValueError, qname, 'is unknown.'


