""" Contains routines that compute auxiliary quantities.
    This shall be updated later to include more complex properties."""

from numpy import pi, arctan

def compute_extra_quantity(qname, f):
    """ Includes extra quantities in the dicitonary data_dict. If present,
        gets the information from the previously open hdf5 file h5file."""


    if qname == 'V':
          return f['Omega'][...]*f['r'][...]

    elif qname == 'D':
          S = f['Shear'][:,:]
          eta_t= f['etat'][:,:]
          alpha = f['alp'][:,:]
          h = f['h'][:,:]

          return alpha * S * (h/1e3)**3 / (eta_t**2)

    elif qname in (r'D_{{\rm crit}}','Dc','Dcrit'):
          if 'R_u' in f:
              Ru = f['R_u']
          else:
              Ru = compute_extra_quantity('R_u', f)
          Cu = 0.25
          return - (pi/2.)**5 * (1. + 4./pi**2 *Cu * Ru)**2

    elif qname == 'R_u':
          h = f['h'][:,:]
          return f['Uz'][...] * (h/1e3) / f['etat'][...]

    elif qname == 'p':
          return arctan(f['Br'][...]/f['Bp'][...])*180/pi

    elif qname == 'q':
          return -f['Shear'][...]/f['Omega'][...],

    elif qname == 'l/h':
          return f['l'][...]/f['h'][...]

    elif qname == 'h/r':
          return f['h'][...]/f['r'][...]/1e3

    elif qname == 'q/(l/h)^2':
          return -f['Shear'][...]/f['Omega'][...]/(f['l'][...]/f['h'][...])**2

    elif qname == r'\tau\Omega':
          return f['tau'][...]*f['Omega'][...]

    else:
        print qname, 'is unknown.'


