#!/usr/bin/env python
""" Contains a script and functions to produce a "portfolio" of a model
    output, i.e. a file containing a series of plots showing radial profiles
    of several quantites for a sample of galaxies. """
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
from astropy.units import Quantity
import matplotlib.pyplot as plt
import numpy as np
import h5py, random
from extra_quantities import compute_extra_quantity


units_dict = {}

formatted_units_dict = {'microgauss':r'\mu{{\rm G}}',
              'cm^-3':r'{{\rm cm}}^{{-3}}',
              'pc':r'{{\rm pc}}',
              'km/s':r'{{\rm km}}\,{{\rm s}}^{{-1}}',
              'km/s/kpc':r'{{\rm km}}\,{{\rm s}}^{{-1}}\,{{\rm kpc}}^{{-1}}',
              'kpc km/s':r'{{\rm kpc}}\,{{\rm km}}\,{{\rm s}}^{{-1}}',
              'Gyr^-1' : r'\rm Gyr^{{-1}}',
              'erg cm^-3' : r'\rm erg\,cm^{{-3}}'
              }

quantities_dict = {'Bp'   : r'\overline{{B}}_\phi',
                   'Beq'  : r'\overline{{B}}_{{\rm eq}}',
                   'Br'   : r'\overline{{B}}_r',
                   'Bzmod':r'|\overline{{B}}_z|',
                   'Uz'   : r'U_z',
                   'Shear': r'S',
                   'Omega':r'\Omega',
                   'tau'  : r'\tau',
                   'alp'  : r'\alpha',
                   'alp_k': r'\alpha_k',
                   'alp_m': r'\alpha_m',
                   'etat' : r'\eta_t',
                   'Bfloor' : r'B_{\rm floor}',
                   'Btot' : r'|\mathbf{{\overline{{B}}}}|',
                   'growth' : r'\gamma',
                   'D_Dc': r'D/D_c',
                   '|Bp|'   : r'|\overline{{B}}_\phi|',
                   'Bpmod'   : r'|\overline{{B}}_\phi|',
                   }

log_quantities = ('Beq','n','h')


def plot_frame(ts, rs, quantity, name='', cmap=plt.cm.YlGnBu,
               ax=None, **args):
    """Plots the time variation of a quantity for a given galaxy """
    if ax is None:
        ax = plt.gcf().add_subplot(111)

    if isinstance(quantity, Quantity):
        unit = r'\;[{0}]'.format(quantity.unit._repr_latex_())
        unit = unit.replace('$','')
    else:
        unit = ''

    for t, r, q in zip(ts, rs.T, quantity.T):
        # Sets the line colour, using the colormap
        color = cmap((t-ts.min())/(ts.max()-ts.min()))

        ax.plot(r, q,color=color, rasterized=True, **args)

        ax.set_xlabel(r'$r\,[\rm kpc]$')

        ax.set_xlim([0.0, r.base.max()])

        # Some advanced formating
        if name in quantities_dict:
            name = quantities_dict[name]
        ax.set_ylabel(r'$ {0} {1} $'.format(name,unit))
        ax.grid(alpha=0.2)

        if name in log_quantities:
          ax.set_yscale('log')

