#!/usr/bin/env python
""" Contains a script and functions to produce a "portfolio" of a model
    output, i.e. a file containing a series of plots showing radial profiles
    of several quantites for a sample of galaxies. """
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import h5py, random

units_dict = {}

formatted_units_dict = {'microgauss':r'\mu{{\rm G}}',
              'cm^-3':r'{{\rm cm}}^{{-3}}',
              'pc':r'{{\rm pc}}',
              'km/s':r'{{\rm km}}\,{{\rm s}}^{{-1}}',
              'km/s/kpc':r'{{\rm km}}\,{{\rm s}}^{{-1}}\,{{\rm kpc}}^{{-1}}',
              'kpc km/s':r'{{\rm kpc}}\,{{\rm km}}\,{{\rm s}}^{{-1}}'}

quantities_dict = {'Bp'   : r'B_\phi',
                   'Beq'  : r'B_{{\rm eq}}',
                   'Br'   : r'B_r',
                   'Bzmod':r'|B_z|',
                   'Uz'   : r'U_z',
                   'Shear': r'S',
                   'Omega':r'\Omega',
                   'tau'  : r'\tau',
                   'alp'  : r'\alpha',
                   'alp_k': r'\alpha_k',
                   'alp_m': r'\alpha_m',
                   'etat' : r'\eta_t'
                   }


def plot_quantity(igal, quantity, data_dict, cmap=plt.cm.YlGnBu,
                  ax=plt.figure().add_subplot(111),ts=None):
    """Plots the time variation of a given quantity for a given galaxy """
    if ts is None:
        ts = data_dict['t'][:]
    data = data_dict[quantity]

    for t in ts:
        it = np.argmin(abs(data_dict['t'][:]-t))
        # Sets the line colour, using the colormap
        # (the factor 0.8 avoids extremely light colours)
        color = cmap((t-ts.min())/(ts.max()-ts.min()))

        r = data_dict['r'][igal,:,it]
        ok = data_dict['Omega'][igal,:,it] > 0

        data = data_dict[quantity][igal,:,it]
        ax.plot(r[ok], data[ok],color=color, rasterized=True)

        ax.set_xlabel(r'$r\,[\rm kpc]$')

        if len(r[ok])==0:
            continue
        ax.set_xlim([0.0, r[ok].max()/2.0])

        # Formating gymnastics...
        if quantity in quantities_dict:
            q = quantities_dict[quantity]
        else:
            q = quantity

        if quantity in units_dict:
          q += r'\,\,['
          unit = units_dict[quantity]
          if unit in formatted_units_dict:
              q += formatted_units_dict[unit]
          else:
              q += unit
          q+=']'

        ax.set_ylabel(r'${0}$'.format(q))
        ax.grid(alpha=0.2)
        if quantity=='alp_m' or quantity=='alp':
          ax.set_ylim([-6,6])
        if quantity in ('Beq','n','h'):
          ax.set_yscale('log')


def generate_portfolio(input_filename, selected_quantities, pdf_filename,
                        output_filename=None, galaxies_per_bin=6,
                        mass_bins = (
                                      [1e8,10**8.75],
                                      [10**8.75,10**9.5],
                                      [10**9.5, 10**10.25],
                                      [10**10.25,10**12]
                                    ),
                        selected_galaxies=None
                        ):

    fi = h5py.File(input_filename,'r')
    if output_filename:
        fo = h5py.File(output_filename,'r')
    else:
        fo = fi

    finput = fi['Input']
    foutput = fo['Output']

    data_dict = {'t' : finput['t'],
                'r' : foutput['r']}
    global units_dict
    units_dict = {}
    for k in foutput:
        if 'Units' in foutput[k].attrs.keys():
            units_dict[k] = foutput[k].attrs['Units'][0]
        if k not in selected_quantities:
            continue
        data_dict[k] = foutput[k]

    # Gets the disk stellar mass in the final redshift
    mstars = finput['Mstars_disk'][:,-1]
    radius = finput['r_disk'][:,-1]
    mgas = finput['Mgas_disk'][:,-1]

    selected_igals = []
    igals = np.arange(mstars.size)

    # Creates a list of galaxy indices satisfying the selection criteria
    if selected_galaxies==None:
        for (m_min, m_max) in mass_bins:
            ok  = mstars > m_min
            ok *= mstars < m_max
            tmp = igals[ok]
            if tmp.size > galaxies_per_bin:
                random.shuffle(tmp)
                tmp = tmp[:galaxies_per_bin]
            selected_igals = np.append(selected_igals,tmp)
        selected_igals = selected_igals.astype(int)
    else:
        selected_igals = np.array(selected_galaxies)

    # The following was meant to be parallelized with parmap
    # however, multiprocessing breaks h5py!
    print 'Producing all figures'
    print
    figures = [single_galaxy_portfolio(igal, data_dict, mstars=mstars, radius=radius, mgas=mgas) for igal in selected_igals]
    pdf = PdfPages(pdf_filename)
    for fig in figures:
        if not fig:
            continue
        pdf.savefig(fig)
    print 'Plots done. Saving file'
    pdf.close()

    return


def single_galaxy_portfolio(igal, data_dict, nrows=5, ncols=3, mstars=None, radius=None, mgas=None):
  if igal==None:
      return
  if mstars is not None:
      #info = r' $-$  $\log(M_{{\star,{{\rm disk}} }}/{{\rm M}}_{{\odot}}) = {0:.2f}$'.format(
        #np.log10(mstars[igal]))
      info = r' $-$  $M_{{\star,{{\rm disk}} }}/{{\rm M}}_{{\odot}} = {0:.3g}$'.format(mstars[igal])
  else:
      info =''
  if radius is not None:
      info += r' $-$  $r_{{\rm disk}} = {0}$'.format(radius[igal])
  if mgas is not None:
      info += r' $-$  $M_{{\rm gas,disk}} = {0:.3g}$'.format(mgas[igal])

  print 'galaxy', igal, '\tmgas = 10^{0:.3}'.format(np.log10(mgas[igal])),
  print '\tmstars = 10^{0:.3}'.format(np.log10(mstars[igal]))
  fig = plt.figure(figsize=(8.268,11.69))
  subplot_idx = 0
  for quantity in data_dict:
      if len(data_dict[quantity].shape)<3 or quantity=='r':
          continue
      subplot_idx += 1
      if subplot_idx > nrows*ncols:
          break
      ax = fig.add_subplot(nrows, ncols, subplot_idx)
      plot_quantity(igal, quantity, data_dict, ax=ax, cmap=plt.cm.viridis)
      fig.tight_layout()
      fig.subplots_adjust(top=0.95, bottom=0.1)
      fig.suptitle('Galaxy {0}{1}'.format(igal, info))

      norm = matplotlib.colors.Normalize(vmin=data_dict['t'][:].min(),
                                    vmax=data_dict['t'][:].max())
      ax = fig.add_axes([.065, .04, .9, .01])
      x = matplotlib.colorbar.ColorbarBase(ax, cmap=plt.cm.viridis, norm=norm,
                                      orientation='horizontal')

      x.set_label(r'$t\,\,[{{\rm Gyr }}]$')
      #x.set_ticks([1.0,1.5,2.0,2.5,3.0])
      #x.set_ticklabels(['$1.0$','$1.5$','$2$','$2.5$','$3.0$'])

  return fig



if __name__ == "__main__"  :
    import argparse

    parser = argparse.ArgumentParser(
             description='Prepares portfolio of galaxies in a Magnetizer run, '
             'i.e. samples from stellar mass bins and plots a summary of '
             'evolution of radial profiles for a set of properties.')

    parser.add_argument("MAGNETIZER_OUTPUT",
                        help="HDF5 output of the Magnetizer run (must contain"
                        " an Input group, or a separate input file must be "
                        "specified).")

    parser.add_argument("PDF_OUTPUT",
                        help="Output filename for the PDF containing the "
                        "portfolio.")

    parser.add_argument('-mi',"--magnetizer_input", help="Name of the Magnetizer "
                        "input file (only needed if input and output hdf5 "
                        "are separate)", default=None)

    parser.add_argument('-n', '--number_of_galaxies_per_bin', default=10,
                        help='Number of galaxies to be sampled per mass bin.'
                        'Default: 10')


    parser.add_argument('-b', '--mass_bins', default='8,8.75,9.5,10.25,12',
                        help='Bin edges in logarithm of galaxy stellar mass.'
                        ' Default: 8,8.75,9.5,10.25,12')

    args = parser.parse_args()

    bins_list = np.array(args.mass_bins.split(',')).astype(float)
    mass_bins = [[10.0**m1,10.0**m2]
                  for m1,m2 in zip(bins_list[:-1],bins_list[1:])]


    pdf_filename = args.PDF_OUTPUT

    if args.magnetizer_input is None:
        input_filename = args.MAGNETIZER_OUTPUT
        output_filename = None
    else:
        input_filename = args.magnetizer_input
        output_filename = args.MAGNETIZER_OUTPUT

    plt.rc( ('axes'),labelsize=8)
    plt.rc( ('xtick','ytick'),labelsize=7)
    #plt.rc("font",size=5)
    plt.rc(("legend"),fontsize=6)
    selected_quantities = ['Beq',
                            'Bp',
                            'Br',
                            'Bzmod',
                            'Omega',
                            'Shear',
                            'Uz',
                            'alp',
                            'alp_k',
                            'alp_m',
                            'delta_r',
                            #'etat',
                            'h',
                            'l',
                            'n',
                            'tau']

    generate_portfolio(input_filename,
                       selected_quantities,
                       pdf_filename,
                       output_filename=output_filename,
                       galaxies_per_bin=int(args.number_of_galaxies_per_bin),
                       mass_bins=mass_bins,
                        selected_galaxies=None
                        )
