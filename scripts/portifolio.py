""" Contains a script and functions to produce a "portifolio" of a model
    output, i.e. a file containing a series of plots showing radial profiles
    of several quantites for a sample of galaxies. """
from matplotlib.backends.backend_pdf import PdfPages
import pylab as P
import numpy as N
import h5py, random

P.rc( ('axes'),labelsize=8)
P.rc( ('xtick','ytick'),labelsize=7)
#P.rc("font",size=5)
P.rc(("legend"),fontsize=6)

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

def plot_quantity(igal, quantity, data_dict, cmap=P.cm.YlGnBu,
                  ax=P.figure().add_subplot(111),ts=None):
    """Plots the time variation of a given quantity for a given galaxy """
    if not ts:
        ts = data_dict['t'][:]
    data = data_dict[quantity]

    for t in ts:
        it = N.argmin(abs(ts-t))


        # Sets the line colour, using the colormap
        # (the factor 0.8 avoids extremely light colours)
        color = cmap((t-ts.min()*0.8)/(ts.max()-ts.min()))


        r = data_dict['r'][igal,:,it]
        ok = data_dict['Omega'][igal,:,it] > 0
        ok *= r>0
        data = data_dict[quantity][igal,:,it]
        ax.plot(r[ok], data[ok], marker='.',color=color)

        ax.set_xlabel(r'$r\,[\rm kpc]$')

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



def generate_portifolio(input_filename, selected_quantities, pdf_filename,
                        output_filename=None, galaxies_per_bin=6,
                        mass_bins = (
                                      [1e8,10**8.75],
                                      [10**8.75,10**9.5],
                                      [10**9.5, 10**10.25],
                                      [10**10.25,10**12]
                                    )):

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
            print foutput[k].attrs['Units'][0]
        if k not in selected_quantities:
            continue
        data_dict[k] = foutput[k]

    # Gets the disk stellar mass in the final redshift
    mstars = finput['Mstars_disk'][:,-1]

    selected_igals = []
    igals = N.arange(mstars.size)

    # Creates a list of galaxy indices satisfying the selection criteria
    for (m_min, m_max) in mass_bins:
        ok  = mstars > m_min
        ok *= mstars < m_max
        tmp = igals[ok]
        if tmp.size > galaxies_per_bin:
            random.shuffle(tmp)
            tmp = tmp[:galaxies_per_bin]
        selected_igals = N.append(selected_igals,tmp)
    selected_igals = selected_igals.astype(int)

    # The following was meant to be parallelized with parmap
    # however, multiprocessing breaks h5py!
    print 'Producing all figures'
    figures = [single_galaxy_portifolio(igal, data_dict, mstars=mstars) for igal in selected_igals]

    pdf = PdfPages(pdf_filename)
    for fig in figures:
        if not fig:
            continue
        pdf.savefig(fig)
    print 'Plots done. Saving file'
    pdf.close()
    return




def single_galaxy_portifolio(igal, data_dict, nrows=5, ncols=3, mstars=None):
  if not igal:
      return
  if mstars != None:
      info = r' $--$  $\log(M_{{\star,{{\rm disk}} }}/{{\rm M}}_{{\odot}}) = {0:.2f}$'.format(
        N.log10(mstars[igal]))
  else:
      infor =''
  print 'galaxy', igal
  fig = P.figure(figsize=(8.268,11.69))
  subplot_idx = 0
  for quantity in data_dict:
      if len(data_dict[quantity].shape)<3 or quantity=='r':
          continue
      subplot_idx += 1
      if subplot_idx > nrows*ncols:
          break
      ax = fig.add_subplot(nrows, ncols, subplot_idx)
      plot_quantity(igal, quantity, data_dict, ax=ax)
      fig.tight_layout()
      fig.subplots_adjust(top=0.95, bottom=0.1)
      fig.suptitle('Galaxy {0}{1}'.format(igal, info))

      norm = P.mpl.colors.Normalize(vmin=data_dict['t'][:].min(),
                                    vmax=data_dict['t'][:].max())
      ax = fig.add_axes([.065, .04, .9, .01])
      x = P.mpl.colorbar.ColorbarBase(ax, cmap=P.cm.YlGnBu, norm=norm,
                                      orientation='horizontal')
      x.set_label(r'$t\,\,[{{\rm Gyr }}]$')
      #x.set_ticks([1.0,1.5,2.0,2.5,3.0])
      #x.set_ticklabels(['$1.0$','$1.5$','$2$','$2.5$','$3.0$'])

  return fig



if __name__ == "__main__"  :
    selected_quantities = [
                            'Beq',
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
                            'tau'
                          ]
    generate_portifolio('example/example_SAM_input.hdf5',
                        selected_quantities, 'example/example_portifolio.pdf',
  output_filename='example/example_magnetized_SAM_output.hdf5')
