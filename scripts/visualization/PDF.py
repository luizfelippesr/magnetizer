import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats as stat
from quantities_dict import quantities_dict
from cycler import cycler

plt.rc( ('axes'),labelsize=8.5)
plt.rc( ('xtick','ytick'),labelsize=8)
plt.rc("text",usetex=True)
plt.rc(("legend"),fontsize=8)

plt.rc( ('axes'),labelsize=11.5)
plt.rc( ('xtick','ytick'),labelsize=11)
plt.rc("text",usetex=True)
plt.rc(("legend"),fontsize=11)
plt.rcParams['axes.prop_cycle'] = cycler('color',['#1f78b4','#a6cee3','#33a02c','#b2df8a',                                       '#e31a1c','#fb9a99','#ff7f00','#fdbf6f',
                                      '#6a3d9a','#cab2d6'])

plt.rcParams['lines.linewidth'] = 1.75



def prepare_PDF(quantity, h5_output, h5_input=None, redshift=0,
                vmin=None, vmax=None, position=0, pdf_type='normal',
                pos_type='relative', fig=plt.figure(),
                mass_bins=np.array([1e7,1e8,1e9,1e10,1e11,1e12]),

                ):
    """
    Plots the Probability Distribution Function (PDF) of a given quantity for
    different galaxy mass bins at a given redshift.
    Optionally, plots the PDF of the log of the quantity.

    Input: quantity -> string containing the th
           h5_output -> hdf5 file containing Magnetizer output data.
           h5_input -> hdf5 file containing Magnetizer input data. If absent,
                       will assume that h5_input = h5_output.
           mass_bins -> array containing the edges of the mass bins.
                        Default: [1e7,1e8,1e9,1e10,1e11,1e12]
           redshift -> redshift at which the PDF should be computed. Default: 0
           vmin, vmax -> maximum and minimum values for the quantity to be
                         considered. If None, the maximum and minium values in
                         the model are used. Default: None
           pdf_type -> 'log' computes P(log(quantity)) while
                       'normal' computes P(quantity). Default: 'normal'
           position -> 'position', in number of half-mass radii, for the
                       calculation of the PDF (for the radial dependent quantities).
           fig -> matplotlib figure which will contain the plot. TODO

    Returns: the matplotlib figure object.
    """
    if h5_input is None:
      h5_input = h5_output

    # Selects mass-related dataset
    Mb = h5_input['Input']['Mstars_bulge']
    Mg = h5_input['Input']['Mstars_disk']
    Md = h5_input['Input']['Mgas_disk']

    # Gets index of the selected redshift
    zs = h5_input['Input']['z'][:]
    iz = np.argmin(zs - redshift)

    # Selects the dataset for the chosen quantity
    if quantity in h5_input['Input']:
        data = h5_input['Input'][quantity]
        output_quantity = False
    elif quantity in h5_output['Output']:
        data = h5_output['Output'][quantity]
        ngals, ngrid, nzs = data.shape
        output_quantity = True

        if pos_type == 'relative':
            rmax_rdisc = 2.25 # TODO read this form the parameters!!!!!
            i_target = int(ngrid/rmax_rdisc*position)
        elif pos_type =='absolute':
            raise NotImplementedError
        else:
            raise ValueError

    for mmin, mmax in zip(mass_bins[:-1], mass_bins[1:]):
        # Creates filter for the current bin
        ok = Md[:,iz]>mmin
        ok *= Md[:,iz]<mmax

        # Skips if there aren't enough galaxies in the bin
        if not len(ok[ok])>1:
            continue

        # Loads the values at the specified mass bin, radius and redshift
        values = data[ok,i_target,iz]
        values = values[values>-100]

        # Sets maximum and minimum values
        if vmax is None:
            values_max = values.max()
        else:
            values_max = vmax
        if vmin is None:
            values_min = values.min()
        else:
            values_min = vmin

        if pdf_type == 'log':
            values = np.log10(values[values>-100])
        elif pdf_type !='normal':
            raise ValueError

        # Ignores mass bin if not enough galaxies with valid values
        if not len(values)>10:
            continue

        # Uses gaussian kernel density estimator to evaluate the PDF
        kernel = stat.gaussian_kde(values)
        x = np.linspace(values_min,values_max,200)
        y = kernel.evaluate(x)

        if pdf_type == 'log':  x = 10**x

        plt.plot(x,y,
                 label=r'$10^{{ {0:.2f} }}<M/M_\odot<10^{{ {1:.2f} }}$, $N={2}$'.format(np.log10(mmin),np.log10(mmax), len(values)))
        plt.title(r'$z={0:.2f}$'.format(abs(zs[iz]))
                  )

    if quantity in quantities_dict:
        name, units = quantities_dict[quantity]
    else:
        name, units = quantity, None

    if pdf_type == 'normal':
        plt.ylabel(r'$P({0})$'.format(name, units))
    else:
        plt.xscale('log')
        if units:
            plt.ylabel(r'$P(\log({0}/{1}))$'.format(name, units))
        else:
            plt.ylabel(r'$P(\log({0}))$'.format(name))
    if units:
        plt.xlabel(r'${0}\,[{1}]$'.format(name, units))
    else:
        plt.xlabel(r'${0}$'.format(name))

    plt.legend(frameon=False)


if __name__ == "__main__"  :

    h = h5py.File('/home/nlfsr/magnetizer_runs/GON9_5000_noB.hdf5','r')
    #h = h5py.File('/home/nlfsr/magnetizer_runs/GON9_5000.hdf5','r')
    #h = h5py.File('/home/nlfsr/magnetizer_runs/GON9_5000_noB_altxi.hdf5','r')

    hi = h5py.File('GON9_5000.hdf5','r')
    prepare_PDF('h', h5_output=h, h5_input=hi, pdf_type='normal',
                position=0.5, vmax=4000)
    plt.show()
