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

default_quantities = [
                        'Omega',
                        'Shear',
                        'Uz',
                        'alp',
                        #'etat',
                        'h',
                        'n',
                        'Beq',
                        'Bp',
                        'Br',
                        'Bzmod',
                        'Btot',
                        'r_disk',
                        ##'Bfloor',
                        ##'growth',
                        #'Dcrit',
                        #'D_Dc',
                        #'P',
                        ##'P2',
                        ]

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
                   'Mstars_bulge': r'M_{{\star,\rm bulge}}',
                   'Mstars_disk': r'M_{{\star,\rm disc}}',
                   'Mgas_disk': r'M_{{\rm gas,disc}}',
                   'r_disk': r'r_{{\rm disc}}',
                   'Bmax': r'B_{{\rm max}}'
                   }

log_quantities = ('Beq','n','h', 'Mstars_disk','Mstars_bulge','Mgas_disk',)


def plot_output(ts, rs, quantity, name='', cmap=plt.cm.YlGnBu,
               ax=None, **args):
    """
    Plots the time variation of a quantity for a given galaxy
    """
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

    if name in log_quantities:
        ax.set_yscale('log')
    # Some advanced formating
    if name in quantities_dict:
        name = quantities_dict[name]

    ax.set_xlim([0.0, rs.base[rs.base>0].max()])
    ax.set_xlabel(r'$r\,[\rm kpc]$')
    ax.set_ylabel(r'$ {0} {1} $'.format(name,unit))
    ax.grid(alpha=0.2)



def plot_input(ts, quantity, name='', zs=None, ax=None, **args):
    """
    Plots the time variation of a quantity
    """
    if ax is None:
        ax = plt.gcf().add_subplot(111)

    if isinstance(quantity, Quantity):
        unit = r'\;[{0}]'.format(quantity.unit._repr_latex_())
        unit = unit.replace('$','')
    else:
        unit = ''

    ax.plot(ts, quantity, **args)

    if name in log_quantities:
        ax.set_yscale('log')
    # Some advanced formating
    if name in quantities_dict:
        name = quantities_dict[name]

    ax.set_xlim([0.0, ts.max()])
    ax.set_xlabel(r'$t\,[\rm Gyr]$')
    ax.set_ylabel(r'$ {0} {1} $'.format(name,unit))
    ax.grid(alpha=0.2)


def plot_mass_summary(igal, run_obj, ax=None, **kwargs):

    if ax is None:
        ax = plt.gcf().add_subplot(111)

    for name in  ('Mstars_disk','Mstars_bulge','Mgas_disk'):
        data = run_obj.get_galaxy(name, igal)
        label = '$'+quantities_dict[name]+'$'
        plot_input(run_obj.times, data, name='M', ax=ax, label=label, **kwargs)
    ax.set_yscale('log')
    ax.legend(frameon=False, loc='lower right', fontsize=7)

def galaxy_portfolio(igal, run_obj, nrows=5, ncols=3, mass_frame=True,
                     selected_quantities=None, cmap=plt.cm.viridis):
    """
    Prepares a page o plots
    """
    if selected_quantities is None:
        # Copies (do not copy reference) the default list of quantities
        selected_quantities = list(default_quantities)
    else:
        # Ensures correct type
        selected_quantities = list(selected_quantities)

    # Reads all the selected galaxy properties (and redshift variations)
    prop_dict = {}
    for name in (selected_quantities +
                 ['Mstars_disk','Mgas_disk','Mstars_bulge','r_disk','r']):
        prop_dict[name] = run_obj.get_galaxy(name, igal)

    # Get galaxy properties closest to the present day and prepares
    # a string with the info
    info = ''
    ok = np.isfinite(prop_dict['Mstars_disk']) # Avoids NaNs
    for q in ('r_disk', 'Mstars_bulge', 'Mstars_disk','Mgas_disk'):
        val = prop_dict[q][ok][-1].base[()] # Ignores units
        if q == 'r_disk':
            val = r'{0:.2f} \,\rm kpc'.format(val)
        else:
            if val==0:
                continue
            val = r'{0:.2g}\times10^{{{1}}}\rm M_\odot'.format(
                val/10**np.floor(np.log10(val)),
                int(np.floor(np.log10(val))))
        info += r' $-$  ${0} = {1}$'.format(quantities_dict[q], val)

    # Prepares the figure
    fig = plt.figure(figsize=(8.268,11.69),dpi=300)
    subplot_idx = 0

    # Plots each of the quantities
    for quantity in selected_quantities:
        subplot_idx += 1
        if subplot_idx > nrows*ncols:
            break
        ax = fig.add_subplot(nrows, ncols, subplot_idx)
        if len(prop_dict[quantity].shape) == 2:
            plot_output(run_obj.times, prop_dict['r'],
                        prop_dict[quantity], name=quantity,ax=ax,
                        cmap=cmap, linewidth=1.5, alpha=0.5)
        else:
            print quantity, len(prop_dict[quantity].shape)
            plot_input(run_obj.times, prop_dict[quantity],
                       name=quantity,ax=ax, zs = run_obj.redshifts,
                       linewidth=1.5)
    if mass_frame:
        subplot_idx += 1
        ax = fig.add_subplot(nrows, ncols, subplot_idx)
        plot_mass_summary(igal, run_obj,ax=ax, linewidth=1.5)

    # Adds title
    fig.suptitle('Galaxy {0}{1}'.format(igal, info))
    # Finds space for the color bar
    fig.tight_layout()
    fig.subplots_adjust(top=0.95, bottom=0.15)
    # Prepares colorbar
    norm = matplotlib.colors.Normalize(vmin=run_obj.times.min(),
                                       vmax=run_obj.times.max())
    ax = fig.add_axes([.065, .04, .9, .01])
    cbar = matplotlib.colorbar.ColorbarBase(ax, cmap=plt.cm.viridis, norm=norm,
                                    orientation='horizontal')
    cbar.set_label(r'$t\,\,[{{\rm Gyr }}]$')
    cbar.ax.xaxis.set_label_position("bottom")


    ticks = run_obj.times[::3]
    tlabels = ['{0:.1f}'.format(x) for x in ticks]
    zlabels = ['{0:.1f}'.format(abs(x)) for x in run_obj.redshifts[::3]]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(tlabels)


    cax2 = ax.twiny()
    cax2.set_xlim(run_obj.times.min(),run_obj.times.max())
    cax2.set_xticks(ticks)
    cax2.set_xticklabels(zlabels)
    cax2.set_xlabel("$z$")

    return fig


def generate_portfolio(run_obj, selected_quantities=None, binning_obj=None,
                       selected_galaxies=None, galaxies_per_bin=10,
                       pdf_filename=None, return_figures=False):
    """
    Prepares portfolios of various galaxies for a given Magnetizer run.
    """

    # Creates a list of galaxy indices satisfying the selection criteria
    if selected_galaxies is None:
        # Dies if neither a binning nor a galaxy numbers were specified
        if binning_obj is None:
            raise ValueError
        # Draws random galaxies based on binning information.
        selected_galaxies = []
        for mask in binning_obj.masks:
            igals = run_obj.gal_id[mask]
            if igals.size > galaxies_per_bin:
                igals = random.sample(igals, galaxies_per_bin)
            selected_galaxies = np.append(selected_galaxies ,igals)

    print 'Producing all figures'

    if return_figures:
        figures = []
    else:
        figures = None

    if pdf_filename is not None:
        pdf = PdfPages(pdf_filename)


    for igal in selected_galaxies:
        fig = galaxy_portfolio(igal, run_obj,
                               selected_quantities=selected_quantities)

        if fig is not None:
            if pdf_filename is not None:
                pdf.savefig(fig)

            if return_figures:
                figures.append(fig)
            else:
                plt.close(fig)

    if pdf_filename is not None:
        print 'Plots done. Saving file'
        pdf.close()

    return figures



def closest_indices(zs, zs_target):
    izs =[]
    for zt in zip(zs_target):
        izs.append(np.abs(zs - zt).argmin())
    return izs

def prepare_mass_bins_list(mag_run, redshifts, **kwargs):
    mass_bins = []
    redshifts = mag_run.redshifts[closest_indices(mag_run.redshifts, redshifts)]
    for z in redshifts:
        mass_bins.append(magnetizer.MassBinningObject(mag_run, z=z, **kwargs))
    return mass_bins


def plot_redshift_evolution(quantity, mag_run,
                            target_redshifts=None, bin_objs=None,
                            minimum_number_per_bin=5,
                            log=True, color='#d95f0e', **kwargs):

    single_binning = no_binning = False

    if bin_objs is None:
        nbins = 1
        bin_dict = None
        no_binning = True
    else:
        if isinstance(bin_objs, magnetizer.BinningObject):
            nbins = bin_objs.nbins
            single_binning = True
        else:
            nbins = bin_objs[0].nbins
            for bin_obj in bin_objs:
                assert bin_obj.nbins == nbins
            bin_dict = {bin_obj.redshift: bin_obj for bin_obj in bin_objs}

    if (no_binning or single_binning) and (target_redshifts is None):
        raise ValueError, 'Must specify either a list of binning objects ' \
          '(bin_objs=[bin_obj1,...]) or a list of redshifts (target_redshifts=[z1,...])'

    if target_redshifts is not None:
        zs = mag_run.redshifts[closest_indices(mag_run.redshifts, target_redshifts)]
    else:
        zs = sorted(bin_dict.keys())

    p15 = [np.empty_like(zs)*np.nan] * nbins
    p85 = [np.empty_like(zs)*np.nan] * nbins
    p50 = [np.empty_like(zs)*np.nan] * nbins

    for j, z in enumerate(zs):
        if no_binning:
            zdata = [mag_run.get(quantity, z=z),]
        else:
            if single_binning:
                zdata = mag_run.get(quantity, z=z, binning=bin_objs)
            else:
                zdata = mag_run.get(quantity, z=z, binning=bin_dict[z])

        for i in range(nbins):
            datum = zdata[i].base
            datum = datum[np.isfinite(datum)]
            if datum.size<minimum_number_per_bin:
                continue

            p15[i][j], p50[i][j], p85[i][j] = np.percentile(datum, [15,50,85])

    for i in range(nbins):
        if (nbins==1):
            plt.subplot(1,1,1)
        else:
            plt.subplot(math.ceil(float(nbins)/2),2,i+1)

        if not no_binning:
            if single_binning:
                bins = bin_objs.bins
            else:
                bins = bin_dict[zs[0]].bins

            plt.title('${0}< \log(M/M_\odot) <{1}$'.format(
                np.log10(bins[i][0].base),
                np.log10(bins[i][1].base)))

        if log:
            for p in (p15, p50, p85):
                p[i] = np.log10(p[i])

        plt.plot(zs, p50[i], color=color, **kwargs)
        plt.plot(zs, p15[i], color=color, linestyle=':')
        plt.plot(zs, p85[i], color=color, linestyle=':')
        plt.fill_between(zs, p15[i], p85[i], color=color, alpha=0.1)

        if quantity in quantities_dict:
            quantity = '$'+quantities_dict[quantity]+'$'
        plt.ylabel(quantity)
        plt.xlabel('$z$')
