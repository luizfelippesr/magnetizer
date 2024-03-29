# Copyright (C) 2018,2019 Luiz Felippe S. Rodrigues, Luke Chamandy
#
# This file is part of Magnetizer.
#
# Magnetizer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Magnetizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Magnetizer.  If not, see <http://www.gnu.org/licenses/>.
#
""" Contains several functions to produce visualizations of Magnetizer data. """
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from astropy.units import Quantity
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import h5py, random, math
import magnetizer
from magnetizer.extra_quantities import compute_extra_quantity
from scipy.interpolate import UnivariateSpline

units_dict = {}


mnras_column_size = 3.32153
mnras_text_size = 6.97522

default_quantities = [
                        'V',
                        'Omega',
                        'Shear',
                        'Uz',
                        #'etat',
                        'h',
                        'n',
                        'alp',
                        'Br',
                        'Bp',
                        'Beq',
                        'Bzmod',
                        'Btot',
                        'D_Dc',
                        'r_disk',
                        ##'Bfloor',
                        ##'growth',
                        #'Dcrit',
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
                   'growth' : r'\Gamma',
                   'growth_max' : r'\Gamma',
                   'D_Dc': r'D/D_c',
                   '|Bp|'   : r'|\overline{{B}}_\phi|',
                   'Bpmod'   : r'|\overline{{B}}_\phi|',
                   'Mstars_bulge': r'M_{{\star,\rm bulge}}',
                   'Mstars_disk': r'M_{{\star,\rm disc}}',
                   'Mgas_disk': r'M_{{\rm gas,disc}}',
                   'r_disk': r'r_{{\rm disc}}',
                   'rmax': r'r_{{\rm max}}',
                   'Bmax': r'B_{{\rm max}}',
                   'Beavg': r'B_{{0}}',
                   'Bavg': r'B_{{s}}',
                   'Bmax_Beq': r'B_{{\rm max}}/B_{{\rm eq}}',
                   'bmax': r'b_{{\rm max}}',
                   'pmax': r'p_{{\rm max}}',
                   'p_at_Bmax': r'p\,[ ^\circ ]',
                   'B2_B2b2_avg': r'\langle \overline{B}^2/(b^2+\overline{B}^2)\rangle',
                   'n': r'n',
                   'h': r'h',
                   'PI_P': r'{\rm PI/I}',
                   }

log_quantities = ('Beq','n','h','Mstars_disk','Mstars_bulge','Mgas_disk')


def PDF(data, name='', plot_histogram=False, ax=None, vmax=None, vmin=None,
        log=False, log0=None, do_not_plot=False, use_seaborn=True, gridsize=300,
        **args):
    """
    Computes a Probability Density Distribution

    Uses scipy's kernel density estimation to compute normalized PDFs of a
    given array.

    Parameters
    ----------
    data : array
        An array or astropy.units.Quantity object containing some data. In the
        case of the latter, units are reported in the axes.
    name : str
        Name of the quantity whose PDF is being plotted
    plot_histogram : bool
        If True, plots an histogram together with the kde, for checking.
    ax : matplotlib.axes.Axes
        The axis where the figures should be plotted (typically, a panel in a
        subplot). If set to None, a new axis is generated.
    vmin, vmax : float
        The minimum and maximum values.
    log : bool
        Whether the ideia is to plot $PDF(data)$ or $PDF(log_{10}(data))$.
    log0 : float or str or None
        Sets the behaviour when using log=True and there is a 0 value in `data`.
        If None, zero values (and infinite logs) in `data` will be propagated.
        If set to 'remove', they are removed before further processing.
        If set to a real number, the 0 is substituted by the number (this is the
        most convenient for plotting).

    Returns
    -------
    x,y : numpy.ndarray
        Values used for the PDF
    """

    unit, values = get_formated_units(data, return_base=True, clean=log)
    values = values[np.isfinite(values)]


    if ax is None and (not do_not_plot):
        ax = plt.subplot(1,1,1)

    if log:
        values = np.log10(values)

        if log0 is not None:
            if log0 == 'remove':
                values = values[np.isfinite(values)]
            else:
                values[~np.isfinite(values)] = log0

    # Sets maximum and minimum values
    if vmax is None:
        values_max = values.max()
    else:
        values_max = vmax
    if vmin is None:
        values_min = values.min()
    else:
        values_min = vmin

    if use_seaborn:
        import seaborn as sns
        if plot_histogram:
            sns.distplot(values, ax=ax, kde_kws={'gridsize': gridsize}, **args)
        else:
            sns.kdeplot(values, ax=ax, gridsize=gridsize, **args)
    else:
        import scipy.stats as stat

        # Uses gaussian kernel density estimator to evaluate the PDF
        kernel = stat.gaussian_kde(values)

        x = np.linspace(values_min,values_max, gridsize)
        y = kernel.evaluate(x)

        if not do_not_plot:
            ax.plot(x,y, **args)
            if plot_histogram:
                ax.hist(values, normed=True)


    if name in quantities_dict:
        if not log:
            if unit != '':
                quantitytxt = r'${0}{1}$'.format(quantities_dict[name],unit)
            else:
                quantitytxt = r'${0}$'.format(quantities_dict[name])
            ylabel = r'${{\rm \mathcal{{P}} }}({0})$'.format(quantities_dict[name])
        else:
            quantitytxt = r'$\log({0}/{1})$'.format(quantities_dict[name],
                                                    unit)
            ylabel = r'${{\rm \mathcal{{P}} }} [ {0} ]$'.format(quantitytxt.replace(
              '$',''))
    else:
        quantitytxt = name
        ylabel = 'PDF'
    if not do_not_plot:
        plt.ylabel(ylabel)
        plt.xlabel(quantitytxt)

    if not use_seaborn:
        return x, y


def plot_output(ts, rs, quantity, name='', cmap=plt.cm.YlGnBu,
               ax=None, **args):
    """
    Plots the time variation of a quantity for a given galaxy
    """
    if ax is None:
        ax = plt.gcf().add_subplot(111)

    unit = get_formated_units(quantity)

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


def plot_input(ts, quantity, name='', zs=None, ax=None, log=False, **args):
    """
    Plots the time variation of a quantity
    """
    if ax is None:
        ax = plt.gcf().add_subplot(111)

    unit = get_formated_units(quantity)

    ax.plot(ts, quantity, **args)

    if name in log_quantities or log:
        ax.set_yscale('log')
    # Some advanced formating
    if name in quantities_dict:
        name = quantities_dict[name]

    ax.set_xlim([0.0, ts.value.max()])
    ax.set_xlabel(r'$t\,[\rm Gyr]$')
    ax.set_ylabel(r'$ {0} {1} $'.format(name,unit))
    ax.grid(alpha=0.2)


def plot_mass_summary(igal, ivol, run_obj, ax=None, **kwargs):

    if ax is None:
        ax = plt.gcf().add_subplot(111)

    for name in  ('Mstars_disk','Mstars_bulge','Mgas_disk'):
        data = run_obj.get_galaxy(name, igal, ivol)
        label = '$'+quantities_dict[name]+'$'
        plot_input(run_obj.times, data, name='M', ax=ax, label=label, **kwargs)
    ax.set_yscale('log')
    ax.legend(frameon=False, loc='lower right', fontsize=7)


def plot_B_summary(igal, ivol, run_obj, ax=None, **kwargs):

    if ax is None:
        ax = plt.gcf().add_subplot(111)

    for name in  ('Bmax','Beavg', 'Bavg'):
        data = run_obj.get_galaxy(name, igal, ivol)
        label = '$'+quantities_dict[name]+'$'
        plot_input(run_obj.times, data, name='B', ax=ax, label=label, **kwargs)
    #ax.set_yscale('log')
    ax.legend(frameon=False, loc='lower right', fontsize=7)

def galaxy_portfolio(igal, ivol, run_obj, nrows=5, ncols=3, mass_frame=True,
                     B_frame=True, selected_quantities=None,
                     cmap=plt.cm.viridis, verbose=False):
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
        prop_dict[name] = run_obj.get_galaxy(name, igal, ivol)

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
        if verbose:
            print(igal, 'Working on ', quantity)
        subplot_idx += 1
        if subplot_idx > nrows*ncols:
            break
        ax = fig.add_subplot(nrows, ncols, subplot_idx)
        if len(prop_dict[quantity].shape) == 2:
            plot_output(run_obj.times, prop_dict['r'],
                        prop_dict[quantity], name=quantity,ax=ax,
                        cmap=cmap, linewidth=1.5, alpha=0.5)
        else:
            plot_input(run_obj.times, prop_dict[quantity],
                       name=quantity,ax=ax, zs = run_obj.redshifts,
                       linewidth=1.5)
    if mass_frame:
        subplot_idx += 1
        ax = fig.add_subplot(nrows, ncols, subplot_idx)
        plot_mass_summary(igal, ivol, run_obj,ax=ax, linewidth=1.5)

    if B_frame:
        subplot_idx += 1
        ax = fig.add_subplot(nrows, ncols, subplot_idx)
        plot_B_summary(igal, ivol, run_obj,ax=ax, linewidth=1.5)

    # Adds title
    fig.suptitle('Galaxy {0},{1}{2}'.format(igal, ivol, info))
    # Finds space for the color bar
    fig.tight_layout()
    fig.subplots_adjust(top=0.95, bottom=0.15)
    # Prepares colorbar
    norm = matplotlib.colors.Normalize(vmin=run_obj.times.value.min(),
                                       vmax=run_obj.times.value.max())
    ax = fig.add_axes([.065, .04, .9, .01])
    cbar = matplotlib.colorbar.ColorbarBase(ax, cmap=plt.cm.viridis, norm=norm,
                                    orientation='horizontal')
    cbar.set_label(r'$t\,\,[{{\rm Gyr }}]$')
    cbar.ax.xaxis.set_label_position("bottom")


    ticks = run_obj.times[::3].value
    tlabels = ['{0:.1f}'.format(x) for x in ticks]
    zlabels = ['{0:.1f}'.format(abs(x)) for x in run_obj.redshifts[::3]]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(tlabels)


    cax2 = ax.twiny()
    cax2.set_xlim(run_obj.times.value.min(),run_obj.times.value.max())
    cax2.set_xticks(ticks)
    cax2.set_xticklabels(zlabels)
    cax2.set_xlabel("$z$")

    return fig


def generate_portfolio(run_obj, selected_quantities=None, binning_obj=None,
                       selected_galaxies=None, selected_ivols=None,
                       galaxies_per_bin=10, B_frame=True,
                       pdf_filename=None, return_figures=False,
                       verbose=False):
    """
    Prepares portfolios of various galaxies for a given Magnetizer run.
    """

    # Creates a list of galaxy indices satisfying the selection criteria
    if selected_galaxies is None:
        # Dies if neither a binning nor a galaxy numbers were specified
        if binning_obj is None:
            raise ValueError
        # Draws random galaxies based on binning information.
        selected_galaxies = np.array([], dtype=int) # Avoids later concatenate
        selected_ivols = np.array([], dtype=int)    # type problem
        for mask in binning_obj.masks:
            igals = run_obj.gal_id[mask]
            ivols = run_obj.ivol[mask]
            if igals.size > galaxies_per_bin:
                idx = np.random.random_integers(0,igals.size-1,galaxies_per_bin)
                igals = igals[idx]
                ivols = ivols[idx]
            selected_galaxies = np.append(selected_galaxies, igals)
            selected_ivols = np.append(selected_ivols, ivols)
    else:
        if selected_ivols is None:
            selected_ivols = [0]*len(selected_galaxies)

    print('Producing all figures')

    if return_figures:
        figures = []
    else:
        figures = None

    if pdf_filename is not None:
        pdf = PdfPages(pdf_filename)


    for igal, ivol in zip(selected_galaxies, selected_ivols):

        fig = galaxy_portfolio(igal, ivol, run_obj, B_frame=B_frame,
                               selected_quantities=selected_quantities,
                               verbose=verbose)

        if fig is not None:
            if pdf_filename is not None:
                pdf.savefig(fig)

            if return_figures:
                figures.append(fig)
            else:
                plt.close(fig)

    if pdf_filename is not None:
        print('Plots done. Saving file')
        pdf.close()

    return figures


def plot_evolution_column(t, name, keypos=None, limits=None, x_limits=None,
                          t_ticks=None, z_ticks=None, use_z=False,
                          fig=None, figsize=None, color='b',**kwargs):
    if fig is None:
        first = True
        if figsize is None:
            figsize=(mnras_text_size/4.11,mnras_text_size*0.59)
        fig, ax = plt.subplots(t.nbins, sharex=True, figsize=figsize)
    else:
        first = False
        ax = fig.axes
    if not use_z:
        zs_or_ts = t.times.value

        if x_limits is None:
            x_limits = (t.times.min()/u.Gyr, t.times.max()/u.Gyr)

        if t_ticks is None:
            t_ticks = np.linspace(t.times.max()/u.Gyr, t.times.min()/u.Gyr, 5)
        else:
            t_ticks = np.array(t_ticks)

        if z_ticks is not None:
            # Find corresponding redshifts with spline interpolation
            z_to_t_spline_converter = UnivariateSpline(t.run.redshifts[::-1],
                                                      t.run.times[::-1]/u.Gyr,
                                                      k=5, s=0)

            z_tick_pos = z_to_t_spline_converter(z_ticks)
        else:
            z_tick_pos = t_ticks
            t_to_z_spline_converter = UnivariateSpline(t.run.times/u.Gyr,
                                                      t.run.redshifts, k=5, s=0)

            z_ticks = t_to_z_spline_converter(z_tick_pos)

        tlabels = ['{0:.1f}'.format(x) for x in t_ticks]
        zlabels = ['{0:.1f}'.format(abs(x)) for x in z_ticks]
    else:
        zs_or_ts = t.zs

        if x_limits is None:
            x_limits = (t.zs.min(), t.zs.max())

        if z_ticks is None:
            z_ticks = np.linspace(t.zs.min(), t.zs.max(), 5)
        else:
            z_ticks = np.array(z_ticks)

        if t_ticks is not None:
            # Find corresponding redshifts with spline interpolation
            t_to_z_spline_converter = UnivariateSpline(t.run.times/u.Gyr,
                                                      t.run.redshifts, k=5, s=0)

            t_tick_pos = t_to_z_spline_converter(t_ticks)
        else:
            t_tick_pos = z_ticks

            z_to_t_spline_converter = UnivariateSpline(t.run.redshifts[::-1],
                                                      t.run.times[::-1]/u.Gyr,
                                                      k=5, s=0)
            t_ticks = z_to_t_spline_converter(t_tick_pos)

        tlabels = ['{0:.1f}'.format(x) for x in t_ticks]
        zlabels = ['{0:.1f}'.format(abs(x)) for x in z_ticks]


    if limits is None:
        limits = (np.nanmin(t.lower), np.nanmax(t.upper))
    try:
        bins = t.bins.bins
    except:
        bins = t.bins[0].bins

    for i, b in enumerate(bins):
        if keypos is not None:
            ax[i].text(keypos[0],keypos[1],format_log_mass(b), size=5.85)

        ax[i].plot(zs_or_ts, t.med[i], color=color,**kwargs)
        ax[i].fill_between(zs_or_ts, t.lower[i], t.upper[i],
                           color=color, alpha=0.1, linewidth=0.05)

        if first:
            ax[i].set_ylabel(name)
            ax[i].set_ylim(limits)
            ax[i].set_yticks(ax[i].get_yticks()[:-1])
    if first:
        ax[0].set_xlim(*x_limits)
        if not use_z:
            ax[0].xaxis.set_ticks(t_ticks)
            ax[0].xaxis.set_ticklabels(tlabels)
            ax[-1].set_xlabel(r"$t\;[\rm Gyr]$")
        else:
            ax[0].xaxis.set_ticks(z_ticks)
            ax[0].xaxis.set_ticklabels(zlabels)
            ax[-1].set_xlabel(r"$z$")

        ax2 = ax[0].twiny()

        if not use_z:
            ax2.set_xlim(*x_limits)
            ax2.xaxis.set_ticks(z_tick_pos)
            ax2.xaxis.set_ticklabels(zlabels)
            ax2.set_xlabel(r"$z$")
        else:
            ax2.set_xlim(*x_limits)
            ax2.xaxis.set_ticks(t_tick_pos)
            ax2.xaxis.set_ticklabels(tlabels)
            ax2.set_xlabel(r"$t\;[\rm Gyr]$")

    fig.subplots_adjust(hspace=0,
                        left=0.28,
                        right=0.96,
                        top=0.92,
                        bottom=0.1)
    return fig

def plot_evolution_row(t, name, keypos=None, limits=None, x_limits=None,
                          t_ticks=None, z_ticks=None, use_z=False,
                          fig=None, figsize=None, color='b',**kwargs):
    if fig is None:
        first = True
        if figsize is None:
            figsize=(mnras_text_size*0.59, mnras_text_size/4.11)
        fig, ax = plt.subplots(ncols=t.nbins,
                               sharey=True, sharex=True,
                               figsize=figsize)
    else:
        first = False
        ax = fig.axes
    if not use_z:
        zs_or_ts = t.times.value

        if x_limits is None:
            x_limits = (t.times.min()/u.Gyr, t.times.max()/u.Gyr)

        if t_ticks is None:
            t_ticks = np.linspace(t.times.max()/u.Gyr, t.times.min()/u.Gyr, 5)
        else:
            t_ticks = np.array(t_ticks)

        if z_ticks is not None:
            # Find corresponding redshifts with spline interpolation
            z_to_t_spline_converter = UnivariateSpline(t.run.redshifts[::-1],
                                                      t.run.times[::-1]/u.Gyr,
                                                      k=5, s=0)

            z_tick_pos = z_to_t_spline_converter(z_ticks)
        else:
            z_tick_pos = t_ticks
            t_to_z_spline_converter = UnivariateSpline(t.run.times/u.Gyr,
                                                      t.run.redshifts, k=5, s=0)

            z_ticks = t_to_z_spline_converter(z_tick_pos)

        tlabels = ['{0:.1f}'.format(x) for x in t_ticks]
        zlabels = ['{0:.1f}'.format(abs(x)) for x in z_ticks]
    else:
        zs_or_ts = t.zs

        if x_limits is None:
            x_limits = (t.zs.min(), t.zs.max())

        if z_ticks is None:
            z_ticks = np.linspace(t.zs.min(), t.zs.max(), 5)
        else:
            z_ticks = np.array(z_ticks)

        if t_ticks is not None:
            # Find corresponding redshifts with spline interpolation
            t_to_z_spline_converter = UnivariateSpline(t.run.times/u.Gyr,
                                                      t.run.redshifts, k=5, s=0)

            t_tick_pos = t_to_z_spline_converter(t_ticks)
        else:
            t_tick_pos = z_ticks

            z_to_t_spline_converter = UnivariateSpline(t.run.redshifts[::-1],
                                                      t.run.times[::-1]/u.Gyr,
                                                      k=5, s=0)
            t_ticks = z_to_t_spline_converter(t_tick_pos)

        tlabels = ['{0:.1f}'.format(x) for x in t_ticks]
        zlabels = ['{0:.1f}'.format(abs(x)) for x in z_ticks]


    if limits is None:
        limits = (np.nanmin(t.lower), np.nanmax(t.upper))
    try:
        bins = t.bins.bins
    except:
        bins = t.bins[0].bins

    for i, b in enumerate(bins):
        if keypos is not None:
            ax[i].text(keypos[0],keypos[1],format_log_mass(b), size=5.85)

        ax[i].plot(zs_or_ts, t.med[i], color=color,**kwargs)
        ax[i].fill_between(zs_or_ts, t.lower[i], t.upper[i],
                           color=color, alpha=0.1, linewidth=0.05)

        if first:
            ax[i].set_ylim(limits)
            ax[i].set_xticks(ax[i].get_xticks()[:-1])


            ax2 = ax[i].twiny()

            if not use_z:
                ax2.set_xlim(*x_limits)
                ax2.xaxis.set_ticks(z_tick_pos[:-1])
                ax2.xaxis.set_ticklabels(zlabels[:-1])
                ax2.set_xlabel(r"$z$")
            else:
                ax2.set_xlim(*x_limits)
                ax2.xaxis.set_ticks(t_tick_pos[:-1])
                ax2.xaxis.set_ticklabels(tlabels[:-1])
                ax2.set_xlabel(r"$t\;[\rm Gyr]$")



    fig.subplots_adjust(wspace=0,
                        left=0.28,
                        right=0.96,
                        top=0.92,
                        bottom=0.1)
    return fig


def plot_profile_grid(ap, keypos=None, name=None, unit=None, ylim=None, xlim=None,
                      fig=None, ax=None, color=None, log=False, **kwargs):
    nbins, nzs, nrs = ap.med.shape
    figsize = (mnras_text_size,mnras_text_size*0.7)
    if fig is None:
        fig, ax = plt.subplots(nrows=nbins, ncols=nzs, sharex=True, sharey=True, figsize=figsize)
        first = True
    else:
        first = False

    r = np.linspace(0,ap.run.parameters.grid['P_RMAX_OVER_RDISK'],ap.ngrid)
    if name is None:
        name = ap.quantity
    if unit is None:
        unit = ap.unit
    for ib in range(nbins):

        for iz in range(nzs):
            ignore_me = ((iz==nzs-1) and (ib==nbins-1))
            if keypos is not None and not ignore_me:
                mass_interval = ap.bins[0].bins[ib]
                masstxt = '${0}< \log(M_\star/M_\odot) <{1}$'.format(
                    np.log10(mass_interval[0]/u.Msun),np.log10(mass_interval[1]/u.Msun))
                masstxt += '\n'+r'$z={0:.1f}$'.format(abs(ap.zs[iz]))
                ax[ib, iz].annotate(masstxt, keypos, fontsize=7)


            if color is None:
                c = colors[ib]
            else:
                c = color
            if not ignore_me:
                ax[ib, iz].plot(r,ap.med[ib,iz,:], color=c, **kwargs)
                ax[ib, iz].fill_between(r, ap.lower[ib,iz,:], ap.upper[ib,iz,:], alpha=0.25, color=c, linewidth=0.05)
                ax[ib, iz].set_xlim(0,r.max())
            else:
                ax[ib, iz].plot(r,ap.med[ib,iz,:]*np.NaN, color=c, **kwargs)

            if ylim is not None:
                ax[ib, iz].set_ylim(ylim)
            if log and not ignore_me:
                ax[ib, iz].set_yscale('log')
            if xlim is not None:
                ax[ib, iz].set_xlim(xlim)



    ylabel = r'${}{}$'.format(name, unit)

    if first:
        for ib in range(nbins):


            ax[ib, 0].set_ylabel(ylabel)
            ax[ib, -1].set_ylabel(ylabel)

            ax[ib, -1].yaxis.set_label_position("right")
            ax[ib, -1].yaxis.tick_right()

        for iz in range(nzs):
            ax[-1,iz].set_xlabel(r'$r/r_{50}$')
            ax[0,iz].set_xlabel(r'$r/r_{50}$')
            ax[0, iz].xaxis.set_label_position("top")
            ax[0, iz].xaxis.tick_top()

        ax[0, 0].set_yticks(ax[0, 0].get_yticks()[:-1])
        ax[0, 0].set_xticks(ax[0, 0].get_xticks()[:-1])

    fig.subplots_adjust(wspace=0, hspace=0)
    return fig, ax


def plot_redshift_evolution(quantity, mag_run, position=None,
                            target_redshifts=None, bin_objs=None,
                            minimum_number_per_bin=5, keypos=None,
                            ignore_zero_valued=False, zero_value_tol = 1e-4,
                            log=True, color='#d95f0e', log0=None,
                            colors = ['#1f78b4','#33a02c','#b2df8a',
                                      '#fb9a99','#fdbf6f','#cab2d6'],
                            single_panel=False,
                            show_t_and_z=False, use_t=False,
                            limits=None,
                            cache_intermediate=True, **kwargs):
    """
    Plots the redshift evolution of the 25th, 50th and 85 percentiles of a
    given quantity. This can be a single panel or a set of panels showing
    different mass bins.

    Parameters
    ----------
    quantity : str
        Name of the quantity to be plotted (e.g. 'Br').
    mag_run : MagnetizerRun
        A MagnetizerRun object containing the run.
    position : float
        For radial dependent quantities, plot_redshift_evolution will show the
        redshift evolution at a radius r=`position`*r_disk, where r_disk is the
        half-mass radius of the disc. For scalar/global galaxy properties
        (e.g. Mgas_disk) this parameter must be set to None.
    target_redshifts: array_like
        If present, specifies which redshifts should be used for the plot.
    bin_objs : BinningObject or list
        If a BinningObject is specified, use it to construct different panels
        for each of the redshifts in the target_redshifts parameter.
        If a list of BinningObjects is specified, use the redshift associated
        with each BinningObject.
    minimum_number_per_bin : int
        Minimum number of galaxies per bin (smaller number are not plotted).
        Default: 5
    ignore_zero_valued : bool
        Ignore galaxies where the `quantity` is less than zero_value_tol.
        Default: False
    zero_value_tol : float
        See ignore_zero_valued. Default: 1e-4
    log : bool
        Plots the log of the percentiles. Default: True
    color : str
        Color of the lines and shades (same colour for all the panels).
    single_panel : bool
        Plots different bins in the same panel.
    colors : list
        List of colours to use, if in single_panel mode.
    show_t_and_z : bool
        Shows an extra x-axis with the time coordinate.
    use_t : bool
        Use time instead of redshift
    cache_intermediate : bool
        Caches intermediate quantities used in computation of quantity.
        Default: True
    """

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
        raise ValueError('Must specify either a list of binning objects '
          '(bin_objs=[bin_obj1,...]) or a list of redshifts (target_redshifts=[z1,...])')

    if target_redshifts is not None:
        zs = mag_run.redshifts[closest_indices(mag_run.redshifts, target_redshifts)]
    else:
        zs = np.array(sorted(bin_dict.keys()))

    p15 = np.empty((nbins, zs.size))*np.nan
    p85 = np.empty((nbins, zs.size))*np.nan
    p50 = np.empty((nbins, zs.size))*np.nan
    ngals = np.empty((nbins, zs.size))*np.nan

    zs = np.array(zs)

    idx = closest_indices(mag_run.redshifts,zs)
    times = mag_run.times[idx]

    if use_t:
        zs_or_ts = times/u.Gyr
    else:
        zs_or_ts = zs

    unit = ''
    for j, z in enumerate(zs):
        if no_binning:
            zdata = [mag_run.get(quantity, z=z, position=position,
                                 cache_intermediate=cache_intermediate),]
        else:
            if single_binning:
                zdata = mag_run.get(quantity, z=z, position=position,
                                    cache_intermediate=cache_intermediate,
                                    binning=bin_objs)
            else:
                zdata = mag_run.get(quantity, z=z, position=position,
                                    cache_intermediate=cache_intermediate,
                                    binning=bin_dict[z])

        for i in range(nbins):
            datum = zdata[i]
            # Checks whether it has units and strips them away
            unit, datum = get_formated_units(datum, return_base=True, clean=log)
            datum = datum[np.isfinite(datum)]

            if ignore_zero_valued:
                datum = datum[np.abs(datum)>zero_value_tol]

            if datum.size<minimum_number_per_bin:
                continue

            p15[i][j], p50[i][j], p85[i][j] = np.percentile(datum, [15,50,85])
            ngals[i][j] = datum.size



    for i in range(nbins):
        if (nbins==1 or single_panel):
            ax = plt.subplot(1,1,1)
        else:
            ax = plt.subplot(math.ceil(float(nbins)/2),2,i+1)

        for p in (p15, p50, p85):
            if log:
                p[i] = np.log10(p[i])

            if log:
                # In a plot of the log of a quantity, it is not uncommon to
                # have p15 = 0, which means log(p15)=-inf.
                # The following lines tackle this problem, allowing the user to
                # set the value of "log(0)" manually.
                if log0 is not None:
                    p15[i][~np.isfinite(p15[i])] = log0

        if not no_binning:
            if single_binning:
                bins = bin_objs.bins
            else:
                bins = bin_dict[zs[0]].bins

        if single_panel:
            color = colors[i]


        plt.plot(zs_or_ts, p50[i], color=color, **kwargs)
        plt.plot(zs_or_ts, p15[i], color=color, linestyle=':')
        plt.plot(zs_or_ts, p85[i], color=color, linestyle=':')
        plt.fill_between(zs_or_ts, p15[i], p85[i], color=color, alpha=0.1)

        if quantity in quantities_dict:
            if not log:
                if unit != '':
                    quantitytxt = r'${0}{1}$'.format(quantities_dict[quantity],unit)
                else:
                    quantitytxt = r'${0}$'.format(quantities_dict[quantity],unit)
            else:
                quantitytxt = r'$\log({0}/{1})$'.format(quantities_dict[quantity],
                                                        unit)
        else:
            quantitytxt = quantity
        if limits is not None:
            plt.axis(limits)

        plt.ylabel(quantitytxt)
        if not use_t:
            plt.xlabel('$z$')
        else:
            plt.xlabel(r"$t\;[\rm Gyr]$")

        if not (no_binning or single_panel):
            masstxt = '${0}< \log(M_\star/M_\odot) <{1}$'.format(
              np.log10(bins[i][0].base), np.log10(bins[i][1].base))
            if keypos is None:
                #ok = np.isfinite(p50[i])
                #ypos = (p50[i][ok].max()-p50[i][ok].min())*0.9+p50[i][ok].min()
                #keypos = (zs.max()*0.3, ypos)
                 ax.tick_params(right=False)
                 ax_label = ax.twinx()
                 ax_label.tick_params(right=False)
                 plt.setp(ax_label.get_yticklabels(), visible=False)
                 ax_label.set_ylabel(masstxt)
            else:
                plt.annotate(masstxt, (keypos[0], keypos[1]))

        if show_t_and_z:

            if not use_t:
                if limits is not None:
                    zmin, zmax = limits[0], limits[1]
                else:
                    zmin, zmax = zs.min(), zs.max()
                redshifts = np.linspace(zmin, zmax, 7)
                idx = closest_indices(mag_run.redshifts,redshifts)
            else:
                if limits is not None:
                    tmin, tmax = limits[0], limits[1]
                else:
                    tmin, tmax = times.min()/u.Gyr, times.max()/u.Gyr
                ltimes = np.linspace(tmin, tmax, 7)
                idx = closest_indices(mag_run.times/u.Gyr,ltimes)

            redshifts = mag_run.redshifts[idx]
            ltimes = mag_run.times[idx]/u.Gyr

            tlabels = ['{0:.1f}'.format(x) for x in ltimes]
            zlabels = ['{0:.1f}'.format(abs(x)) for x in redshifts]

            if not use_t:
                ax.set_xlim(zmin,zmax)
                ax.xaxis.set_ticks(redshifts)
                ax.xaxis.set_ticklabels(zlabels)

                ax2 = ax.twiny()
                ax2.set_xlim(zmin,zmax)
                ax2.xaxis.set_ticks(redshifts)
                ax2.xaxis.set_ticklabels(tlabels)
                ax2.set_xlabel(r"$t\;[\rm Gyr]$")
            else:
                ax.set_xlim(tmin, tmax)
                ax.xaxis.set_ticks(ltimes)
                ax.xaxis.set_ticklabels(tlabels)

                ax2 = ax.twiny()
                ax2.set_xlim(tmin, tmax)
                ax2.xaxis.set_ticks(ltimes)
                ax2.xaxis.set_ticklabels(zlabels)
                ax2.set_xlabel(r"$z$")


    return ngals


def closest_indices(zs, zs_target):
    izs =[]
    for zt in zip(zs_target):
        izs.append(np.abs(zs - zt).argmin())
    return izs


def prepare_mass_bins_list(mag_run, redshifts,
                           filter_quantity=None,
                           filter_threshold=None,
                           filter_greater=False,
                           spirals_only=False,
                           fixed_binning_redshift=None,
                           binning_obj=magnetizer.MassBinningObject,
                           extra_filter_list=None,
                           **kwargs):
    """
    Helper function that prepares a list of `MassBinningObject`s computed at
    each redshift.



    Parameters
    ----------
    mag_run


    """
    mass_bins = []
    redshifts = mag_run.redshifts[closest_indices(mag_run.redshifts, redshifts)]
    for i, z in enumerate(redshifts):
        if filter_quantity is None:
            filt = None
        else:
            if filter_threshold is None:
                raise ValueError
            if (isinstance(filter_threshold, Quantity) or
                type(filter_threshold) == type(1.0)) :
                filt = filter_threshold
            else:
                filt = filter_threshold[i]

            if not filter_greater:
                filt = mag_run.get(filter_quantity, z) > filt
            else:
                filt = mag_run.get(filter_quantity, z) < filt
        if spirals_only:
            spiral_filter = mag_run.get('Mstars_bulge',z) / (
                                                 mag_run.get('Mstars_bulge', z)
                                               + mag_run.get('Mstars_disk', z)
                                              ) < 0.5
            if filt is None:
                filt = spiral_filter
            else:
                filt = spiral_filter * filt

        if extra_filter_list is not None:
            if filt is None:
                filt = extra_filter_list[i]
            else:
                filt *= extra_filter_list[i]

        if fixed_binning_redshift is None:
            bin_obj = binning_obj(mag_run, extra_filter=filt, z=z, **kwargs)
        else:
            bin_obj = binning_obj(mag_run, extra_filter=filt,
                                  z=fixed_binning_redshift, **kwargs)
            bin_obj.redshift = z # Overwrites the redshift
        mass_bins.append(bin_obj)

    return mass_bins


def get_formated_units(quantity, return_base=False, clean=False):
    if isinstance(quantity, Quantity):
        unit = r'{0}'.format(quantity.unit._repr_latex_())
        unit = unit.replace('$','')
        if not clean:
            unit = r'\;[{0}]'.format(unit)
        if return_base:
            base = quantity.base
    else:
        unit = ''
        if return_base:
            base = quantity
    if return_base:
        return unit, base
    else:
        return unit



def format_log_mass(v):
    # Need to make this better later?
    masstxt = '${0}< \log(M_\star/M_\odot) \leq {1}$'.format(
              np.log10(v[0].base), np.log10(v[1].base))
    return masstxt
