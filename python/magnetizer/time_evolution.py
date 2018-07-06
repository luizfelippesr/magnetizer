# Copyright (C) 2018  Luiz Felippe S. Rodrigues, Luke Chamandy
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
import magnetizer
import numpy as np
import astropy.units as u

class TimeEvolution(object):
    """
    Stores the redshift/time evolution of 3 percentiles of a
    given quantity.

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
    cache_intermediate : bool
        Caches intermediate quantities used in computation of quantity.
        Default: True
    percentiles : list
        The three percentiles to be computed. Default: [15,50,85]

    Attributes
    ----------
    zs : array
        Redshifts
    times : astropy.units.Quantity
        Times
    med : array
        50th percentile (or the value at percentiles[1] if set)
    lower : array
        15th percentile (or the value at percentiles[0] if set)
    upper : array
        85th percentile (or the value at percentiles[2] if set)
    percentiles : list
        Available percentiles (set with the corresponding parameter)
    bins : BinningObject or list
        Link to binning used in the computation
    nbins : int
        Number of bins
    ngals : array
        Number of galaxies used for the calculation at each time/redshift

    """
    def __init__(self, quantity, mag_run, position=None, bin_objs=None,
                 target_redshifts=None, minimum_number_per_bin=5,
                 use_t=False, log0=None,
                 ignore_zero_valued=False, log=True, percentiles=[15,50,85],
                 cache_intermediate=True, cache=True):
        self.quantity = quantity
        self.run = mag_run
        self.percentiles = percentiles
        single_binning = no_binning = False

        if bin_objs is None:
            self.bins = None
            self.nbins = 1
            bin_dict = None
            no_binning = True
        else:
            self.bins = bin_objs
            if isinstance(bin_objs, magnetizer.BinningObject):
                self.nbins = bin_objs.nbins
                single_binning = True
            else:
                self.nbins = bin_objs[0].nbins
                for bin_obj in bin_objs:
                    assert bin_obj.nbins == self.nbins
                bin_dict = {bin_obj.redshift: bin_obj for bin_obj in bin_objs}

        if (no_binning or single_binning) and (target_redshifts is None):
            raise ValueError, 'Must specify either a list of binning objects ' \
              '(bin_objs=[bin_obj1,...]) or a list of redshifts (target_redshifts=[z1,...])'

        if target_redshifts is not None:
            self.zs = mag_run.redshifts[closest_indices(mag_run.redshifts,
                                                        target_redshifts)]
        else:
            self.zs = np.array(sorted(bin_dict.keys()))

        self.lower = np.empty((self.nbins, self.zs.size))*np.nan
        self.upper = np.empty((self.nbins, self.zs.size))*np.nan
        self.med = np.empty((self.nbins, self.zs.size))*np.nan
        self.ngals = np.empty((self.nbins, self.zs.size))*np.nan

        self.zs = np.array(self.zs)

        idx = closest_indices(mag_run.redshifts,self.zs)
        self.times = mag_run.times[idx]

        if use_t:
            zs_or_ts = self.times/u.Gyr
        else:
            zs_or_ts = self.zs

        self.unit = ''
        for j, z in enumerate(self.zs):
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

            for i in range(self.nbins):
                datum = zdata[i]
                # Checks whether it has units and strips them away
                self.unit, datum = get_formated_units(datum, return_base=True,
                                                      clean=log)
                datum = datum[np.isfinite(datum)]

                if ignore_zero_valued:
                    datum = datum[np.abs(datum)>zero_value_tol]

                if datum.size<minimum_number_per_bin:
                    continue

                self.lower[i][j], self.med[i][j], self.upper[i][j] = \
                                            np.percentile(datum, percentiles)
                if log:
                    for p in (self.lower, self.med, self.upper):
                        p[i][j] = np.log10(p[i][j])

                self.ngals[i][j] = datum.size

        if log and (log0 is not None):
            invalid = ~np.isfinite(self.lower)
            self.lower[invalid] = log0

def get_formated_units(quantity, return_base=False, clean=False):
    if isinstance(quantity, u.Quantity):
        unit = r'{0}'.format(quantity.unit._repr_latex_())
        unit = unit.replace('$','')
        if not clean:
            if unit is not '':
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

def closest_indices(zs, zs_target):
    izs =[]
    for zt in zip(zs_target):
        izs.append(np.abs(zs - zt).argmin())
    return izs
