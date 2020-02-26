#!/usr/bin/env python
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
from os.path import basename
import h5py, numpy as np
import sys
import astropy.units as u
from magnetizer.parameters import Parameters
import magnetizer.extra_quantities as eq
from importlib import reload

class MagnetizerRun(object):
    """
    Provides an interface to access a given Magnetizer run.

    Parameters
    ----------
    output_path : str or list
        Path to hdf5 file containing Magnetizer output data (or, optionally,
        also the input, see input_path). If a list of filenames is provided, the
        contents of all the files are concatenated on demand.
    input_path : str or list
        Path to hdf5 file containing Magnetizer input data. If a list of
        filnames is provided, the contents of all the files are used, as in the
        output_path case (N.B make sure the lists output_path and input_path
        are compatible). If None, assumes the same as output_file_path.
        Default: None
    z_tolerance : float
        Level of tolerance in redshift to be used when MagnetizerRun.get(quantity, z)
        is invoked. Default: 0.01
    verbose: bool
        Verboses use of object methods. Default: False

    Attributes
    ----------
    verbose : bool
        verbosity
    redshifts : array
        Redshift available in this run
    times : array
        Times available in this run
    parameters : parameters.Parameters
        Object containing parameters used in the present run.
    name : str
        Name of the present run.
    ngals : int
        Number of galaxies
    ngrid : int
        Number of grid points in the profiles
    nz : int
        Number of reshifts
    gal_id : array
        Identifiers for each galaxy (consistent with the the main Magnetizer
        program's gal_ids)
    gal_id_orig : array
        Galform's GalaxyID associated with each galaxy.
    sample_z : array
        Redshift at which each galaxy was selected to be included in the sample
    ivol : array
        Number associated with the file from each galaxy was read.
    """
    def __init__(self, output_path, input_path=None,
                 z_tolerance=0.01, verbose=False):

        if isinstance(output_path, str):
            output_path = [output_path]
        if isinstance(input_path, str):
            input_path = [input_path]

        # Open files for reading
        houts = [h5py.File(path,'r') for path in output_path]

        if input_path is not None:
            hins = [h5py.File(path,'r') for path in input_path]
        else:
            hins = houts

        # Construct a dictionaries of references to the hdf5 files
        # (merges all groups...)
        self._data = []
        for hout, hin in zip(houts, hins):
            if 'Output' not in hout:
                print('Problems with file', hout.filename)
                continue
            if 'Input' not in hin:
                print('Problems with file', hin.filename)
                continue
            data_dict = {x: hout['Output'][x] for x in hout['Output']}
            data_dict.update({x: hin['Input'][x] for x in hin['Input']})
            data_dict.update({x: hout['Log'][x] for x in hout['Log']})
            self._data.append(data_dict)

        self.verbose = verbose
        self._hin = hins
        self._hout = houts
        # Uses the first hdf5 as reference
        self.redshifts = self._data[0]['z'][:]
        self.times = self._data[0]['t'][:] * u.Gyr
        self.parameters = Parameters(houts[0])
        self._cache = {}
        self._galaxies_cache = {}
        self.name = self.parameters.name.decode()
        self._z_tolerance = z_tolerance

        self.ngrid, self.nz = self._data[0]['h'].shape[1:]

        self.gal_id = []
        self.gal_id_orig = []
        self.sample_z = []
        self.ivol = []
        self._completed = []

        for i, data_dict in enumerate(self._data):
            # Creates a masks for finished runs
            completed = data_dict['completed'][:]>0.
            completed = completed[:,0] # Eliminates bogus axis
            self._completed.append(completed)

            # Stores indices and ivolumes of all galaxies
            gal_id = np.arange(completed.size)[completed]
            ivol = np.ones_like(gal_id)*i # Later this will be substituted
                                          # by the actual ivolume
            self.gal_id.append(gal_id)
            self.ivol.append(ivol)

            if 'sample_z' in data_dict: # Backwards compatibility
                self.gal_id_orig.append(data_dict['IDs'][completed])
                self.sample_z.append(data_dict['sample_z'][completed])

        for name in ('gal_id','ivol'):
            # Concatenates list attributes
            attrlist = getattr(self, name)
            setattr(self,name, np.concatenate(attrlist))

        if 'sample_z' in data_dict: # Backwards compatibility
            for name in ('gal_id_orig','sample_z'):
                # Concatenates list attributes
                attrlist = getattr(self, name)
                setattr(self,name, np.concatenate(attrlist))

        # Convenience
        self.ngals = self.gal_id.size

        # Place-holders used elsewhere
        self._valid = None
        self._valid_scalar = None


    def get(self, quantity, z, position=None, binning=None, cache=True,
            cache_intermediate=True):
        """
        Loads a given quantity at specified redshift

        Parameters
        ----------
        quantity : str
            Name of the quantiy one is interested. E.g. "Bmax"
        z : float
            Redshift
        position : float, optional
            Number of half-mass radii at which the output will be returned.
        binning : BinningObject, optional
            Specifies binning to be applied. See BinningObject
        cache : bool
            Caches the requested quantity
        cache_intermediate : bool
            Caches intermediate results needed for the computation of the
            requested quantity.

        Returns
        -------
        data : array or list
            If `binning` is None, the result is an array containing the
            values of `quantity` for all (completed) galaxies. If `quantity` is
            a scalar (e.g. Mstars_disk) or if `position` is not None, the
            result will be a 1D-array. If `quantity` is a profile (e.g. Bp) and
            `position` is None, the result will be a 2D-array, containing the
            profiles of each galaxy.
            If `binning` is provided, a list containing the results for each bin
            is returned.
        """
        reload(eq) # Allows for on-the-fly changes to this module!
        iz = self._closest_redshift_idx(z)
        key = (quantity, iz)

        if cache == False:
            cache_intermediate = False

        # If this quantity at this redshift was not already loaded, load or
        # computes it (and save to the cache)
        if key in self._cache:
            result = self._cache[key]
        else:
            if quantity == 'status':
                result = status[self._completed, iz]
                raise NotImplementedError
            elif quantity in self._data[0]:
                if quantity in multiple_entry_datasets:
                    quantity_type = 'multiple'
                else:
                    if len(self._data[0][quantity].shape)==3:
                        quantity_type = 'profile'
                    else:
                        quantity_type = 'scalar'

                if 'Units' in self._data[0][quantity].attrs:
                    unit = self._data[0][quantity].attrs['Units']
                    if len(unit)==1:
                        unit = unit[0] # unpacks array..
                    unit = units_dict[unit.decode("ascii")]
                else:
                    unit = 1

                if self.verbose:
                    print('Loading {0} at z={1}'.format(quantity,
                                                        self.redshifts[iz]))
                data = []
                for i, data_dict in enumerate(self._data):
                    data.append(self._clean(data_dict[quantity],iz,i,quantity_type))
                if quantity_type == 'multiple':
                    result = np.concatenate(data, axis=1)
                else:
                    result = np.concatenate(data)

                if unit is not None:
                    result = result*unit
            else:
                if self.verbose:
                    print('Computing {0} at z={1}'.format(quantity,
                                                          self.redshifts[iz]))

                result = eq.compute_extra_quantity(quantity,self,z=z,
                                                   cache=cache_intermediate)
            if cache:
                self._cache[key] = result

        if position is None:
            # Returns the cached quantity
            return_data = result
        else:
            # Returns the cached quantity at selected radius
            rmax_rdisc = self.parameters.grid['P_RMAX_OVER_RDISK']
            target_pos = int(self.ngrid/rmax_rdisc*position)
            if isinstance(result, u.Quantity):
                return_data = result.base[:,target_pos]*result.unit
            else:
                return_data = result[:,target_pos]


        if binning is None:
            return return_data

        else:
            return_list = [None]*binning.nbins

            for i in range(binning.nbins):
                return_list[i] = return_data[binning.masks[i]]

            return return_list

    def get_galaxy(self, quantity, gal_id, ivol=0):
        """
        Loads a given quantity for a specific galaxy

        Parameters
        ----------
        quantity : str
            Name of the quantiy one is interested. E.g. "Bmax"

        gal_id : int
            ID of the galaxy one wants to select. If inspecting closer a galaxy
            previously shown with `MagnetizerRun.get()`, use
            `MagnetizerRun.gal_id` to translate the index into an ID.

        ivol : int
            Index of the output file where the galaxy is to be found (typically,
            each output file will contain an ivolume). If inspecting closer a
            galaxy previously shown with `MagnetizerRun.get()`, use
            `MagnetizerRun.ivol` to translate the index into an ivol. Default: 0

        Returns
        -------
        data : array
            An array contaning the time evolution of the specified quanatity
            (if a non-scalar quantity is selected, last axis will correspond to
            different times).
        """
        reload(eq) # Allows for on-the-fly changes to this module!

        key = (quantity, gal_id, ivol)
        data = self._data[ivol]
        # If this quantity at this redshift was not already loaded, load or
        # computes it (and save to the cache)
        if key not in self._galaxies_cache:
            if quantity in data:
                profile = len(data[quantity].shape)==3

                if 'Units' in data[quantity].attrs:
                    unit = data[quantity].attrs['Units']
                    if len(unit)==1:
                        unit = unit[0] # unpacks array..
                    unit = units_dict[unit.decode("ascii")]
                else:
                    unit = 1

                if self.verbose:
                    print('Loading {0} for galaxy {1}'.format(quantity, gal_id))

                self._galaxies_cache[key] = self._clean_gal(data[quantity],
                                                            gal_id, ivol,
                                                            quantity, profile)
                if unit is not None:
                    self._galaxies_cache[key] = self._galaxies_cache[key]*unit

            else:
                if self.verbose:
                    print('Computing {0} for galaxy {1}'.format(quantity, gal_id))


                self._galaxies_cache[key] = eq.compute_extra_quantity(quantity,
                                                self, gal_id=gal_id, ivol=ivol)

        return self._galaxies_cache[key]

    @property
    def cache_size(self):
        size = 0
        for k in self._cache:
            size += self._cache[k].nbytes

        for k in self._galaxies_cache:
            size += self._galaxies_cache[k].nbytes

        return size/1e6 * u.Mbyte

    def show_outputs(self):
        print('{0:12s}{1:13s}{2}'.format('Quantity','Units','Description'))
        print('{0:12s}{1:13s}{2}'.format('-'*9,'-'*10, '-'*35))
        for k in self._hout[0]['Output']:
            try:
                unit = self._hout[0]['Output'][k].attrs['Units'][0].decode()
            except:
                unit = ''
            try:
                description = self._hout[0]['Output'][k].attrs['Description'][0].decode()
            except:
                description = ''
            print('{0:12s}{1:13s}{2}'.format(k, unit, description))

    def reset_galaxy_cache(self, gal_id=None):
        if gal_id is None:
            self._galaxies_cache = {}
        else:
            keys = self._galaxies_cache.keys()
            for key in keys:
                if key[1] == gal_id:
                    del self._galaxies_cache[key]

    def reset_cache(self, quantity=None, iz=None):
        if quantity is None:
            self._cache = {}
        else:
            keys = self._cache.keys()
            for key in keys:
                if key[0] == quantity:
                    if iz is None or key == (quantity,iz):
                        del self._cache[key]

    def _clean(self, dataset, iz = slice(None), ivol=0, quantity_type=True,
               pre_selected_z=False):
        """
        Removes incomplete and invalid data from a dataset.

        Parameters
        ----------
        dataset :

        iz :
             (Default value = slice(None)

        Returns
        -------

        """

        if self._valid is None:
            self._valid = {}

        if ivol not in self._valid:
            self._valid[ivol]= {}

        if iz not in self._valid[ivol]:
            # Pre-loads heights, which will be used to distinguish between valid
            # and invalid outputs.
            try:
                self._valid[ivol][iz]=self._data[ivol]['h'][self._completed[ivol],:,iz]>0
            except OSError:
                print('Problem with ivol:', ivol,'iz:',iz)
                raise

        valid = self._valid[ivol][iz]
        completed = self._completed[ivol]

        if quantity_type =='profile':
            if pre_selected_z:
                # If the redshift was previously selected
                # and if corresponds to a profile
                return np.where(valid[:,:], dataset[completed,:], np.NaN)
            return np.where(valid[:,:], dataset[completed,:,iz], np.NaN)
        elif quantity_type == 'multiple':
            if pre_selected_z:
                output = [ np.where(valid[:,0], dataset[completed,iq], np.NaN)
                          for iq in range(dataset.shape[-1])]
            else:
                output = [ np.where(valid[:,0], dataset[completed,iq,iz], np.NaN)
                          for iq in range(dataset.shape[-2])]
            if len(output)==1:
                output = output[0]
            else:
                output = np.array(output)
            # These kind of datasets can also have some extra invalid entries,
            # due to, e.g., the disk size being to small for a meaningful
            # line-of-sight integration. The following takes care of this
            return np.where( output>-1e10, output, np.NaN)
        else:
            # If doesn't correspond to a profile
            if pre_selected_z:
                # If the redshift was previously selected
                # and if doesn't correspond to a profile
                return np.where(valid[:,0], dataset[completed], np.NaN)

            return np.where(valid[:,0], dataset[completed,iz], np.NaN)


    def _clean_gal(self, dataset, gal_id, ivol, quantity, profile=True,
                   pre_selected_gal_id=False):
        """
        Removes incomplete and invalid data from a dataset.

        Parameters
        ----------
        dataset :

        gal_id :

        ivol :

        profile :
             (Default value = True)
        pre_selected_gal_id :
             (Default value = False)

        Returns
        -------

        """
        if profile:
            # Uses the scaleheight as the proxy for a valid profile
            valid = self._data[ivol]['h'][gal_id, :, :] > 0
            if not pre_selected_gal_id:
                return np.where(valid, dataset[gal_id,:,:], np.NaN)
            else:
                return np.where(valid, dataset[:,:], np.NaN)
        else:
            # Uses the disk stellar mass as proxy for a valid profile
            # (with some tolerance)
            valid = self._data[ivol]['Mstars_disk'][gal_id, :] > -0.1
            if quantity not in self._hin:
                # This deals with quantities as Bmax
                valid *= self._data[ivol]['h'][gal_id, 0, :] > 0
            if not pre_selected_gal_id:
                return np.where(valid, dataset[gal_id,:], np.NaN)
            else:
                return np.where(valid, dataset, np.NaN)


    def _closest_redshift_idx(self, z):
        """
        Auxiliary function: gets index of closest redshift
        """
        iz = np.abs(self.redshifts - z).argmin()
        z_actual = self.redshifts[iz]
        if abs(z_actual-z) > self._z_tolerance:
            print(z, z_actual, z_actual-z)
            raise ValueError("Requested redshift not available (check tolerance).")

        return iz


units_dict = {
              '': 1, # For dimensionless quantities
              'Gyr' : u.Gyr,
              'Mpc^-3' : u.Mpc**-3,
              'Msun' : u.Msun,
              'Msun/yr' : u.Msun/u.yr,
              'cm^-3' : u.cm**-3,
              '1/cm^3' : u.cm**-3,
              'erg cm^-3' : u.erg*u.cm**-3,
              'km/s' : u.km/u.s,
              'km/s/kpc': u.km/u.s/u.kpc,
              'kpc' : u.kpc,
              'kpc km/s' : u.kpc*u.km/u.s,
              'microgauss' : u.microgauss,
              'pc' : u.pc,
              's' : u.s,
              'arbitrary' : 1.0,
              'radians' : u.radian,
              'rad/m^2' : u.radian/u.m/u.m,
             }

# Quantities with multiple entries per redshift
multiple_entry_datasets = {'RM','RM_LoS_y','RM_LoS_z','column_density','theta',
                           'FRB_RM','FRB_RM_LoS_y','FRB_RM_LoS_z','FRB_column_density','FRB_theta',}
