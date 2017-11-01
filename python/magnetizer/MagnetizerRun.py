#!/usr/bin/env python
from os.path import basename
import h5py, numpy as np
import sys
from astropy.units import Quantity
from parameters import Parameters
import extra_quantities as eq

class MagnetizerRun(object):
    """
    Provides an interface to access a given Magnetizer run.

    output_file_path -> path to hdf5 file containing Magnetizer output data.
    input_file_path -> path hdf5 file containing Magnetizer input data.
                       If absent, will assume the same as output_file_path.
    """
    def __init__(self, output_path, input_path=None,
                 z_tolerance=0.01, verbose=False):

        if isinstance(output_path, basestring):
            output_path = [output_path]
        if isinstance(input_path, basestring):
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
            data_dict = {x: hout['Output'][x] for x in hout['Output']}
            data_dict.update({x: hin['Input'][x] for x in hin['Input']})
            data_dict.update({x: hout['Log'][x] for x in hout['Log']})
            self._data.append(data_dict)

        self.verbose = verbose
        self._hin = hins
        self._hout = houts
        # Uses the first hdf5 as reference
        self.redshifts = self._data[0]['z'][:]
        self.times = self._data[0]['t'][:]
        self.parameters = Parameters(houts[0])
        self._cache = {}
        self._galaxies_cache = {}
        self.name = self.parameters.name
        self._z_tolerance = z_tolerance

        self.ngals, self.ngrid, self.nz = self._data[0]['h'].shape

        self.gal_id = []
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

        self.gal_id = np.concatenate(self.gal_id)
        self.ivol = np.concatenate(self.ivol)

        # Place-holders used elsewhere
        self._valid = None
        self._valid_scalar = None


    def get(self, quantity, z, position=None, binning=None):
        reload(eq) # Allows for on-the-fly changes to this module!
        iz = self.__closest_redshift_idx(z)
        key = (quantity, iz)

        # If this quantity at this redshift was not already loaded, load or
        # computes it (and save to the cache)
        if key not in self._cache:
            if quantity == 'status':
                self._cache[key] = status[self._completed, iz]
                raise NotImplementedError
            elif quantity in self._data[0]:
                profile = len(self._data[0][quantity].shape)==3

                if 'Units' in self._data[0][quantity].attrs:
                    unit = self._data[0][quantity].attrs['Units']
                    if len(unit)==1:
                        unit = unit[0] # unpacks array..
                    unit = eq.units_dict[unit]
                else:
                    unit = 1

                if self.verbose:
                    print 'Loading {0} at z={1}'.format(quantity,
                                                        self.redshifts[iz])
                data = []
                for i, data_dict in enumerate(self._data):
                      data.append(self._clean(data_dict[quantity],iz,i,profile))
            else:
                if self.verbose:
                    print 'Computing {0} at z={1}'.format(quantity,
                                                          self.redshifts[iz])
                data = []
                for i, data_dict in enumerate(self._data):
                    new_data, unit = eq.compute_extra_quantity(quantity,
                                                               data_dict,
                                                               select_z=iz,
                                                               return_units=True)
                    profile = len(new_data.shape)==2

                    data.append( self._clean(new_data, iz, i, profile,
                                             pre_selected_z=True) )

            self._cache[key] = np.concatenate(data)
            if unit is not None:
                self._cache[key] = self._cache[key]*unit

        if position is None:
            # Returns the cached quantity
            return_data = self._cache[key]
        else:
            # Returns the cached quantity at selected radius
            rmax_rdisc = self.parameters.grid['P_RMAX_OVER_RDISK']
            target_pos = int(self.ngrid/rmax_rdisc*position)
            if isinstance(self._cache[key], Quantity):
                return_data = self._cache[key].base[:,target_pos]*self._cache[key].unit
            else:
                return_data = self._cache[key][:,target_pos]


        if binning is None:
            return return_data

        else:
            return_list = [None]*binning.nbins

            for i in range(binning.nbins):
                return_list[i] = return_data[binning.masks[i]]

            return return_list

    def get_galaxy(self, quantity, gal_id, ivol=0):

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
                    unit = eq.units_dict[unit]
                else:
                    unit = 1

                if self.verbose:
                    print 'Loading {0} for galaxy {1}'.format(quantity, gal_id)

                self._galaxies_cache[key] = self._clean_gal(data[quantity],
                                                            gal_id,
                                                            ivol, profile)
            else:
                if self.verbose:
                    print 'Computing {0} for galaxy {1}'.format(quantity, gal_id)
                new_data, unit = eq.compute_extra_quantity(quantity, data,
                                                        select_gal=gal_id,
                                                        return_units=True)

                profile = len(new_data.shape)==2

                self._galaxies_cache[key] = self._clean_gal(new_data, gal_id,
                                        ivol, profile, pre_selected_gal_id=True)

            if unit is not None:
                self._galaxies_cache[key] = self._galaxies_cache[key]*unit

        return self._galaxies_cache[key]


    def _clean(self, dataset, iz = slice(None), ivol=0, profile=True,
               pre_selected_z=False):
        """
        Removes incomplete and invalid data from a dataset.
        """

        if self._valid is None:
            self._valid = {}

        if ivol not in self._valid:
            # Pre-loads heights, which will be used to distinguish between valid
            # and invalid outputs.
            self._valid[ivol] = self._data[ivol]['h'][self._completed[ivol], :, :] > 0

        valid = self._valid[ivol]
        completed = self._completed[ivol]

        if profile:
            if pre_selected_z:
                # If the redshift was previously selected
                # and if corresponds to a profile
                return np.where(valid[:,:,iz], dataset[completed,:], np.NaN)

            return np.where(valid[:,:,iz], dataset[completed,:,iz], np.NaN)
        else:
            # If doesn't correspond to a profile
            if pre_selected_z:
                # If the redshift was previously selected
                # and if doesn't correspond to a profile
                return np.where(valid[:,0,iz], dataset[completed], np.NaN)

            return np.where(valid[:,0,iz], dataset[completed,iz], np.NaN)


    def _clean_gal(self, dataset, gal_id, ivol, profile=True,
                   pre_selected_gal_id=False):
        """
        Removes incomplete and invalid data from a dataset.
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
            if not pre_selected_gal_id:
                return np.where(valid, dataset[gal_id,:], np.NaN)
            else:
                # This deals with quantities as Bmax
                valid *= self._data[ivol]['h'][gal_id, 0, :] > 0
                return np.where(valid, dataset, np.NaN)


    def __closest_redshift_idx(self, z):
        iz = np.abs(self.redshifts - z).argmin()
        z_actual = self.redshifts[iz]
        if abs(z_actual-z) > self._z_tolerance:
            print z, z_actual, z_actual-z
            raise ValueError("Requested redshift not available (check tolerance).")

        return iz


