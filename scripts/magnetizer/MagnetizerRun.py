#!/usr/bin/env python
from os.path import basename
import h5py, numpy as np
import sys
from astropy.units import Quantity
from parameters import Parameters
from extra_quantities import compute_extra_quantity, units_dict

class MagnetizerRun(object):
    """
    Provides an interface to access a given Magnetizer run.

    output_file_path -> path to hdf5 file containing Magnetizer output data.
    input_file_path -> path hdf5 file containing Magnetizer input data.
                       If absent, will assume the same as output_file_path.
    """
    def __init__(self, output_file_path, input_file_path=None,
                 z_tolerance=0.01, verbose=False):

        # Open files for reading
        hout = h5py.File(output_file_path,'r')
        if input_file_path is not None:
            hin = h5py.File(input_file_path,'r')
        else:
            hin = hout

        # Construct a dictionary of references to the hdf5 files
        # (merges all groups...)
        data_dict = {x: hout['Output'][x] for x in hout['Output']}
        data_dict.update({x: hin['Input'][x]
                              for x in hin['Input']})
        data_dict.update({x: hout['Log'][x]
                              for x in hout['Log']})
        self.verbose = verbose
        self._data = data_dict
        self._hin = hin
        self._hout = hout
        self.redshifts = self._data['z'][:]
        self.times = self._data['t'][:]
        self.parameters = Parameters(hout)
        self._cache = {}
        self.name = self.parameters.name
        self._z_tolerance = z_tolerance

        self.ngals, self.ngrid, self.nz = data_dict['h'].shape

        # Creates a mask for finished runs
        self._completed = self._data['completed'][:]>0.
        self._completed = self._completed[:,0] # Eliminates bogus axis

        self.gal_id = np.arange(self._completed.size)[self._completed]

        # Place-holders used elsewhere
        self._valid = None
        self._valid_scalar = None


    def get(self, quantity, z, position=None, binning=None):
        iz = self.__closest_redshift_idx(z)
        key = (quantity, iz)

        # If this quantity at this redshift was not already loaded, load or
        # computes it (and save to the cache)
        if key not in self._cache:
            if quantity == 'status':
                self._cache[key] = status[self._completed, iz]
                raise NotImplementedError
            elif quantity in self._data:
                profile = len(self._data[quantity].shape)==3

                if 'Units' in self._data[quantity].attrs:
                    unit = self._data[quantity].attrs['Units']
                    if len(unit)==1:
                        unit = unit[0] # unpacks array..
                    unit = units_dict[unit]
                else:
                    unit = 1

                if self.verbose:
                    print 'Loading {0} at z={1}'.format(quantity,
                                                        self.redshifts[iz])

                self._cache[key] = self._clean(self._data[quantity],iz,
                                                   profile)
            else:
                if self.verbose:
                    print 'Computing {0} at z={1}'.format(quantity,
                                                          self.redshifts[iz])
                new_data, unit = compute_extra_quantity(quantity, self._data,
                                                        select_z=iz,
                                                        return_units=True)

                profile = len(new_data.shape)==2

                self._cache[key] = self._clean(new_data, iz, profile,
                                                   pre_selected_z=True)

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

    def get_galaxy(self, quantity, gal_id):

        key = ('gal', quantity, gal_id)

        # If this quantity at this redshift was not already loaded, load or
        # computes it (and save to the cache)
        if key not in self._cache:
            if quantity in self._data:
                profile = len(self._data[quantity].shape)==3

                if 'Units' in self._data[quantity].attrs:
                    unit = self._data[quantity].attrs['Units']
                    if len(unit)==1:
                        unit = unit[0] # unpacks array..
                    unit = units_dict[unit]
                else:
                    unit = 1

                if self.verbose:
                    print 'Loading {0} for galaxy {1}'.format(quantity, gal_id)

                self._cache[key] = self._clean_gal(self._data[quantity],gal_id,
                                                   profile)
            else:
                if self.verbose:
                    print 'Computing {0} for galaxy {1}'.format(quantity, gal_id)
                new_data, unit = compute_extra_quantity(quantity, self._data,
                                                        select_gal=gal_id,
                                                        return_units=True)

                profile = len(new_data.shape)==2

                self._cache[key] = self._clean_gal(new_data, gal_id, profile,
                                                   pre_selected_gal_id=True)

            if unit is not None:
                self._cache[key] = self._cache[key]*unit

        return self._cache[key]


    def _clean(self, dataset, iz = slice(None, None, None), profile=True,
               pre_selected_z=False):
        """
        Removes incomplete and invalid data from a dataset.
        """

        if self._valid is None:
            # Pre-loads heights, which will be used to distinguish between valid
            # and invalid outputs.
            self._valid = self._data['h'][self._completed, :, :] > 0


        if profile:
            if pre_selected_z:
                # If the redshift was previously selected
                # and if corresponds to a profile
                return np.where(self._valid[:,:,iz],
                                dataset[self._completed,:],
                                np.NaN)

            return np.where(self._valid[:,:,iz],
                            dataset[self._completed,:,iz],
                            np.NaN)

        else:
            # If doesn't correspond to a profile
            if pre_selected_z:
                # If the redshift was previously selected
                # and if corresponds to a profile
                return np.where(self._valid[:,0,iz],
                                dataset[self._completed],
                                np.NaN)

            return np.where(self._valid[:,0,iz],
                                dataset[self._completed,iz],
                                np.NaN)

    def _clean_gal(self, dataset, gal_id, profile=True,
                   pre_selected_gal_id=False):
        """
        Removes incomplete and invalid data from a dataset.
        """

        if profile:
            # Uses the scaleheight as the proxy for a valid profile
            valid = self._data['h'][gal_id, :, :] > 0
            if not pre_selected_gal_id:
                return np.where(valid, dataset[gal_id,:,:], np.NaN)
            else:
                return np.where(valid, dataset[:,:], np.NaN)
        else:
            # Uses the disk stellar mass as proxy for a valid profile
            # (with some tolerance)
            valid = self._data['Mstars_disk'][gal_id, :] > -0.1
            if not pre_selected_gal_id:
                return np.where(valid, dataset[gal_id,:], np.NaN)
            else:
                return np.where(valid, dataset, np.NaN)


    def __closest_redshift_idx(self, z):
        iz = np.abs(self.redshifts - z).argmin()
        z_actual = self.redshifts[iz]
        if abs(z_actual-z) > self._z_tolerance:
            print z_actual-z
            raise ValueError("Requested redshift not available (check tolerance).")

        return iz


