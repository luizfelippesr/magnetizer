#!/usr/bin/env python
from os.path import basename
import h5py, numpy as np
import sys

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
                 z_tolerance=0.01):

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

        # Place-holders used elsewhere
        self._completed = None
        self._valid = None
        self._valid_scalar = None



    def get(self, quantity, z, position=None, binning=None):
        iz = self.__closest_redshift_idx(z)
        keypair = (quantity, iz)

        # If this quantity at this redshift was not already loaded, load or
        # computes it (and save to the cache)
        if keypair not in self._cache:
            if quantity == 'status':
                self._cache[keypair] = status[self._completed, iz]

            elif quantity in self._data:
                profile = len(self._data[quantity].shape)==3

                if 'Units' in self._data[quantity].attrs:
                    unit = self._data[quantity].attrs['Units']
                    if len(unit)==1:
                        unit = unit[0] # unpacks array..
                    unit = units_dict[unit]
                else:
                    unit = 1
                self._cache[keypair] = self._clean(self._data[quantity],iz,
                                                   profile)
            else:
                new_data, unit = compute_extra_quantity(quantity, self._data,
                                                        select_z=iz,
                                                        return_units=True)

                profile = len(new_data.shape)==2

                self._cache[keypair] = self._clean(new_data, iz, profile,
                                                   pre_selected_z=True)

            if unit is not None:
                self._cache[keypair] = self._cache[keypair]*unit

        if position is None:
            # Returns the cached quantity
            return_data = self._cache[keypair]
        else:
            # Returns the cached quantity at selected radius
            rmax_rdisc = self.parameters.grid['P_RMAX_OVER_RDISK']
            target_pos = int(self.ngrid/rmax_rdisc)
            print target_pos, target_pos
            return_data = self._cache[keypair].base[:,target_pos]*self._cache[keypair].unit

        if binning is None:
            return return_data

        else:
            return_list = [None]*binning.nbins

            for i in range(binning.nbins):
                return_list[i] = return_data[binning.masks[i]]

            return return_list


    def _clean(self, dataset, iz = slice(None, None, None), profile=True,
               pre_selected_z=False):
        """
        Removes incomplete and invalid data from a dataset.
        """

        if self._completed is None:
            # Creates a mask for finished runs
            self._completed = self._data['completed'][:]>0.
            self._completed = self._completed[:,0] # Eliminates bogus axis

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


    def __closest_redshift_idx(self, z):
        iz = np.abs(self.redshifts - z).argmin()
        z_actual = self.redshifts[iz]
        if abs(z_actual-z) > self._z_tolerance:
            print z_actual-z
            raise ValueError("Requested redshift not available (check tolerance).")

        return iz


