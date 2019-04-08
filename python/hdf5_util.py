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
import h5py
import numpy as N

def add_dataset(storage, dataset_name, dataset, compression=False):
    """ Concatenates a matrix (dataset) called dataset_name to an openned HDF5 file (storage), creating if necessary """

    # If the dataset doesn't exist, create it!
    if dataset_name not in storage:
        shape = list(N.shape(dataset))
        shape[0]=None
        if compression:
            storage.create_dataset(dataset_name, data=dataset, maxshape=shape, compression='gzip')
        else:
            storage.create_dataset(dataset_name, data=dataset, maxshape=shape)
    # Otherwise, reshapes the dataset in the HDF5 storage and add the new data
    else:
        old_length = storage[dataset_name].shape[0]
        storage[dataset_name].resize(old_length+len(dataset), axis=0)

        storage[dataset_name][old_length:] = dataset

    return
