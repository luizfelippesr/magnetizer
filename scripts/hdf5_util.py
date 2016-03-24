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
