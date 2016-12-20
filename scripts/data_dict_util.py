""" Useful routines to deal with SAM_data dictionaries """

def filter_dict(ok, data_dict, ignore_list=[]):
    """
    Filters a dictionary full of numpy (or similar) arrays using 'ok'
    filter array.
    
    Returns a new up-to-date dictionary
    """
    new_dict = dict()

    for key in data_dict:
        if key not in ignore_list:
            new_dict[key] = data_dict[key][ok]

    return new_dict


def filter_dictionary_inplace(ok, dictionary, ignore_list=[]):
    """
    Filters a dictionary full of numpy (or similar) arrays using 'ok'
    filter array.
    
    NB the dictionary is altered in the process. 
    """
    for k in dictionary:
        if k not in ignore_list:
            dictionary[k] = dictionary[k][ok]
    return

def dictionary_concatenate_inplace(dict_a, dict_b):
    """
    Concatenates all the keys of two dictionaries.
    """
    import numpy as np
    keys = dict_a.keys()
    for k in keys:
        if k not in dict_b:
            print 'Deleting key: ', k
            del dict_a[k]
        else:
            dict_a[k] = np.concatenate((dict_a[k],dict_b[k]))
    return

