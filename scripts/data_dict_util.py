""" Useful routines to deal with SAM_data dictionaries """

def filter_dict(data_dict, ok, ignore_list=[]):
    """ Filters a dictionary full of numpy (or similar) arrays using 'ok' 
    filter array.
    
    Returns a new up-to-date dictionary
    """
    new_dict = dict()
    for key in data_dict:
        if k not in ignore_list:
            new_dict[key] = data_dict[key][ok]
    return new_dict


def filter_dictionary_inplace(ok, dictionary, ignore_list=[]):
    """ Filters a dictionary full of numpy (or similar) arrays using 'ok' 
    filter array.
    
    NB the dictionary is altered in the process. 
    """
    for k in dictionary:
        if k not in ignore_list:
            dictionary[k] = dictionary[k][ok]
    return

