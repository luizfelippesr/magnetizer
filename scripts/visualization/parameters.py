""" Contains a class for dealing with the parameters in the magnetizer files """
import re

class parameters(object):
    """
    A simple class that carries dictionaries of parameter values as attributes

    The parameters objects receive an openned Magnetizer output hdf5 file as
    as argument.

    Example:
        # Opens the file
        h5file = h5py.File('test.hdf5')
        # Initializes the object
        params = parameters(h5file)
        # Prints a particular parameter values
        print params.ISM_and_disk['P_NX_REF']
        print params.dynamo['ALG_QUENCH']
    """
    def __init__(self,h5file):

        self.name = h5file.attrs['Model name']
        self.date = h5file.attrs['Run date']
        self.file = h5file

        for p in ('ISM_and_disk','dynamo', 'grid','io','outflow','run'):
            setattr(self, p, self.__parse_parameters__(p+'_parameters'))

    def __parse_parameters__(self,attribute):
        """
        Returns a dictionary with the parameters contained in a Magnetizer
        output file.
        """
        txt = self.file.attrs[attribute][0]
        d = dict()
        for match in re.finditer(r'(\w+)=\"?\s*([\w\.\/]+) *\"?,', txt):
            value = match.group(2)
            try:
                value = float(value)
            except:
                pass
            if value == 'T':
                value = True
            elif value == 'F':
                value = False
            d[match.group(1)] = value

        return d
