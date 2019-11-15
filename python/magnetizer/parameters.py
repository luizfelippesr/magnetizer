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
""" Contains a class for dealing with the parameters in the magnetizer files """
import re

class Parameters(object):
    """
    A simple class that carries dictionaries of parameter values as attributes

    The parameters objects receive an open Magnetizer output hdf5 file as
    as argument.

    Usage
    -----
    # Opens the file
    h5file = h5py.File('test.hdf5')
    # Initializes the object
    params = parameters(h5file)
    # Prints a particular parameter values
    print params.ISM_and_disk['P_NX_REF']
    print params.dynamo['ALG_QUENCH']

    Parameters
    ----------
    h5file : h5py.File
        A previously open HDF5 file

    Attributes
    ---------
    name : str
        Name of the model
    file : h5py.File
        The same as h5file parameter
    ISM_and_disk : dict

    dynamo : dict
        Global dynamo parameters
    grid : dict
        Global grid parameters
    io : dict
        Global io parameters
    outflow : dict
        Global outflow parameters
    run : dict
        Global run parameters
    txt : str
        Input data converted into a nicely formated text file
    """
    def __init__(self,h5file):

        self.name = h5file.attrs['Model name'][0]
        #self.date = h5file.attrs['Run date'][0]
        self.file = h5file
        self.txt = ''
        for p in ('ISM_and_disk','dynamo', 'grid','io','outflow','run'):
            setattr(self, p, self.__parse_parameters__(p+'_parameters'))
        #self.txt.replace('/\n\n&', '\n/\n\n&')

    def __parse_parameters__(self,attribute):
        """
        Returns a dictionary with the parameters contained in a Magnetizer
        output file.
        """
        d = {}; out = ''
        # Reads the attributes
        txt = self.file.attrs[attribute][0]
        # Extracts a list of string associated with this name list
        match = re.match('(&\w+)\s+(.+)',txt)
        nml, vals = match.group(1), match.group(2)
        print nml
        vals = vals.split(',')
        vals_str = ['']

        # One needs to account for list of strings in the namelist
        for val in vals:
            # Locates a new parameter
            if '=' in val:
                vals_str[-1] += ','
                # Adjusts the spaces
                val = re.sub('\s+','', val)
                val = val.replace('=',' = ')
                vals_str.append('  '+val)
            else:
                # If there is a list of entries, include them
                vals_str[-1] += ', '+val.replace(' ','')

        # Includes this in the nicely formated text file
        out += nml + '\n'.join(vals_str)

        # Goes through values preparing the dictionary
        for txt in vals_str:
            if '=' not in txt:
                continue
            match = re.search(r'(\w+) = (.+),.*$', txt)
            if match is None:
                continue
            param = match.group(1)
            value = match.group(2)
            try:
                value = float(value)
            except:
                pass
            if value == 'T':
                value = True
            elif value == 'F':
                value = False
                d[param] = value
        # Ajudsts spacing in formated text file
        self.txt += '\n'+out.replace(', /',',\n/')

        return d

    def dump(self, filename=None):
        """
        Saves the parameters to an input file

        Parameters
        ----------
        filname : str
            Name of the file where the input parameters should be saved`
        """
        if filename is None:
            filename = self.name + '.in'
        print('Saving parameters to input file: {}'.format(filename))
        with open(filename, 'w') as f:
            f.write(self.txt)
        return


