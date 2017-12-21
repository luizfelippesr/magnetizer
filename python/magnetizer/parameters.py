# Copyright (C) 2018  Luiz Felippe S. Rodrigues, Luke Chamandy
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

        self.name = h5file.attrs['Model name'][0]
        #self.date = h5file.attrs['Run date'][0]
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
