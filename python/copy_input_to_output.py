#!/usr/bin/env python
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
import sys



def copy_input_to_output(input_file, output_file, final_file=None):
    input_h5 = h5py.File(input_file,'r')
    output_h5 = h5py.File(output_file)

    if final_file:
        final_h5 = h5py.File(final_file)
        for x in ('Output', 'Log'):
            if x in final_h5:
              del final_h5[x]
            output_h5.copy(x, final_h5)
    else:
        final_h5 = output_h5

    for x in ('Input', 'Galform Parameters'):
        if x in final_h5:
          del final_h5[x]
        input_h5.copy(x, final_h5)

    output_h5.close()
    input_h5.close()
    if final_file:
        final_h5.close()
    return


if __name__ == "__main__"  :

    import argparse
    parser = argparse.ArgumentParser(description='Copies the input group from '
                                     'one Magnetizer hdf5 file into another.')
    parser.add_argument("INPUT_FILE", help="Magnetizer input HDF5 file.")
    parser.add_argument("OUTPUT_FILE",help="Magnetizer output HDF5 file.")
    parser.add_argument('-o', '--alternate_output', default=None,
                        help='If present, the previously specified input and '
                        'output will both be written into this third file.')
    args = parser.parse_args()
  
    copy_input_to_output(args.INPUT_FILE, 
                         args.OUTPUT_FILE, 
                         final_file=args.alternate_output)

