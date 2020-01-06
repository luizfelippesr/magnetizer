#! /usr/bin/env python3
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
"""
Print to screen summary statistics about errors in a particular output
"""
import h5py
import numpy as np

def error_stats(h5file, error_codes=['e','g','h','H','s','i','p'], quiet=False):
    """
    Prints to screen summary statistics about errors in a particular output.

    Input:  h5file -> h5py file object
            error_codes, optional -> list with exit codes considered to be errors.
                                     Default: ['e','g','h','H','s','i','p']
    """

    name = h5file.attrs['Model name']
    status = h5file['Log']['status']
    completed = h5file['Log']['completed'][:]
    runtime = h5file['Log']['runtime'][:]

    # Initializes dictionaries
    err_number = {c: 0 for c in error_codes}
    err_gals = {c: [] for c in error_codes}
    nice = 0

    completed_indices = np.where(completed>0)
    ngal = completed_indices[0].size

    for i in completed_indices[0]:
        # Reads the status
        gal_status = status[i,:]

        # Counts errors
        ok = True
        for code in error_codes:
            if code in gal_status:
                err_number[code]+=1
                err_gals[code].append(str(i))
                ok = False

        # Counts nice galaxies
        if ok:
            nice+=1
    if not quiet:
        print('\nModel name:',name[0])
        print('File: ',h5file.filename)
        print('Total number of galaxies', completed.size)
        print('Number of completed galaxies', ngal)
        print('\nSummary of errors')
        print(' Code\tFrac \t N\tExamples')
    for c in err_number:
        max_idx = min(len(err_gals[c]),9)
        if not quiet:
            print(' {0} \t{1:.1%}\t {2} \t{3}'.format(c, err_number[c]/ngal,
                                  err_number[c],','.join(err_gals[c][:max_idx])))

    if not quiet:
        print(' ok\t{1:.1%}\t {2}'.format(c, nice/ngal, nice))

        print('\nAverage CPU time per galaxy:', runtime[completed_indices].sum()/ngal)
        print()

    return nice, ngal, completed.size


if __name__ == "__main__"  :
    import argparse

    parser = argparse.ArgumentParser(description='Prints to screen summary '
                  'statistics about errors in a particular Magnetizer output')

    parser.add_argument("MAGNETIZER_OUTPUT",help="Name(s) of the Magnetizer"
                        "output file(s) to be examined.", nargs="+")

    parser.add_argument('-k', "--skip_errors", action="store_true",
                        help="If present, files with errors will be skipped.")

    parser.add_argument('-s', "--only_summary", action="store_true",
                        help="Shows only the summary.")

    args = parser.parse_args()

    report = {}

    total_number_of_galaxies = 0.
    total_nice_galaxies = 0.
    for path in args.MAGNETIZER_OUTPUT:
        nice, ngals, total = 0, 0, 0
        try:
            with h5py.File(path,'r') as f:
                nice, ngals, total = error_stats(f, quiet=args.only_summary)
        except:
            if args.skip_errors:
                print('Error\n')
            else:
                raise
        total_nice_galaxies += nice
        total_number_of_galaxies += ngals
        report[path] = (ngals, total)

    if len(report)>4 or args.only_summary:
        print('\n','-'*23,'Summary','-'*23,'\n')
        maxsize = max([len(k) for k in report.keys()])

        for p in sorted(report.keys()):
            ngals, total = report[p]
            name = p + ' '*(maxsize-len(p))
            print(' {0}\t{1}\t{2}\t{3:.2%}'.format(name,int(ngals),total,ngals/float(total)))

    print('\nTotal number of completed galaxies: {0:d}'.format(int(total_number_of_galaxies)))
    print('Galaxies with errors: {0:.2%}'.format(1.0
                                - total_nice_galaxies/total_number_of_galaxies))
