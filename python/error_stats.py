#! /usr/bin/env python
"""
Print to screen summary statistics about errors in a particular output
"""
from __future__ import division
import h5py
import numpy as np

def error_stats(h5file, error_codes=['e','g','h','H','s','i','p']):
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
    print 'Model name:',name
    print 'File path:',h5file.filename
    print 'Total number of galaxies', completed.size
    print 'Number of completed galaxies', ngal
    print '\nSummary of errors'
    print ' Code\tFrac \t N\tExamples'
    for c in err_number:
        max_idx = min(len(err_gals[c]),9)
        print ' {0} \t{1:.1%}\t {2} \t{3}'.format(c, err_number[c]/ngal,
                                  err_number[c],','.join(err_gals[c][:max_idx]))

    print ' ok\t{1:.1%}\t {2}'.format(c, nice/ngal, nice)

    print '\nAverage CPU time per galaxy:', runtime[completed_indices].sum()/ngal
    print

    return


if __name__ == "__main__"  :
    import argparse

    parser = argparse.ArgumentParser(description='Prints to screen summary '
                  'statistics about errors in a particular Magnetizer output')

    parser.add_argument("MAGNETIZER_OUTPUT",help="Name(s) of the Magnetizer"
                        "output file(s) to be examined.", nargs="+")

    args = parser.parse_args()

    for path in args.MAGNETIZER_OUTPUT:
        with h5py.File(path,'r') as f:
            error_stats(f)
