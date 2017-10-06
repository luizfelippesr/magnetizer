#!/usr/bin/env python
import argparse
import numpy as np
import magnetizer
import astropy.units as u
import matplotlib.pyplot as plt

if __name__ == "__main__"  :

    parser = argparse.ArgumentParser(
             description='Prepares portfolio of galaxies in a Magnetizer run, '
             'i.e. samples from stellar mass bins and plots a summary of '
             'evolution of radial profiles for a set of properties.')

    parser.add_argument("MAGNETIZER_OUTPUT",
                        help="HDF5 output of the Magnetizer run (must contain"
                        " an Input group, or a separate input file must be "
                        "specified).")

    parser.add_argument("PDF_OUTPUT",
                        help="Output filename for the PDF containing the "
                        "portfolio.")

    parser.add_argument('-mi',"--magnetizer_input", help="Name of the Magnetizer "
                        "input file (only needed if input and output hdf5 "
                        "are separate)", default=None)

    parser.add_argument('-n', '--number_of_galaxies_per_bin', default=10,
                        help='Number of galaxies to be sampled per mass bin.'
                        'Default: 10')


    parser.add_argument('-b', '--mass_bins', default='8,8.75,9.5,10.25,12',
                        help='Bin edges in logarithm of galaxy stellar mass.'
                        ' Default: 8,8.75,9.5,10.25,12')

    args = parser.parse_args()

    bins_list = np.array(args.mass_bins.split(',')).astype(float)
    mass_bins = [[10.0**m1,10.0**m2]
                  for m1,m2 in zip(bins_list[:-1],bins_list[1:])]

    magnetizer_run = magnetizer.MagnetizerRun(
                                      input_file_path=args.magnetizer_input,
                                      output_file_path=args.MAGNETIZER_OUTPUT)

    binning = magnetizer.MassBinningObject(magnetizer_run, z=0.0,
                                           bins=mass_bins*u.Msun)
    ngals_bin = int(args.number_of_galaxies_per_bin)

    # Some adjustments in the output settings
    plt.rc( ('axes'),labelsize=8)
    plt.rc( ('xtick','ytick'),labelsize=7)
    plt.rc(("legend"),fontsize=6)

    # Prepares the portfolio
    magnetizer.portfolio.generate_portfolio(run_obj=magnetizer_run,
                                            binning_obj=binning,
                                            galaxies_per_bin=ngals_bin,
                                            pdf_filename=args.PDF_OUTPUT)
