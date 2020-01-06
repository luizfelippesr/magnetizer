"""
Utilities interfacing between magnetizer and pandas.
"""
import pandas as pd
import numpy as np


def prepare_DataFrame(quantities, magnetizer_run, redshifts, verbose=False,
                      cache_intermediate=False, cache=True):
    """
    Prepares pandas DataFrame objects from the magnetizer run

    Parameters
    ----------
    quantities : list
        List of strings with the names of the requested quantites
    magnetizer_run : MagnetizerRun
        A MagnetizerRun object containing the run.
    redshifts : list
        List of redshifts to be loaded into the DataFrame

    Returns
    -------
    df : pandas.DataFrame
        Dataframe containing all the selected quantities in selected redshifts.
    """
    if verbose:
      print('Preparing a pandas DataFrame object')
    df = pd.DataFrame()

    for quantity in quantities:
        if verbose:
            print('  Loading ',quantity)
        if 'z' not in df:
            zs = []
        data = []
        for z in redshifts:
            if verbose:
                print('   Reading redshift ', z)

            datum = magnetizer_run.get(quantity, z,
                                       cache=cache,
                                       cache_intermediate=cache_intermediate)
            data = np.append(data, datum)

            if 'z' not in df:
                these_zs = np.ones_like(datum)*z
                zs = np.append(zs, these_zs)

        df[quantity] = data
        if 'z' not in df:
            df['z'] = zs
    return df
