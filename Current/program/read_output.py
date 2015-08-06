""" Contains functions to read the output files."""
import numpy as N

def read_ts_output(filepath, output_quantities=None):
    """ Reads the filepath which contains the full time series output data
        associated with a single galaxy.
        Optional: output_quantities contains a list of the quantities different
        in lines of the output files _in_order_.
        Returns a dictionary, containing the time series of various quantities,
        with quantity names as keys.
    """

    # Use standard settings if no output_quantities are specified
    if output_quantities==None:
        output_quantities = [
                             'ts_t', # This has less elements...
                             'ts_Br','ts_Bp','ts_alpm','ts_Bzmod'
                             'ts_h', 'ts_om', 'ts_G', 'ts_l', 'ts_v', 'ts_etat',
                             'ts_tau', 'ts_alp_k', 'ts_Uz', 'ts_Ur', 'ts_n',
                             'ts_Beq','ts_rmax',
                             'ts_delta_r', 'ts_alp' # This has less elements...
                            ]
    output_dict = dict()
    with open(filepath) as output_file:
        for output_line, quantity in zip(output_file,output_quantities):
            # lfsr: Should we use N.genfromtxt instead of parsing it manually?
            # I guess the difference is small, anyway...
            output_dict[quantity] = N.array(output_line.split()).astype(float)

    return output_dict

def read_output(filepath):
    """ Reads the filepath which contains the full time series output data
        associated with a single galaxy.
        Returns a dictionary, containing the time series of various quantities,
        with quantity names as keys.
        NB this is actually a wrapper of read_ts_output with a particular
        choice of output_quantities.
    """
    output_quantities = [
                         't',
                         'Br','Bp',
                         'Fr','alp_m',
                         'Bzmod', 'h', 'om', 'G', 'l', 'v',
                         'etat', 'tau', 'alp_k', 'Uz', 'Ur',
                         'n','Beq'
                        ]
    return read_ts_output(filepath, output_quantities=output_quantities)


if __name__ == "__main__"  :
    test_dict = read_ts_output('../output/test/ts_00000001.out')
    for k in test_dict:
        print('{0}: {1} , shape: {2}'.format(k, '',test_dict[k].shape))
    print
    test_dict = read_output('../output/test/run_00000001.out')
    for k in test_dict:
        print('{0}: {1} , shape: {2}'.format(k, '',test_dict[k].shape))
    print




