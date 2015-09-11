""" Contains functions to read the output files."""
import numpy as N


def read_parameter_file(filepath):
    """ Reads one output parameter file.
        Returns a dictionary with parameter values.
    """

    param_dict = {}

    # Recycled contents of rcm.py...
    with open(filepath,'r') as paramf:
        words=paramf.readline().split() #reads line1 into string list
        param_dict['t0_Gyr'] =  float(words[ 0])
        param_dict['t0_kpcskm'] =  float(words[ 1])
        param_dict['h0_kpc'] =  float(words[ 2])
        param_dict['etat0_cm2s'] =  float(words[ 3])
        param_dict['n0_cm3'] =  float(words[ 4])
        param_dict['B0_mkG'] =  float(words[ 5])
        words=paramf.readline().split() #reads line2 into string list
        param_dict['nvar'] =  float(words[ 0])
        param_dict['dt'] =  float(words[ 1])
        param_dict['n1'] =  int  (words[ 2])
        param_dict['n2'] =  int  (words[ 3])
        param_dict['dx'] =  float(words[ 4])
        param_dict['nxphys'] =  int  (words[ 5])
        param_dict['nxghost'] =  int  (words[ 6])
        param_dict['nx'] =  int  (words[ 7])
        words=paramf.readline().split() #reads line3 into string list
        param_dict['l_sol_kpc'] =  float(words[ 0])
        param_dict['r_l_kpc'] =  float(words[ 1])
        param_dict['v_sol_kms'] =  float(words[ 2])
        param_dict['r_v_kpc'] =  float(words[ 3])
        param_dict['n_sol_cm3'] =  float(words[ 4])
        param_dict['r_n_kpc'] =  float(words[ 5])
        param_dict['Uz_sol_kms'] =  float(words[ 6])
        param_dict['r_Uz_kpc'] =  float(words[ 7])
        param_dict['h_sol_kpc'] =  float(words[ 8])
        param_dict['r_h_kpc'] =  float(words[ 9])
        param_dict['Uphi_sol_kms'] =  float(words[10])
        param_dict['r_om_kpc'] =  float(words[11])
        param_dict['R_kappa'] =  float(words[12])
        words=paramf.readline().split() #reads line4 into string list
        param_dict['r_in'] =  float(words[ 0])
        param_dict['r_disk_kpc'] =  float(words[ 1])
        param_dict['r_sol_kpc'] =  float(words[ 2])
        param_dict['r1_kpc'] =  float(words[ 3])
        param_dict['ctau'] =  float(words[ 4])
        param_dict['nn'] =  int  (words[ 5])
        param_dict['lam'] =  float(words[ 6])
        param_dict['C_alp'] =  float(words[ 7])
        param_dict['alpceil'] =  float(words[ 8])
        param_dict['Rm_inv'] =  float(words[ 9])
        words=paramf.readline().split() #reads line5 into string list
        param_dict['etat_sol_cm2s']= float(words[ 0])
        param_dict['td_sol_Gyr'] =  float(words[ 1])
        param_dict['tau_sol_Gyr'] =  float(words[ 2])
        param_dict['etat_sol'] =  float(words[ 3])
        param_dict['td_sol'] =  float(words[ 4])
        param_dict['tau_sol'] =  float(words[ 5])

    return param_dict



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
                             'ts_Br','ts_Bp','ts_alpm','ts_Bzmod',
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

    # If there are multiple times, the arrays need reshaping...
    if 'ts_t' in output_dict:
        n_times = output_dict['ts_t'].size

        for quantity in output_dict:
            if quantity in ('ts_delta_r', 'ts_t'):
                continue
            qsize = output_dict[quantity].size
            output_dict[quantity] = output_dict[quantity].reshape(qsize/n_times,
                                                         n_times)

    return output_dict

def read_final_output(filepath):
    """ Reads the filepath which contains the final output data
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

    param_dict = read_parameter_file('../output/test/param_00000001.out')
    for k in param_dict:
        print('{0}: {1}'.format(k, param_dict[k]))
    print
    test_dict = read_ts_output('../output/test/ts_00000001.out')
    for k in test_dict:
        print('{0}: {1} , shape: {2}'.format(k, '',test_dict[k].shape))
    print
    test_dict = read_final_output('../output/test/run_00000001.out')
    for k in test_dict:
        print('{0}: {1} , shape: {2}'.format(k, '',test_dict[k].shape))
    print
