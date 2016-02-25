""" Contains routines which prepare galaxy input files for the magnetic field
    evolution code. """
import computes_quantities as cq
from read_data import read_time_data

rs_to_r50 = 1.678 # half-mass to scale-radius conversion
Mpc_to_kpc = 1000

def turbulent_scale(*args):
    """ Returns a fit to l = l_sol * exp(r/r_l)
        WARNING At the moment, this is a dummy for testing.
        Output: l_sol (kpc), r_l (kpc)
    """
    l_sol = 0.1 # kpc
    r_l = 1000 # Something big to kill the spatial dependence...
    return l_sol, r_l


def turbulent_speed(*args):
    """ Returns a fit to v = v_sol * exp(r/r_v)
        WARNING At the moment, this is a dummy for testing.
        Output: v_sol (km/s), r_v (kpc)
    """
    v_sol = 10 # km\s
    r_v = 1000 # Something big to kill the spatial dependence...
    return v_sol, r_v


if __name__ == "__main__"  :

    model_dir = 'test_SAM_output'
    odir = '../input/test'

    model_dir = '/home/nlfsr/galform_models/FON'
    odir = '../input/test'

    number_of_r50 = 5.0

    data_dict = read_time_data(model_dir,
                                maximum_final_B_over_T=0.5,
                                minimum_final_stellar_mass=1e10,
                                minimum_final_gas_mass=1e7,
                                number_of_galaxies=20,
                                minimum_final_disk_size=5e-4, # Gpc, i.e. 0.5 kpc
                                empirical_disks=False,
                                ivol_dir='ivol0')
    ts = data_dict['tout']
    IDs = data_dict[ts[0]]['ID']

    for i, ID in enumerate(IDs):
        f_gal_id = '{0:8d}'.format(i+1).replace(' ','0')
        t_dep_input_file = '{1}/time_dep_params_{0}.in'.format(f_gal_id, odir)
        t_indep_input_file='{1}/time_indep_params_{0}.in'.format(f_gal_id, odir)
        header =('time        |'
                 'l_sol_kpc   |'
                 'r_l_kpc     |'
                 'r_disk_kpc  |'
                 'v_disk_kpc  |'
                 'r_bulge_kpc |'
                 'v_bulge_kpc |'
                 'r_halo_kpc  |'
                 'v_halo_kpc  |'
                 'nfw_cs1     |'
                 'Mgas_disk   |'
                 'Mstars_disk |'
                 'SFR\n')

        with open(t_dep_input_file, 'w+') as fdep:
            t_dep_input = header
            r_50_max = 0
            for t in sorted(ts):
                if t<7:
                    continue
                select = data_dict[t]['ID'] == ID
                # Skips missing times..
                if not select.any():
                    continue

                # Loads what is necessary
                M_stars = data_dict[t]['mstars_disk'][select][0]
                M_gas = data_dict[t]['mcold'][select][0]
                SFR = data_dict[t]['mstardot'][select][0]
                r_50 = data_dict[t]['rdisk'][select][0]
                v_50 = data_dict[t]['vdisk'][select][0]
                # NB mstars_bulge and mcold_burst not considered
                #    they would need different profiles...

                r_50_max = max(r_50, r_50_max)

                # Computes and stores

                l_sol, r_l = turbulent_scale()
                v_sol, r_v = turbulent_speed()

                r_disk = r_50 * Mpc_to_kpc
                v_disk = v_50

                r_bulge = data_dict[t]['rbulge'][select][0] * Mpc_to_kpc
                v_bulge = data_dict[t]['vbulge'][select][0]

                r_halo = data_dict[t]['halo_r_virial'][select][0] * Mpc_to_kpc
                v_halo = data_dict[t]['vchalo'][select][0]

                nfw_cs1 = data_dict[t]['strc'][select][0]


                t_dep_input += '{0:.3e}    '.format(t)
                t_dep_input += '{0:.3e}    '.format(l_sol)
                t_dep_input += '{0:.3e}    '.format(r_l)
                t_dep_input += '{0:.3e}    '.format(r_disk)
                t_dep_input += '{0:.3e}    '.format(v_disk)
                t_dep_input += '{0:.3e}    '.format(r_bulge)
                t_dep_input += '{0:.3e}    '.format(v_bulge)
                t_dep_input += '{0:.3e}    '.format(r_halo)
                t_dep_input += '{0:.3e}    '.format(v_halo)
                t_dep_input += '{0:.3e}    '.format(nfw_cs1)
                t_dep_input += '{0:.3e}    '.format(M_gas)
                t_dep_input += '{0:.3e}    '.format(M_stars)
                t_dep_input += '{0:.3e}    '.format(SFR)
                t_dep_input += '\n'

            # Writes the file
            fdep.write(t_dep_input)
        with open(t_indep_input_file, 'w+') as findep:
            t_indep_input_file = 'r_max_kpc\n'
            t_indep_input_file += '{0:.3e}\n'.format(
                                             number_of_r50*r_50_max*Mpc_to_kpc)

            findep.write(t_indep_input_file)






