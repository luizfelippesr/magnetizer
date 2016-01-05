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


def density(r_50, M_stars, M_gas, vt, alpha=4.0):
    """ Computes the density at a given radius in g cm^{-3}.
        Input: r_50 -> half mass radius (Mpc), M_stars -> stellar mass (Msun),
               M_gas -> gas mass (Msun), vt -> turbulent speed (km/s)
        Output: rho_sol, r_n (kpc)
        WARNING This is valid if vt is a constant...
    """
    r_n = r_50/(2.0*rs_to_r50)
    rho_sol = cq.density_cgs(r_n, r_50, M_stars, M_gas, sigma=vt, alpha=4.0)

    return rho_sol, r_n * Mpc_to_kpc

def height(r_50, M_stars, M_gas, vt, alpha=4.0):
    """ Computes the density at a given radius in cm^{-3}.
        Input: r_50 -> half mass radius (Mpc), M_stars -> stellar mass (Msun),
               M_gas -> gas mass (Msun), vt -> turbulent speed (km/s)
        Output: h_sol (kpc), r_n (kpc)
        WARNING This is valid if vt is a constant...
    """
    r_n = r_50/rs_to_r50 # half-mass to scale-radius conversion
    h_sol = cq.scaleheight_kpc(r_n, r_50, M_stars, M_gas, sigma=vt, alpha=4.0)

    return h_sol, r_n * Mpc_to_kpc


def Brandt_template(r, Omega_0, r_Omega, n):
    """ Generalized Brandt profile (Brandt -> n=2). Template for fitting.
    """
    return Omega_0  / (1.0+(r/r_Omega)**n)**(1./n)

def Brandt_rotation_curve(r_50, v_50):
    """ From the values of r_50 (in the future r_50, 2*r_50, 0.5*r_50, 0.25*r50)
        and v_50 (likewise) fits the parameters of a Brandt rotation curve.
        WARNING This is, at the moment, only a dummy that returns r_Omega = r_50
        and Uphi = V_disk, for testing purposes only!
        Output: U_phi (km/s), r_Omega (kpc), n
    """
    U_phi = v_50
    r_Omega = r_50* Mpc_to_kpc
    n = 2
    return U_phi, r_Omega, n


def outflow_velocity(SFR_Msun_per_Gyr, height, r_50, density_cgs):
    """ Computes the vertical velocity
        WARNING At the moment, this simply scales the previous advection speed
        Output: Uz_sol (km/s), r_Uz (kpc)
    """
    r_Uz = r_50/rs_to_r50
    Uz_sol = cq.V_ad_avg_km_per_s(SFR_Msun_per_Gyr, height, r_50, density_cgs,
                                  attenuation=True)

    return Uz_sol, r_Uz * Mpc_to_kpc


if __name__ == "__main__"  :

    model_dir = 'test_SAM_output'
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
                 #'v_sol_kms   |'
                 #'r_v_kpc     |'
                 #'n_sol_cm3   |'
                 #'r_n_kpc     |'
                 'Uz_sol_kms  |'
                 'r_Uz_kpc    |'
                 #'h_sol_kpc   |'
                 #'r_h_kpc     |'
                 'r_disk_kpc  |'
                 'v_disk_kpc  |'
                 'r_bulge_kpc |'
                 'v_bulge_kpc |'
                 'r_halo_kpc  |'
                 'v_halo_kpc  |'
                 'nfw_cs1     |'
                 'Mgas_disk   |'
                 'Mstars_disk\n')

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
                #    they will need different profiles...

                r_50_max = max(r_50, r_50_max)

                # Computes and stores

                l_sol, r_l = turbulent_scale()

                v_sol, r_v = turbulent_speed()

                vt = v_sol # Currently, a constant turbulent velocity is assumed

                rho_cgs, r_n = density(r_50, M_stars, M_gas, vt, alpha=4.0)

                m = 1.67372e-24 # g (used hydrogen mass... to be updated...)
                n_sol = rho_cgs/m

                h_sol, r_h = height(r_50, M_stars, M_gas, vt, alpha=4.0)
                Uz_sol, r_Uz = outflow_velocity(SFR, h_sol, r_50, rho_cgs)
                Uz_sol = 2
                
                r_disk = r_50 * Mpc_to_kpc
                v_disk = v_50

                r_bulge = data_dict[t]['rbulge'][select][0] * Mpc_to_kpc
                v_bulge = data_dict[t]['vbulge'][select][0]

                r_halo = data_dict[t]['halo_r_virial'][select][0] * Mpc_to_kpc
                v_halo = data_dict[t]['vchalo'][select][0]


                nfw_cs1 = data_dict[t]['strc'][select][0]

                t_dep_input += '{0:.3e}    '.format(t)
                t_dep_input += '{0:.3e}    {1:.3e}    '.format(l_sol,r_l)
                #t_dep_input += '{0:.3e}    {1:.3e}    '.format(v_sol,r_v)
                #t_dep_input += '{0:.3e}    {1:.3e}    '.format(n_sol,r_n)
                t_dep_input += '{0:.3e}    {1:.3e}    '.format(Uz_sol,r_Uz)
                #t_dep_input += '{0:.3e}    {1:.3e}    '.format(h_sol,r_h)
                t_dep_input += '{0:.3e}    {1:.3e}    '.format(r_disk, v_disk)
                t_dep_input += '{0:.3e}    {1:.3e}    '.format(r_bulge, v_bulge)
                t_dep_input += '{0:.3e}    {1:.3e}    '.format(r_halo, v_halo)
                t_dep_input += '{0:.3e}    '.format(nfw_cs1)
                t_dep_input += '{0:.3e}    '.format(M_gas)
                t_dep_input += '{0:.3e}'.format(M_stars)
                t_dep_input += '\n'

            # Writes the file
            fdep.write(t_dep_input)
        with open(t_indep_input_file, 'w+') as findep:
            t_indep_input_file = 'r_max_kpc\n'
            t_indep_input_file += '{0:.3e}\n'.format(
                                             number_of_r50*r_50_max*Mpc_to_kpc)

            findep.write(t_indep_input_file)






