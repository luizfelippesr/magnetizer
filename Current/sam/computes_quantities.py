""" Routines that compute intermediate physical quantites from Galform outputs.
"""



def h_kpc(r, r_50, M_stars, M_gas, sigma=10.0, alpha=4.0):
    """ Computes the height at a given radius in kpc.
        Input units: r       -> Mpc
                     r_50    -> Mpc
                     M_stars -> Msun
                     M_gas   -> Msun
                     sigma   -> km/s
                     alpha   -> dimensionless
        """

    return 2.622e12 * N.exp(1.68*r/r_50) * sigma**2 * r_50**2 \
            * (M_stars + M_gas)**(-1) * (alpha/4.0)


def density_cgs(r, r_50, M_stars, M_gas, sigma=10.0, alpha=4.0):
    """ Computes the density at a given radius in g/cm^3. """
    return N.float64(3.269e-50)* M_gas/2.0 * (M_gas +M_stars) \
            * sigma**(-2) * r_50**(-4) * (alpha/4)* N.exp(-3.36*r/r_50)


def V_ad_km_per_s(SFR_Msun_per_Gyr, h_kpc, R_kpc, rho_cgs): # WARNING
    """ Computes the advection speed. It is still not clear how this should
        scale with r (in the initial paper this was already computed from
        averaged quantities). The main problem is the radial dependence of the
        SFR.
    """

    rho_hot_phase = 1.7e-27 # g/cm^3
    hydrogen_mass = 1.67372e-24 # g

    R_kpc = r_50 * 1000 # r_50 is usually in Mpc

    attenuation_factor = N.minimum(rho_hot_phase/rho_cgs, N.ones(len(rho_cgs)))

    V_ad = 3.29e-7 * (SFR_Msun_per_Gyr) * (avg_h_kpc)**(4.0/3.0) * (R_kpc)**2.0 \
             * (rho_cgs/hydrogen_mass)**(-1.0/3.0) * attenuation_factor

    return V_ad

  
#-----------------------------------------------------------------------------
#   Average quantities, used in the first paper of the series
# -----------------------------------------------------------------------------
def h_avg_kpc(r_50, M_stars, M_gas, sigma=10.0, alpha=4.0):
    """ Computes the average height in kpc.
        Input units: r_50    -> Mpc
                     M_stars -> Msun
                     M_gas   -> Msun
                     sigma   -> km/s
                     alpha   -> dimensionless
        """
    conversion_factor = 3.28998 # how the central value relates to the average
    h_0 = h_kpc(0.0, r_50, M_stars, M_gas, sigma=sigma, alpha=4.0)
    return h_0 * conversion_factor


def rho_avg_cgs(r_50, M_stars, M_gas, sigma=10.0, alpha=4.0):
    """ Computes the average density """
    conversion_factor = 0.216 # how the central value relates to the average
    rho_0 = rho_cgs(0.0, r_50, M_stars, M_gas, sigma=sigma, alpha=4.0)
    return rho_0 * conversion_factor


def V_ad_avg_km_per_s(SFR_Msun_per_Gyr, avg_h_kpc, r_50, rho_cgs):
    """ Computes the helicity advection speed """

    rho_hot_phase = 1.7e-27 # g/cm^3
    hydrogen_mass = 1.67372e-24 # g

    R_kpc = r_50 * 1000 # r_50 is usually in Mpc

    attenuation_factor = N.minimum(rho_hot_phase/rho_cgs, N.ones(len(rho_cgs)))

    V_ad = 3.29e-7 * (SFR_Msun_per_Gyr) * (avg_h_kpc)**(4.0/3.0) * (R_kpc)**2.0 \
             * (rho_cgs/hydrogen_mass)**(-1.0/3.0) * attenuation_factor

    return V_ad
