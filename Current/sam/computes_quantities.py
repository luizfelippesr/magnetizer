""" Routines that compute intermediate physical quantites from Galform outputs.
"""
from scipy.optimize import curve_fit
import numpy as N



def scaleheight_kpc(r, r_50, M_stars, M_gas, sigma=10.0, alpha=4.0):
    """ Computes an approximation for the height at a given radius in kpc.
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


def V_ad_km_per_s(r, r_50, SFR_Msun_per_Gyr, scaleheight_kpc,  density_cgs,
                  scaling_type='realistic', nmax=5): # WARNING
    """ Computes the advection speed. It is still not clear how this should
        scale with r (in the initial paper this was already computed from
        averaged quantities). The main problem is the radial dependence of the
        SFR.
    """

    rho_hot_phase = 1.7e-27 # g/cm^3
    hydrogen_mass = 1.67372e-24 # g

    R_kpc = r_50 * 1000 # r_50 is usually in Mpc

    if SFR_scaling=='gas':
        # Assumes SFR scales with density (but with a truncation)
        # rmax -> the maximum possible radius with star formation
        rmax = nmax*r_50
        norm = 2.0
        r_SFR_Msun_per_Gyr = SFR_Msun_per_Gyr * N.exp(-1.678*r/r_50) * norm

    elif SFR_scaling=='realistic':
        # Takes into account the dependence of molecular fraction on r
        # (using e.g. Blitz-Rosolowsky expression)
        exit('Error in V_ad_km_per_s: Realistic SFR scaling not yet implemented.')
    else:
        exit('Error in V_ad_km_per_s: Unrecognized option SFR_scaling={0}'.format(
             SFR_scaling))

    # See equation (19) of Paper 1.
    attenuation_factor = N.minimum(rho_hot_phase/density_cgs,
                                   N.ones(len(density_cgs)))

    V_ad = 3.29e-7 * (r_SFR_Msun_per_Gyr) * (scaleheight_kpc)**(4.0/3.0) \
             * (R_kpc)**2.0 * (density_cgs/hydrogen_mass)**(-1.0/3.0) \
             * attenuation_factor 


    return V_ad



# -----------------------------------------------------------------------------
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
    h_0 = scaleheight_kpc(0.0, r_50, M_stars, M_gas, sigma=sigma, alpha=4.0)
    return h_0 * conversion_factor


def rho_avg_cgs(r_50, M_stars, M_gas, sigma=10.0, alpha=4.0):
    """ Computes the average density """
    conversion_factor = 0.216 # how the central value relates to the average
    rho_0 = density_cgs(0.0, r_50, M_stars, M_gas, sigma=sigma, alpha=4.0)
    return rho_0 * conversion_factor


def V_ad_avg_km_per_s(SFR_Msun_per_Gyr, avg_h_kpc, r_50, density_cgs,
                      attenuation=True):
    """ Computes the helicity advection speed """

    rho_hot_phase = 1.7e-27 # g/cm^3
    hydrogen_mass = 1.67372e-24 # g

    R_kpc = r_50 * 1000 # r_50 is usually in Mpc
    if attenuation:
        attenuation_factor = N.minimum(rho_hot_phase/density_cgs,
                                       N.ones(density_cgs.size))
    else:
        attenuation_factor = 1.0

    V_ad = 3.29e-7 * (SFR_Msun_per_Gyr) * (avg_h_kpc)**(4.0/3.0) * (R_kpc)**2.0 \
             * (density_cgs/hydrogen_mass)**(-1.0/3.0) * attenuation_factor

    return V_ad




# ---------------------------------------------------------------------------
if __name__ == "__main__"  :
    r_50 = N.arange(5.0)

    t =  fit_Brandt_profile(teste, r_50, 1.)
    print t.shape
    print t



