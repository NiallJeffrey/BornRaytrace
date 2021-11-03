import numpy as np
import scipy as sp


def E_sq(z, om0):
    """
    A function giving Hubble's law for flat cosmology
    :param z: redshift value
    :param om0: matter density
    :return: A value for the Hubble parameter
    """
    return om0 * (1 + z) ** 3 + 1 - om0


def f_integrand(z, om0):
    """
    A function for the redshift integrand in the intrinsic alignment calculation
    :param z: redshift value
    :param om0: matter density
    :return: redshift integrand
    """
    return (z + 1) / (E_sq(z, om0)) ** 1.5


def D_single(z, om0):
    """
    Provides the normalised linear growth factor
    :param z: single redshift value
    :param om0: matter density
    :return: normalised linear growth factor
    """
    first_integral = sp.integrate.quad(f_integrand, z, np.inf, args=(om0))[0]
    second_integral = sp.integrate.quad(f_integrand, 0, np.inf, args=(om0))[0]

    return (E_sq(z, om0) ** 0.5) * first_integral / second_integral


def D_1(z, om0):
    """
    Provides the normalised linear growth factor
    :param z: single redshift value or array values
    :param om0: matter density
    :return: normalised linear growth factor
    """
    
    if (isinstance(z, float)) or (isinstance(z, int)):
        D_values = D_single(z, om0)
    else:
        z = list(z)
        D_values = [D_single(z[i], om0) for i in range(len(z))]
        D_values = np.array(D_values)
    
    return D_values


def F_nla(z, om0, A_ia, rho_c1, eta=0., z0=0., lbar=0., l0=1e-9, beta=0.):
    """
    The remaining function after solving for the weighting used in the intrinsic alignment
    :param z: redshift value
    :param om0: matter density
    :param A_ia: strength of the intrinsic alignment signal. A free parameter.
    :param rho_c1: rho_crit x C1 (C1 approx 1.508e+27 cm3 / g)
    :param eta: redshift dependence
    :param z0: arbitrary redshift pivot parameter
    :param lbar: average luminosity of source galaxy population
    :param l0: arbitrary luminosity pivot parameter
    :param beta: luminosity dependence
    :return: NLA F(z) amplitude
    """
    
    prefactor = - A_ia * rho_c1 * om0 
    inverse_linear_growth = 1. / D_1(z, om0)
    redshift_dependence = ((1+z)/(1+z0))**eta
    luminosity_dependence = (lbar/l0)**beta
    
    return prefactor * inverse_linear_growth * redshift_dependence * luminosity_dependence

