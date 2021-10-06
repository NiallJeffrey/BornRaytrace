import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from astropy.io import fits


from scipy.integrate import trapz
from scipy.integrate import simps
import os

from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy import constants as const
from scipy.interpolate import interp1d
from astropy.cosmology import z_at_value

import numpy as np
import scipy as sp
from scipy import integrate
from astropy import units as u

def kappa_prefactor(H0, Om0, length_unit = 'Mpc'):
    bit_with_units = H0.to(u.s ** -1)/const.c.to(str(length_unit + '/s'))
    return 1.5 * Om0 * bit_with_units * bit_with_units


def raytrace_integration(kappa_prefactor, overdensity_array, a_centre, comoving_edges, mask=None):
    """
    This function is meant to give the identical result to raytrace_overdensity, but is meant to be cleaner code
    :param kappa_prefactor: defined as the output of the function kappa_prefactor
    :param overdensity_array: an 2D array of overdensity healpix maps in radial shells
    :param a_centre: scale factor at comoving centre of shells
    :param mask: healpix map where 1 is observed and 0 is mask
    :param comoving_edges: comoving distance to edges of shells
    :return: convergence kappa map
    """
    
    assert overdensity_array.shape[1] + 1 == comoving_edges.shape[0]

    dr_array = comoving_edges[1:] - comoving_edges[:-1]
    comoving_max = comoving_edges[-1]
    comoving_centre = 0.5*(comoving_edges[:-1] + comoving_edges[1:])
    
    comoving_prefactors = dr_array * (comoving_max - comoving_centre) * comoving_centre / (comoving_max * a_centre)
    comoving_prefactors *= kappa_prefactor
    
    if mask is not None:
        mask = np.where(mask>0.5,1.,0.).T
        overdensity_array = mask * overdensity_array
        
    return np.sum(comoving_prefactors * overdensity_array,axis=1).value



def W_kernel(r_array, z_array, nz, simpsons=False):
    """
    lensing kernel W s.t.  kappa = prefactor * integral  W(r) * overdensity(r)  dr
    :param r_array: comoving distances array
    :param z_array: redshift array matching r_array (cosmology dependent)
    :param nz: source redshift distribution
    :param simpsons: boolean to use simpsons integratio
    :return: W = r * q /r
    """

    # normalised redshift distribution nr
    if simpsons == True:
        normalisation = simps(nz, r_array)
    else:
        normalisation = trapz(nz, r_array)

    nr = nz / normalisation

    q = np.empty(r_array.shape)  # q efficiency  eq. 24 in Kilbinger 15
    for i in range(len(r_array)):
        r = r_array[i]
        integrand = np.multiply(np.divide(r_array[i:] - r, r_array[i:]), nr[i:])
        if simpsons == True:
            q[i] = simps(integrand, r_array[i:])
        else:
            q[i] = trapz(integrand, r_array[i:])

    return q * r_array * (1. + z_array)


def rotate_mask_approx(mask, rot_angles, flip=False):
    """
    rotate healpix mask on sphere
    :param mask: healpix map of ones and zeros
    :param rot_angles: rotation on the sphere (e.g. [ 45.91405291 ,150.72092269 , 46.34505909])
    :param flip: boolean, mirror the mask 
    :return: rotated map
    """
    nside = hp.npix2nside(len(mask))
    alpha, delta = hp.pix2ang(nside, np.arange(len(mask)))
    alpha = alpha[np.where(mask>0.)]
    delta = delta[np.where(mask>0.)]
    rot = hp.rotator.Rotator(rot=rot_angles, deg=True)
    rot_alpha, rot_delta = rot(alpha, delta)
    rot_i = hp.ang2pix(nside, rot_alpha, rot_delta*(1-2*float(flip)))
    rot_map = mask*0.
    rot_map[rot_i] += 1 
    return rot_map


def E_sq(z, Om0):
    """
    A function giving Hubble's law for flat cosmology
    :param z: redshift value
    :param Om0: matter density
    :return: A value for the Hubble parameter
    """
    return Om0 * (1 + z) ** 3 + 1 - Om0


def f(z, Om0):
    """
    A function for the redshift integrand in the intrinsic alignment calculation
    :param z: redshift value
    :param H0: Hubble's constant
    :param Om0: matter density
    :return: redshift integrand
    """
    return (z + 1) / (E_sq(z, Om0)) ** 1.5


def D(z, Om0):
    """
    Provides the normalised linear growth factor
    :param z: redshift value
    :param H0: Hubble's constant
    :param Om0: matter density
    :return: normalised linear growth factor
    """
    first_integral = sp.integrate.quad(f, z, np.inf, args=(Om0))[0]
    second_integral = sp.integrate.quad(f, 0, np.inf, args=(Om0))[0]

    return (E_sq(z, Om0) ** 0.5) * first_integral / second_integral


def omega(c1,h,crit_den):
    """

    :param c1: normalisation constant
    :param h: reduced Hubble constant
    :param crit_den: critical density of the universe
    :return: a unit-correct c1* crit_den
    """
    return (crit_den/h**2)*(c1*h**2)


def F(z, Om0, A_ia, eta=None, z0=None, lbar=None, l0=None, beta=None):
    """
    The remaining function after solving for the weighting used in the intrinsic alignment
    :param z: redshift value
    :param Om0: matter density
    :param A_ia: strength of the intrinsic alignment signal. A free parameter.
    :param eta: redshift dependence
    :param z0: arbitrary redshift pivot parameter
    :param lbar: average luminosity of source galaxy population
    :param l0: arbitrary luminosity pivot parameter
    :param beta: luminosity dependence
    :return: a complete function of z, A_ia which can be multiplied by an overdensity to reproduce IA effect
    """

    if eta and z0 is not None:
        return -A_ia*omega(c1)*Om0/D(z, Om0)*((1+z)/(1+z0))**eta
    if lbar and l0 and beta is not None:
        return -A_ia*omega(c1)*Om0/D(z, Om0)*(lbar/l0)**beta
    if eta and z0 and lbar and l0 and beta is not None:
        return -A_ia*omega(c1)*Om0/D(z, Om0)*((lbar/l0)**beta)*(((1+z)/(1+z0))**eta)
    else:
        return -A_ia*omega(c1)*Om0/D(z, Om0)
<<<<<<< HEAD
=======


>>>>>>> 0be1c4892bbf9364448dc7f34153d55f58f6d624
    
    
def shear2kappa(shear_map, lmax=None):
    """
    Performs Kaiser-Squires on the sphere with healpy spherical harmonics
    :param shear_map: healpix format complex shear map
    :param lmax:
    :return: kappa map
    """

    nside = hp.npix2nside(len(shear_map))
    alms = hp.map2alm([shear_map.real, shear_map.real, shear_map.imag], lmax=lmax, pol=True)
    ell, emm = hp.Alm.getlm(lmax=lmax)
    almsE = alms[1] * 1. * ((ell * (ell + 1.)) / ((ell + 2.) * (ell - 1))) ** 0.5

    almsE[ell == 0] = 0.0
    almsE[ell == 1] = 0.0

    kappa_E = hp.alm2map(almsE, nside=nside, lmax=lmax, pol=False)
    kappa_B = hp.alm2map(almsE, nside=nside, lmax=lmax, pol=False)

    return kappa_E + 1j*kappa_B


def shear2kappa(kappa_map, lmax=None):
    """
    Performs inverse Kaiser-Squires on the sphere with healpy spherical harmonics
    :param kappa_map: healpix format complex convergence (kappa) map
    :param lmax:
    :return: complex shear map (gamma1 + 1j * gamma2)
    """

    nside = hp.npix2nside(len(kappa_map))
    # alms = hp.map2alm([kappa_map.real, kappa_map.real, kappa_map.imag], lmax=lmax, pol=True)
    # ell, emm = hp.Alm.getlm(lmax=lmax)

    # almsE = alms[1] / (1. * ((ell * (ell + 1.)) / ((ell + 2.) * (ell - 1))) ** 0.5)
    # almsB = alms[2] / (1. * ((ell * (ell + 1.)) / ((ell + 2.) * (ell - 1))) ** 0.5)
    # almsE[ell == 0] = 0.0
    # almsB[ell == 0] = 0.0

    alms = hp.map2alm(kappa_map.real, lmax=lmax, pol=False)
    ell, emm = hp.Alm.getlm(lmax=lmax)
    kalmsE = alms / (1. * ((ell * (ell + 1.)) / ((ell + 2.) * (ell - 1))) ** 0.5)

    kalmsE[ell == 0] = 0.0

    alms = hp.map2alm(kappa_map.imag, lmax=lmax, pol=False)
    ell, emm = hp.Alm.getlm(lmax=lmax)
    kalmsB = alms / (1. * ((ell * (ell + 1.)) / ((ell + 2.) * (ell - 1))) ** 0.5)

    kalmsB[ell == 0] = 0.0

    _, gamma1, gamma2 = hp.alm2map([kalmsE, kalmsE, kalmsB], nside=nside, lmax=lmax, pol=True)
    return gamma1 + 1j*gamma2
<<<<<<< HEAD

=======
>>>>>>> 0be1c4892bbf9364448dc7f34153d55f58f6d624
