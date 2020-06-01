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


def kappa_prefactor(H0, Om0, length_unit = 'Mpc'):
    bit_with_units = H0.to(u.s ** -1)/const.c.to(str(length_unit + '/s'))
    return 1.5 * Om0 * bit_with_units * bit_with_units

def raytrace_overdensity(cosmo, h, om, overdensity_shells, bin_centres, a_array, max_index, z_max, mask):

    print('z_max should be roughly: ', 1/a_array[max_index - 1] - 1, z_max)
    kappa_map = np.zeros(overdensity_shells.shape[1])
    comoving_max = h * cosmo.comoving_distance(z_max).value

    comoving_prefactors = ((comoving_max - bin_centres[:]) * bin_centres / comoving_max) / 0.7

    kappa_prefactor_val = kappa_prefactor(cosmo.H(0), om, 'Mpc').value
    print(kappa_prefactor_val)

    shell_width = (bin_centres[1] - bin_centres[0]) / h
    print(shell_width)

    counter = 0
    length = len(np.where(mask > 0.5)[0])
    print(length)

    for i in np.where(mask > 0.5)[0]:
        if counter % int(length / 20) == 0:
            print(counter, length)

        integrand = (comoving_prefactors * overdensity_shells[:, i] / a_array)[:max_index]

        integral_val = np.sum(integrand) * shell_width
        kappa_map[i] = kappa_prefactor_val * integral_val
        counter += 1

    return np.where(mask > 0.5, kappa_map, hp.UNSEEN)


def sum_kappa(kappa_prefactor, overdensity_array, comoving_centre, a_centre, comoving_max, mask, comoving_edges=None):
    """
    This function is meant to give the identical result to raytrace_overdensity, but is meant to be cleaner code
    :param kappa_prefactor: defined as the output of the function kappa_prefactor
    :param overdensity_array: an 2D array of overdensity healpix maps in radial shells
    :param comoving_centre: centre of radial shells
    :param a_centre: scale factor at centre of shells
    :param comoving_max: maximum comoving distance
    :param mask: healpix map where 1 is observed and 0 is mask
    :param comoving_edges: edges of comoving bins
    :return: kappa map
    """

    kappa = np.zeros(overdensity_array[:, 0].shape)

    if comoving_edges is not None:
        dr_array = comoving_edges[1:] - comoving_edges[:-1]
    else:
        dr_array = (comoving_centre[1]-comoving_centre[0]) * np.ones(len(comoving_centre))

    print(dr_array.shape, comoving_max)
    comoving_prefactors = dr_array * (comoving_max - comoving_centre) * comoving_centre / (comoving_max * a_centre)
    print(comoving_prefactors.shape)

    for i in np.where(mask > 0.5)[0]:
        if i % 100000 == 0:
            print(i)
        kappa[i] = np.sum(kappa_prefactor * overdensity_array[i, :] * comoving_prefactors[:])

    return kappa


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