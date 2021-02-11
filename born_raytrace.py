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
        
    return np.sum(comoving_prefactors * overdensity_array,axis=1)



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