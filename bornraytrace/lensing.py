import healpy as hp
from scipy.integrate import trapz
from scipy.integrate import simps
from astropy import constants as const
import numpy as np
from astropy import units as u


def kappa_prefactor(H0, om0, length_unit='Mpc'):
    """
    Gives prefactor (3 H_0^2 Om0)/2

    :param H0: Hubble parameter with astropy units
    :param om0: Omega matter
    :param length_unit: for H0 (default Mpc)
    :return: prefactor for lensing

    """

    bit_with_units = H0.to(u.s ** -1)/const.c.to(str(length_unit + '/s'))

    return 1.5 * om0 * bit_with_units * bit_with_units


def raytrace_integration(kappa_prefactor, overdensity_array, a_centre, comoving_edges, mask=None):
    """
    This function evaluates the Born weak lensing integral

    :param kappa_prefactor: defined as the output of the function kappa_prefactor
    :param overdensity_array: an 2D array of overdensity healpix maps in radial shells
    :param a_centre: scale factor at comoving centre of shells
    :param comoving_edges: comoving distance to edges of shells
    :param mask: healpix map where 1 is observed and 0 is mask
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
        overdensity_array = (mask * overdensity_array.T).T

    return np.sum(comoving_prefactors * overdensity_array,axis=1).value


def raytrace(H0, om0, overdensity_array, a_centre, comoving_edges, mask=None, Hubble_length_unit = 'Mpc'):
    """
    Evaluate weak lensing convergence map using Born approximation

    :param H0: Hubble parameter with astropy units
    :param om0: Omega matter
    :param overdensity_array: an 2D array of overdensity healpix maps in radial shells
    :param a_centre: scale factor at comoving centre of shells
    :param comoving_edges: comoving distance to edges of shells
    :param mask: healpix map where 1 is observed and 0 is mask
    :param length_unit: for H0 (default Mpc)
    :return: convergence kappa map
    """

    kappa_pref_evaluated = kappa_prefactor(H0, om0, length_unit = Hubble_length_unit)

    kappa_raytraced = raytrace_integration(kappa_pref_evaluated, overdensity_array, a_centre, comoving_edges, mask)

    return kappa_raytraced


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
    if simpsons:
        normalisation = simps(nz, r_array)
    else:
        normalisation = trapz(nz, r_array)

    nr = nz / normalisation

    q = np.empty(r_array.shape)  # q efficiency  eq. 24 in Kilbinger 15
    for i in range(len(r_array)):
        r = r_array[i]
        integrand = np.multiply(np.divide(r_array[i:] - r, r_array[i:]), nr[i:])
        if simpsons:
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


def shear2kappa(shear_map, lmax, upsample=False):
    """
    Performs Kaiser-Squires on the sphere with healpy spherical harmonics

    :param shear_map: healpix format complex shear map
    :param lmax: maximum ell multipole
    :param upsample: option to upsample map for transforms to mitigate numerical errors
    :return: kappa map
    """

    nside = hp.npix2nside(len(shear_map))
    
    if upsample==True:
        nside=nside*2
        shear_map = hp.ud_grade(shear_map, nside)
        if lmax<=nside:
            print('Upsample warning: lmax = ' + str(lmax) + ' and nside = ' + str(nside))
        
    alms = hp.map2alm([shear_map.real, shear_map.real, shear_map.imag], lmax=lmax, pol=True)
    ell, emm = hp.Alm.getlm(lmax=lmax)

    almsE = alms[1]*((ell*(ell+1.))/((ell+2.)*(ell-1)))**0.5
    almsB = alms[2]*((ell*(ell+1.))/((ell+2.)*(ell-1)))**0.5

    almsE[ell==0] = 0.0
    almsB[ell==0] = 0.0

    almsE[ell==1] = 0.0
    almsB[ell==1] = 0.0

    kappa_E = hp.alm2map(almsE, nside=nside, lmax=lmax, pol=False)
    kappa_B = hp.alm2map(almsB, nside=nside, lmax=lmax, pol=False)
    
    if upsample==True:
        nside=int(nside/2)
        kappa_E = hp.ud_grade(kappa_E, nside)
        kappa_B = hp.ud_grade(kappa_B, nside)

    return kappa_E + 1j*kappa_B


def kappa2shear(kappa_map, lmax, upsample=False):
    """
    Performs inverse Kaiser-Squires on the sphere with healpy spherical harmonics

    :param kappa_map: healpix format complex convergence (kappa) map
    :param lmax: maximum multipole
    :param upsample: option to upsample map for transforms to mitigate numerical errors
    :return: complex shear map (gamma1 + 1j * gamma2)
    """

    nside = hp.npix2nside(len(kappa_map))
    
    if upsample==True:
        nside=nside*2
        kappa_map = hp.ud_grade(kappa_map, nside)
        if lmax<=nside:
            print('Upsample warning: lmax = ' + str(lmax) + ' and nside = ' + str(nside))

    almsE = hp.map2alm(kappa_map.real, lmax=lmax, pol=False)
    almsB = hp.map2alm(kappa_map.imag, lmax=lmax, pol=False)
    ell, emm = hp.Alm.getlm(lmax=lmax)

    kalmsE = almsE / (1. * ((ell * (ell + 1.)) / ((ell + 2.) * (ell - 1))) ** 0.5)
    kalmsE[ell == 0] = 0.0

    kalmsB = almsB / (1. * ((ell * (ell + 1.)) / ((ell + 2.) * (ell - 1))) ** 0.5)
    kalmsB[ell == 0] = 0.0

    _, gamma1, gamma2 = hp.alm2map([kalmsE*0., kalmsE, kalmsB], nside=nside, lmax=lmax, pol=True)
    
    if upsample==True:
        nside=int(nside/2)
        gamma1 = hp.ud_grade(gamma1, nside)
        gamma2 = hp.ud_grade(gamma2, nside)
        
    return gamma1 + 1j*gamma2


def recentre_nz(z_sim_edges, z_samp_centre, nz_input):
    """
    Takes input n(z) sampled at z_samp_centre
    and evaluates interpolated n(z) at new z values
    to match a simulation at z_sim_edges

    :param z_sim_edges: new z values for n(z)
    :param z_samp_centre: original z values for n(z)
    :param nz_input: original n(z)
    :return: new n(z)
    """

    nz_input = np.interp(z_sim_edges[1:],z_samp_centre, nz_input)

    return nz_input/np.sum(nz_input*(z_sim_edges[1:]-z_sim_edges[:-1]))


def get_neighbour_array(nside):
    """
    array of indices labelling the 8 neighbouring pixels for each pixel

    :param nside: nside of map
    :return: neighbour indices array
    """

    neighbour_array = np.empty((hp.nside2npix(nside), 8), dtype=int)
    for i in range(hp.nside2npix(nside)):
        neighbour_array[i] = hp.get_all_neighbours(nside,i)

    return neighbour_array


def peak_find(map_input, nside, neighbour_array=None):
    """
    Find peaks (local maxima) for a given input map

    :param map_input: input map
    :param nside: nside of map
    :param neighbour_array: optional array of indices labelling the 8 neighbouring pixels for each pixel
    :return: list of pixel indices for the peaks
    """

    if neighbour_array==None:
        neighbour_array = get_neighbour_array(nside)

    peak_loc = []
    for i in range(hp.nside2npix(nside)):
        if map_input > np.max(map_input[neighbour_array[i]]):
            peak_loc.append(i)

    return peak_loc
