import numpy as np
import healpy as hp
from astropy.io import fits

import pyccl as ccl

cosmo_taka = ccl.Cosmology(
    Omega_c=0.233, Omega_b=0.046, h=0.70, sigma8=.82, n_s=0.97,
    transfer_function='bbks')

# Define a galaxy number density n(z)

z_bin_edges = np.genfromtxt('takahashi_mock_example/z.txt')
source_edge_index = 20

pz = (z_bin_edges*0.)
pz[source_edge_index] = 1.

dNdz = 1e4 * pz # Number density distribution
b = np.ones(z_bin_edges.shape) # Galaxy bias (constant with scale and z)

ell_kappa_ccl = np.arange(1, 2048)
lens1 = ccl.WeakLensingTracer(cosmo_taka, dndz=(z_bin_edges, dNdz))
cl_kappa_ccl = ccl.angular_cl(cosmo_taka, lens1, lens1, ell_kappa_ccl)

np.save('takahashi_mock_example/cl_example.npy', np.array([ell_kappa_ccl,cl_kappa_ccl]))
