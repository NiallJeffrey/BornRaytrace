import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import os, sys, gc

from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.cosmology import z_at_value

sys.path = ['../'] + sys.path
import born_raytrace as br

index = int(sys.argv[1])
print(index)

cosmo_par = np.genfromtxt('/global/cscratch1/sd/dominikz/DES_Y3_PKDGRAV_SIMS/grid_run_1/cosmo.par').T

# make list of files
filenames = []
om_list = []
s8_list = []
for om, s8 in zip(cosmo_par[0],cosmo_par[6]):
    for i in range(5):
        filenames.append('cosmo_Om=' + str(om) + '_num=' + str(i) + '_s8=' + str(s8))
        om_list.append(om)
        s8_list.append(s8)

# inputs

om = om_list[index]
s8 = s8_list[index]
h = 0.6736
nside = 1024
zmax= 3.
print(om,s8,filenames[index], flush=True)

directory_input = '/global/cscratch1/sd/dominikz/DES_Y3_PKDGRAV_SIMS/grid_run_1/' + str(filenames[index])
directory_output = '/global/cscratch1/sd/ucapnje/DES_Y3_PKDGRAV_kappa/grid_run_1/' + str(filenames[index])
output_name_base = 'kappa_DES-Y3-shell_z'

# generate filenames and z_bin_edges
log = np.genfromtxt(os.path.join(directory_input, 'DES-Y3.log'))
z_bin_edges = log.T[1][::-1]
z_bin_edges = z_bin_edges[np.where(z_bin_edges<zmax)]

sim_filenames = [filename for filename in os.listdir(directory_input) if '.fits' in filename]
sim_filenames = sorted(sim_filenames)
sim_filenames = sim_filenames[:-1]

# cosmo code

cosmo_fiducial = FlatLambdaCDM(H0= h * 100. * u.km / u.s / u.Mpc, Om0=om)

kappa_pref_evaluated = br.kappa_prefactor(cosmo_fiducial.H0, cosmo_fiducial.Om0, length_unit = 'Mpc')

overdensity_array = np.zeros((len(sim_filenames),hp.nside2npix(nside)), dtype=np.float32)
for i in range(overdensity_array.shape[0]):
    map_read = hp.read_map(os.path.join(directory_input, sim_filenames[i]), verbose=False).astype(np.float32)
    print(sim_filenames[i], z_bin_edges[i], flush=True)
    overdensity_array[i] = hp.ud_grade(map_read/np.mean(map_read)-1.,nside)
    

comoving_edges =  cosmo_fiducial.comoving_distance(z_bin_edges)

z_centre = np.empty(len(comoving_edges)-1)
for i in range(len(comoving_edges)-1):
    z_centre[i] = z_at_value(cosmo_fiducial.comoving_distance,
                             0.5*(comoving_edges[i]+comoving_edges[i+1]))
    
for source_edge_index in np.arange(1,len(z_bin_edges)):
    print(z_bin_edges[source_edge_index], flush=True)
    map_temp = br.raytrace_integration(kappa_prefactor=kappa_pref_evaluated,
                                                        overdensity_array=overdensity_array[:source_edge_index].T,
                                                        a_centre=1./(1.+z_centre[:source_edge_index]),
                                                        comoving_edges=comoving_edges[:source_edge_index+1])
    try:
        hp.write_map(os.path.join(directory_output,str(output_name_base)+str(z_bin_edges[source_edge_index])+'.fits'),
                 map_temp, overwrite=True)
    except:
        print(str(os.path.join(directory_output) + ' does not exist - mkdir command'))
        os.mkdir(os.path.join(directory_output))    
        hp.write_map(os.path.join(directory_output,str(output_name_base)+str(z_bin_edges[source_edge_index])+'.fits'),
                 map_temp, overwrite=True)
        
    gc.collect()
