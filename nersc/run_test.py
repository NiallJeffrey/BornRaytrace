import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import os, sys, gc

from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.cosmology import z_at_value


cosmo_par = np.genfromtxt('/global/cscratch1/sd/dominikz/DES_Y3_PKDGRAV_SIMS/grid_run_1/cosmo.par').T

filenames = []
om_list = []
s8_list = []
for om, s8 in zip(cosmo_par[0],cosmo_par[6]):
    for i in range(5):
        filenames.append('cosmo_Om=' + str(om) + '_num=' + str(i) + '_s8=' + str(s8))
        om_list.append(om)
        s8_list.append(s8)
        
print(len(filenames))

bad_index = []
for index in range(len(filenames)):
    test = os.listdir('/global/cscratch1/sd/ucapnje/DES_Y3_PKDGRAV_kappa/grid_run_1/' + str(filenames[index]))
    test = sorted(test)
    z_max_temp = float(test[-1].split("z",1)[1][:-5])
    
    log = np.genfromtxt(os.path.join('/global/cscratch1/sd/dominikz/DES_Y3_PKDGRAV_SIMS/grid_run_1/' + str(filenames[index]), 'DES-Y3.log'))
    z_bin_edges = log.T[1][::-1]
    z_bin_edges = z_bin_edges[np.where(z_bin_edges<3.)]
    if z_max_temp<z_bin_edges[-1]:
        print(index, filenames[index], z_max_temp, z_bin_edges[-1])
        bad_index.append(index)
        
print(bad_index,sep=',')

if len(bad_index)==0:
    print('All files completed')

for bad_index_val in bad_index:
    print('sbatch --array=' + str(bad_index_val) + ' born_job.sh')
    
