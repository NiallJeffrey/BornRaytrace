# BornRaytrace
### Simulating weak gravitational lensing effects: 
+ Raytrace through overdensity Healpix maps to return a convergence map
+ Include shear-kappa transformation on the full sphere
+ Include intrinsic alignments (NLA model)

Requirements (python3):
``` numpy```
``` scipy```
``` astropy```
``` healpy```

If you find this code useful, please cite: "Likelihood-free inference with neural compression of DES SV weak lensing map statistics", [Jeffrey, Alsing, Lanusse 2020]


# Born weak lensing formalism:

The weak lensing convergence κ is given by a weighted projection of the density along the line of sight from the observer to a point with radial comoving distance χ and angular position φ on the sky

<img src="https://render.githubusercontent.com/render/math?math=\kappa({\phi}, \chi ) = \frac{3 H_0^2 \Omega_m}{2} \int_0^\chi  \frac{\chi' (\chi - \chi')}{\chi} \frac{\delta({\phi}, \chi')}{a(\chi')} \ \mathrm{d} \chi' ">

where H_0 is the present value of the Hubble parameter, a is the cosmological scale factor, Ω_m is the matter density parameter, δ is the overdensity, and the speed of light c=1. We have assumed flatness, such that the cosmological global curvature is zero, K=0.

For a radial (redshift) distribution n(χ ) of lensed source galaxies, the convergence is given by

<img src="https://render.githubusercontent.com/render/math?math=\kappa ({\phi}) = \int_0^\infty n(\chi) \kappa({\phi}, \chi ) \mathrm{d} \chi = \frac{3 H_0^2 \Omega_m}{2} \int_0^\infty \mathrm{d} \chi' f(\chi')  \chi' \frac{\delta({\phi}, \chi')}{a(\chi')} ">

where 

<img src="https://render.githubusercontent.com/render/math?math=f(\chi') = \int^{\infty}_{\chi^\prime}  \left( \frac{\chi - \chi^\prime}{\chi}\right)n(\chi) \mathrm{d} \chi  "> 

The convergence for the distribution of source galaxies at angular position φ on the sky is therefore given by

<img src="https://render.githubusercontent.com/render/math?math=\kappa({\phi}) = \frac{3 H_0^2 \Omega_m}{2} \int_0^\infty  \Big[ \int_0^\chi\frac{\chi' (\chi - \chi')}{\chi} \frac{\delta({\phi}, \chi')}{a(\chi')}  \mathrm{d} \chi'  \Big] n(\chi) \mathrm{d} \chi "> 

# Validation:

Included in the repository are example Healpix simulated density maps described by (https://arxiv.org/abs/1706.01472) available: http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/

The BornRaytrace code is validated against the results from GRayTrace (http://th.nao.ac.jp/MEMBER/hamanatk/GRayTrix/index.html) using these simulations.

