.. BornRaytrace documentation master file, created by
   sphinx-quickstart on Fri Nov 12 11:32:19 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

BornRaytrace documentation
========================================
https://github.com/NiallJeffrey/BornRaytrace


.. toctree::
   :maxdepth: 1
   :caption: Contents:
   
   modules

* :ref:`genindex`
* :ref:`search`


Simulating weak gravitational lensing effects
================================================================

+ Raytrace through overdensity Healpix maps to return a convergence map

+ Include shear-kappa transformation on the full sphere

+ Include intrinsic alignments (NLA model)


* `Simple example notebook <https://github.com/NiallJeffrey/BornRaytrace/blob/main/demo/simple_lens_demo.ipynb>`_
* `Full example notebook <https://github.com/NiallJeffrey/BornRaytrace/blob/main/demo/full_lensing_demo.ipynb>`_

===========================================

.. image:: https://github.com/NiallJeffrey/BornRaytrace/blob/dev/demo/demo_plot.jpg
    :alt: READMEimage


Requirements (python3): numpy, scipy, astropy, healpy

===========================================

Citation:

If you find this code useful, please cite: "Likelihood-free inference with neural compression of DES SV weak lensing map statistics", Jeffrey, Alsing, Lanusse 2020::
article{2020,
   title={Likelihood-free inference with neural compression of DES SV weak lensing map statistics},
   volume={501},
   ISSN={1365-2966},
   url={http://dx.doi.org/10.1093/mnras/staa3594},
   DOI={10.1093/mnras/staa3594},
   number={1},
   journal={Monthly Notices of the Royal Astronomical Society},
   publisher={Oxford University Press (OUP)},
   author={Jeffrey, Niall and Alsing, Justin and Lanusse, François},
   year={2020},
   month={Nov},
   pages={954–969}
}



