#!/usr/bin/env python

from setuptools import setup, find_packages
import sys

import sys, os
print(os.path.dirname(sys.executable), '\n')

setup(name='bornraytrace',
      version='0.2',
      description='Weak gravitational lensing: born raytrace maps, noise and intrinsic alignments',
      author='Niall Jeffrey',
      url='https://github.com/NiallJeffrey/BornRaytrace',
      packages=find_packages(),
      install_requires=[
            "numpy",
            "astropy",
            "healpy",
            "scipy",
      ])
