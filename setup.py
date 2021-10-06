#!/usr/bin/env python

from setuptools import setup, find_packages
import sys

import sys, os
print(os.path.dirname(sys.executable), '\n')

setup(name='bornraytrace',
      version='0.1',
      description='Weak gravitational lensing: born raytrace maps, noise and intrinsic alignments',
      author='Niall Jeffrey',
      url='https://github.com/NiallJeffrey/born_raytrace',
      packages=find_packages(),
      install_requires=[
         "numpy", "astropy", "healpy",
      ])
