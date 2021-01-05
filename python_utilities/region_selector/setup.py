#!/usr/bin/env python

from distutils.core import setup

setup(name='region_selector',
      version='1.0',
      description='Hyperspectral rectangular region selector',
      author='Asgeir Bjorgan',
      packages = ['region_selector'],
      scripts = ['region_selector/region_selector.py', 'region_selector/random_forest_from_regions.py']
     )
