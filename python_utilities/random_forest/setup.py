#!/usr/bin/env python

from distutils.core import setup

setup(name='random_forest',
      version='1.0',
      description='Random forest skin segmentation',
      author='Asgeir Bjorgan',
      author_email='asgeir.bjorgan@iet.ntnu.no',
      py_modules=['random_forest'],
      scripts=['hyper_mask_random_forest'],
     )
