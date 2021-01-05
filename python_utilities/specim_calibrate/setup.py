#!/usr/bin/env python

from distutils.core import setup

setup(name='specim_calibrate',
      version='1.0',
      description='Calibration tool for specim hyperspectral files',
      author='Asgeir Bjorgan',
      py_modules=['specim_calibration'],
      scripts=['specim_hyper_calibrate'],
     )
