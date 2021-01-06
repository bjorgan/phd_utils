#!/usr/bin/env python

from distutils.core import setup

import subprocess
from glob import glob
from pymcml import MCML_EXECUTABLE_NAME, GPUMCML_EXECUTABLE_NAME


#hack: compile GPU-MCML and MCML
subprocess.run(['make'])

#resulting executables
mcml_exec = 'build/' + MCML_EXECUTABLE_NAME
gpumcml_exec = 'build/' + GPUMCML_EXECUTABLE_NAME

setup(name='pymcml',
      version='1.0',
      description='(GPU)MCML wrapper for Python',
      author='Asgeir Bjorgan',
      py_modules=['pymcml'],
      data_files=[('bin', [mcml_exec, gpumcml_exec]),
                  ('share/gpumcml/', ['gpumcml/executable/safeprimes_base32.txt'])]
     )
