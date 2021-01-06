#!/usr/bin/env python

from numpy.distutils.core import setup

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('photon_transport')

    return config

metadata=dict(
      name='photon_transport',
      version='1.0',
      description='Simulate light transport in skin',
      author='Asgeir Bjorgan',
      packages=['photon_transport'],
      configuration=configuration
     )
setup(**metadata)
