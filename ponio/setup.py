#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
from setuptools import setup, find_packages

import ponio
 
setup(
  name = 'ponio' ,
  version = ponio.__version__ ,

  packages = find_packages() ,

  author = "Josselin Massot",
  author_email = "josselin.massot@univ-rennes1.fr" ,

  description = "Python Objects for Numerical IntegratOr",
  long_description = open('README.rst').read() ,

  install_requires = ["matplotlib","numpy","sympy","jinja2"] ,
  include_package_data = True ,

  url = 'http://github.com/kivvix/miMas/ponio' ,

  # https://pypi.python.org/pypi?%3Aaction=list_classifiers.
  classifiers=[
      "Programming Language :: Python",
      "Development Status :: 1 - Planning",
      "License :: OSI Approved",
      "Natural Language :: French",
      "Operating System :: OS Independent",
      "Programming Language :: Python :: 3",
      "Topic :: Scientific/Engineering :: Mathematics",
  ],
  #entry_points = {
  #    'console_scripts': [
  #        'proclame-sm = sm_lib.core:proclamer',
  #    ],
  #},

  license="WTFPL",
)

