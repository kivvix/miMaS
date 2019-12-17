#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
  Module to estimate numericlally CFL number of transport equation. The
  most classical example of this library is estimate of CFL number with
  varius RK(s,n) -- WENO5 couples. In fact the CFL number in that case
  is an estimation with the linear case.

  Some export can also be done in JSON for a Web application.

  A library of various RK(s,n) is include (explicit and implicit).
  There is no code generation for implicit methods.
"""

__version__ = "0.1"

from .utils import *
#from .runge_kutta import *
#from .runge_kutta import runge_kutta
