#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

__all__ = ['contour','contourf']

def contour ( expr , z , Ir , Ii , **kwargs):
  r"""
    function: `contour`
    brief: draw a contour of a complex expression
    params:
      - `expr`:`sympy.Expr` : complex expression to draw
      - `z`:`sympy.Symbol` : symbol on which the expression depends
      - `Ir`:`typle(3)` : real interval to draw
      - `Im`:`typle(3)` : imaginary interval to draw
      - `**kwargs`:`dict` : `pyplot.contour` optional arguments
  """
  import matplotlib.colors as mcolors # for mcolors.TABLEAU_COLORS
  palette = list(mcolors.TABLEAU_COLORS)

  x,y = sp.symbols("x y",real=True)
  f = sp.lambdify((x, y), expr.subs(z,x+sp.I*y),'numpy')
  
  a = np.linspace(*Ir)
  b = np.linspace(*Ii)
  X, Y = np.meshgrid(a, b)
  try:
    lab = kwargs.pop("label")
  except:
    lab = ""
  try:
    c = kwargs["colors"]
  except:
    c = palette[0]
  if len(lab) > 0 :
    plt.plot([0],[0],label=lab,color=c)
  return plt.contour(X,Y, np.real(f(X,Y)) ,**kwargs)

def contourf ( expr , z , Ir , Ii , **kwargs):
  r"""
    function: `contour`
    brief: draw a contour filled of a complex expression
    params:
      - `expr`:`sympy.Expr` : complex expression to draw
      - `z`:`sympy.Symbol` : symbol on which the expression depends
      - `Ir`:`typle(3)` : real interval to draw
      - `Im`:`typle(3)` : imaginary interval to draw
      - `**kwargs`:`dict` : `pyplot.contourf` optional arguments
  """
  import matplotlib.colors as mcolors # for mcolors.TABLEAU_COLORS
  palette = list(mcolors.TABLEAU_COLORS)

  x,y = sp.symbols("x y",real=True)
  f = sp.lambdify((x, y), expr.subs(z,x+sp.I*y),'numpy')
  
  a = np.linspace(*Ir)
  b = np.linspace(*Ii)
  X, Y = np.meshgrid(a, b)
  try:
    lab = kwargs.pop("label")
  except:
    lab = ""
  try:
    c = kwargs["colors"]
  except:
    c = palette[0]
  if len(lab) > 0 :
    plt.plot([0],[0],'s',markersize=1,label=lab,color=c)
  return plt.contourf(X,Y, np.real(f(X,Y)),**kwargs)

def reim(b):
  """
    function to split into 2 arrays real part and imaginary part of complex array
  """
  return ([z.real for z in b],[z.imag for z in b])

def plot(expr,x,I,*args,**kwargs):
  """
    just plot a SymPy expression `expr` of variable `x` on the interval `I`
  """
  X = np.linspace(I[0],I[1],500)
  F = sp.lambdify(x,expr,'numpy')
  plt.plot(X,F(X),*args,**kwargs)

def latex_aligned(list_eq):
  r"""
    function: `latex_aligned`
    brief: return a LaTeX string in `aligned` environnement of a list
     of equations
    params:
      - `list_eq`:`list(sympy.Eq)` : list of `sympy` equations
  """
  return r"\begin{aligned}" + r"\\".join([ sp.latex(eq.lhs) + r"&=" + sp.latex(eq.rhs) for eq in list_eq ]) + r"\end{aligned}"

