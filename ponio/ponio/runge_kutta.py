#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import ponio.utils as poniou

__all__ = ['runge_kutta','butcher','scheme']

def _plot_expression_data (plot_func,Ir,Ii,base64,figsize,**kwargs):
  r"""
    function: `_plot_expression_data`
    brief: helper to manage plot expression if user would a `pyplot`
     object or `png/base64` string
    params:
      - `plot_func`:`lamnd (Ir,Ii,**kwargs)->plt.plot` : plot generator
      - `Ir`:`tuple(3)` : real interval to plot
      - `Ii`:`tuple(3)` : imaginary interval to plot
      - `base64`:`bool` : boolean to choose `png/base64` export
      - `figsize`:`list(2)` : set figsize of
        `matplotlib.pyplot.rcParams['figure.fisize']` parameter
      - `kwargs`:`dict` : see `matplotlib.pyplot.plot` for more
        information
  """
  if not base64:
    return plot_func(Ir,Ii,**kwargs)
  else:
    import io
    import base64
    imgdata = io.BytesIO()
    plt.rcParams['figure.figsize'] = figsize
    plot_func(Ir,Ii,**kwargs)
    plt.axis('equal')
    plt.savefig(imgdata,format='png')
    imgdata.seek(0)
    img = base64.b64encode(imgdata.read())
    plt.close()
    imgdata.close()
    return str(img)[2:-1]


def _name_stages ( u_str , n ):
  u_ = [ sp.symbols(u_str+"^n") ]
  u_.extend([ sp.symbols(u_str+"^{(%i)}"%j) for j in range(1,n) ])
  u_.append(sp.symbols(u_str+"^{n+1}"))
  return u_

class runge_kutta (object) :
  def __init__ (self,n,label=None) :
    r"""
      function: `__init__`
      brief: constructor of a theorical Runge-Kutta method of order $n$
       with $n$ stages.
      params:
        - `n`:`int` : order of a theorical Runge-Kutta method of order
         $n$.
        - `label=None`:`str` : label of method, by default it will be
         `RKn` with `n` the order.
    """
    self.label = label
    self._R = None
    self._data = None
    if self.label == None:
      self.label = "RK%i"%n

  def stability_function (self):
    r"""
      function: `stability_function`
      brief: return the stability function of the method, for an
       explicit Runge-Kutta method of order $n$ with $n$ stages it will
       be the truncate exponential Taylor series:
       $$
         \sum_{k=0}^n \frac{z^k}{k!}.
       $$

       This function could be slow, so result is store to be reuse.
    """
    import math
    if self._R == None:
      z = sp.symbols("z",complex=True)
      self._R = sum([ sp.Rational(1,math.factorial(n))*z**k for k in range(self.order+1) ])
    return self._R

  def data (self,Ir=(-4,2,100),Ii=(-3,3,100)):
    r"""
      function: `data`
      brief: return list of points of the stability domain define as:
       $$
         \mathcal{D} = \left\{ z\in\mathbb{C} , |p(z)|\leq 1 \right\}
       $$
       where $p$ is the stability function of the method.
       
       This function could be slow, so result is store to be reuse.
      params:
       - `Ir=(-4,2,100)`:`tuple(3)` : real interval to explor
        `(xmin,xmax,npts)`
       - `Ii=(-3,3,100)`:`tuple(3)` : imaginary interval to explor
        `(ymin,ymax,npts)`
    """
    if self._data == None:
      fig = plt.figure()
      self._data = poniou.contour(sp.Abs(self.stability_function()),sp.symbols("z",complex=True),Ir,Ii,levels=[1.]).allsegs[0]
      plt.close(fig)
    return self._data

  def stability_domain (self,Ir=(-4,2,100),Ii=(-3,3,100),base64=False,figsize=[6,6],**kwargs):
    r"""
      function: `stability_domain`
      brief: return a plot of stability domain:
       $$
         \mathcal{D}:=\left\{ z\in\mathbb{C}, |p(z)| \leq 1 \right\}
       $$
       By default it returns a plot with `matplotlib.pyplot`, but it
       could be a `base64/png` string.
      params:
       - `Ir=(-4,2,100)`:`tuple(3)` : real interval to explor
        `(xmin,xmax,npts)`
       - `Ii=(-3,3,100)`:`tuple(3)` : imaginary interval to explor
        `(ymin,ymax,npts)`
       - `base64=False`:`bool` : if `True`, returns a `str` of
        `base64\png` representation, else function returns a plot with
        `matplotlib.pyplot`
       - `figsize=[6,6]`:`list(2)` : set figsize of
        `matplotlib.pyplot.rcParams['figure.fisize']` parameter
       - `kwargs`:`dict` : see `matplotlib.pyplot.plot` for more
        information
    """
    def plot(Ir,Ii,**kwargs):
      r"""
        Yet Another Dummy Function
        It is the short version of what we want to plot, argument manager is the goal of `_plot_expression_data`
      """
      return poniou.contour(sp.Abs(self.stability_function()),sp.symbols("z",complex=True),Ir,Ii,levels=[1.],**kwargs)
    return _plot_expression_data(plot_func=plot,Ir=Ir,Ii=Ii,base64=base64,figsize=figsize,**kwargs)


  def order_star (self,Ir=(-5,5,100),Ii=(-5,5,100),base64=False,figsize=[6,6],**kwargs):
    r"""
      function: `order_star`
      brief: return a plot of order star:
       $$
         \mathcal{A}_+:=\left\{ z\in\mathbb{C}, \left|p(z)e^{-z}\right| \geq 1 \right\}
       $$
       By default it returns a plot with `matplotlib.pyplot`, but it
       could be a `base64/png` string.
      params:
       - `Ir=(-5,5,100)`:`tuple(3)` : real interval to explor
        `(xmin,xmax,npts)`
       - `Ii=(-5,5,100)`:`tuple(3)` : imaginary interval to explor
        `(ymin,ymax,npts)`
       - `base64=False`:`bool` : if `True`, returns a `str` of
        `base64\png` representation, else function returns a plot with
        `matplotlib.pyplot`
       - `figsize=[6,6]`:`list(2)` : set figsize of
        `matplotlib.pyplot.rcParams['figure.fisize']` parameter
       - `kwargs`:`dict` : see `matplotlib.pyplot.plot` for more
        information
    """
    def plot(Ir,Ii,**kwargs):
      """
        Yet Another Dummy Function
        It is the short version of what we want to plot, argument manager is the goal of `_plot_expression_data`
      """
      z = sp.symbols("z",complex=True)
      return poniou.contourf(sp.Abs(self.stability_function()*sp.exp(-z)),z,Ir,Ii,levels=[1.,1.e299],**kwargs)
    return _plot_expression_data(plot_func=plot,Ir=Ir,Ii=Ii,base64=base64,figsize=figsize,**kwargs)

  def cfl (self,fourier_symbol,Ir=(-6,4,500),Ii=(-5,5,500),display=True):
    main_domain = max( self.data(Ir,Ii) , key=lambda l:len(l) ) # main domain is the one with the most of points
    main_domain_C = sorted([ z[0]+1j*z[1] for z in main_domain ],key=np.angle) # convert in complex number list of points
    reordered_main_domain_C = np.array([ min(main_domain_C,key=lambda z:np.abs(np.angle(z)-np.angle(-z_fs))) for z_fs in fourier_symbol ]) # for each `z_fs` we find the closest value in term of argument than a boundary element of stability domain

    del main_domain, main_domain_C
    # now we finish to reordonned our data, we can compute the stretching ratio for each value of `fourier_symbol
    dtype = [('angle',float),('sigma_k',float)]
    sigmas = np.array([ (np.angle(z_b),np.abs(z_b/z_fs)) for (z_b,z_fs) in zip(reordered_main_domain_C,fourier_symbol) if np.abs(z_fs)>0.5 and z_b.real<0. ] , dtype=dtype )
    if display:
      sigmas = np.sort(sigmas, order='angle')
      plt.xlabel(r"\varphi_k");plt.ylabel(r"\sigma(\varphi_k)")
      plt.plot(sigmas[:,0],sigmas[:,1])
      plt.show()
    return min(sigmas[:,1])

  def y_max (self):
    y = sp.symbols("y",real=True)

    yset = list(sp.solve(sp.Eq(sp.Abs(self.stability_function().subs(sp.symbols("z",complex=True),sp.I*y))**2,1),y))
    return sp.Abs(max([ sp.re(ysol) for ysol in yset if sp.Abs(sp.im(y))<10**-7 ],key=sp.Abs))



class butcher (runge_kutta) :
  r"""
    `butcher` class represent a Runge-Kutta method from its Butcher
    tableau, with the flowing naming convention:
    $$
      \begin{array}{c|c}
        c & A \\
        \hline
          & b
      \end{array}
    $$
    `A`, `b` and `c` are convert into `sympy.Matrix`. This is the main
    class in `ponio` to represent a Runge-Kutta method beacause we can
    write the underlying Lawson method from it, which is impossible
    from a general prupose of scheme representation.
  """

  def __init__ (self,A,b,c=None,label=None):
    r"""
      function: `__init__`
      brief: constructor method of `butcher` class
      params:
        - `A`:`list(list)|sympy.Matrix` : $A$ matrix
        - `b`:`list|sympy.Matrix` : $b$ row vector
        - `c=None`:`list|sympy.Matrix` : $c$ column vector of time steps
        - `label=None` : label of this method, by default it's build as
          `RK(s,n)`, with `s` the number of stages and `n` the order
    """
    self.A = sp.Matrix(A)
    self.b = sp.Matrix(b)

    if c == None:
      c = [ sum(self.A.row(i)) for i in range(A.cols) ]
    self.c = sp.Matrix(c)

    self.label = label
    if label == None:
      self.label = "RK({},{})".format(self.A.shape[0],butcher.order(self))

    self._R = None           # stability function
    self._scheme = []        # list of stages of the scheme
    self._lawson_scheme = [] # list of stages of the Lawson associated scheme
    self._data = None        # points of border of the stability domain

  def stages (self):
    r"""
      function: `stages`
      brief: return number of stages
    """
    return len(self.c)

  def stability_function (self):
    r"""
      function: `stability_function`
      brief: return the stability_function as a `sympy.Expr` object
    """
    if self._R == None :
      z = sp.symbols("z",complex=True)
      r = range(self.A.shape[0])
      u_s = [ 0 for i in r ]
      for i in r:
        u_s[i] = 1 + sum([ self.A[i,j]*z*u_s[j] for j in r ])
      un1 = 1 + sum([ self.b[j]*z*u_s[j] for j in r ])
      self._R = un1.expand().collect(z)
    return self._R

  def order (self):
    r"""
      function: `order`
      brief: return the order of the method
      description: |
       This method works only for an explicit Runge-Kutta method
       because it works on the stability_function as a polynomial
    """
    import math
    return next(
          ( i-1 for i,r in enumerate(
                [ x-y for (x,y) in zip( reversed(sp.Poly(self.stability_function()).as_list()) , [1./math.factorial(k) for k in range(self.stages()+1)] ) ]
              ) if r ) , # first non nul element in list of (stability_function - exponential_series)
          self.stages() # default value
        )

  def numerical_order (self,tol=1e-5):
    r"""
      function: `numerical_order`
      brief: return the numerical order of the method if coefficients
       are not exact
      description: |
       This method works only for an explicit Runge-Kutta method
       because it works on the stability_function as a polynomial
    """
    import math
    return next(
          ( i-1 for i,r in enumerate(
                [ x-y for (x,y) in zip( reversed(sp.Poly(self.stability_function()).as_list()) , [1./math.factorial(k) for k in range(self.stages()+1)] ) ]
              ) if sp.Abs(r)>tol ) , # first element bigger than `tol` in list of (stability_function - exponential_series)
          self.stages() # default value
        )


  def butcher_tableau (self):
    r"""
      function: `butcher_tableau`
      brief: return a LaTeX string of the Butcher tableau
    """
    tmp = r" \\ ".join([ sp.latex(bi) + " & " +
                         " & ".join([sp.latex(aij) for aij in self.A.row(i)[:i] ]) for i,bi in enumerate(self.b) ])
    tmp += r"\\ \hline & " + " & ".join([ sp.latex(ck) for ck in self.c ])
    return r"\begin{array}{c|"+"c"*len(self.A.row(0))+ r"}" + tmp + r"\end{array}"

  def scheme (self,shu_osher=False) :
    r"""
      function: `scheme`
      brief: return a list of `sympy.Eq` object which represent the
       scheme of the Runge-Kutta method
      params:
        - `shu_osher=False`:`bool` : a boolean to display the Shu-Osher
         scheme (only one call to the $L$ function by stage)
    """
    if self._scheme == []:
      dt = sp.symbols("\\Delta\\ t",real=True)
      L  = sp.Function("L")
      tn = sp.symbols("t^n")
      us = _name_stages("u",len(self.c))
      r = range(len(self.c))
      for i in r:
        u_si = us[0] + dt*sum([ self.A[i,j]*L(tn+self.b[j]*dt,us[j]) for j in r ])
        self._scheme.append(sp.Eq(us[i],u_si))
      un1 = us[0] + dt*sum([ self.c[i]*L(tn+self.b[i]*dt,us[i]) for i in r ])
      self._scheme.append(sp.Eq(us[-1],un1))
      # first line is just u0 = un
      self._scheme = self._scheme[1:]

      if shu_osher :
        L_subs = []
        for i,eq in enumerate(self._scheme[:-1]):
          self._scheme[i] = eq.subs(L_subs).simplify().expand()
          L_subs.append(( L(tn+self.b[i]*dt,us[i]) , sp.solve(self._scheme[i],L(tn+self.b[i]*dt,us[i]))[0] ))
        self._scheme[-1] = self._scheme[-1].subs(L_subs).expand()

    return self._scheme

  def lawson_scheme (self,shu_osher=False) :
    r"""
      function: `lawson_scheme`
      brief: return a list of `sympy.Eq` object which represent the
       scheme of the Lawson method of the underlying Runge-Kutta method
       represented by this object
      params:
        - `shu_osher=False`:`bool` : a boolean to display the Shu-Osher
         scheme (only one call to the $L$ function by stage)
    """
    if self._lawson_scheme == [] :
      dt = sp.symbols("\\Delta\\ t",real=True)
      L = sp.symbols("L",real=True)
      N = sp.Function("N",nargs=1)
      def Nt(t,v):
        """ `Nt` for N tilde, just for `ifrk` scheme """
        return sp.exp(-L*t)*N(sp.exp(L*t)*v)

      tn = sp.symbols("t^n")
      vs = _name_stages("v",len(self.c))
      us = _name_stages("u",len(self.c))
      r = range(len(self.c))
      
      scheme_stages = []
      u_s = [ 0 for i in r ]
      for i in r:
        u_si = vs[0] + dt*sum([ self.A[i,j]*Nt(tn+self.b[j]*dt,vs[j]) for j in r ])
        eq = sp.Eq(vs[i],u_si)
        # on ne prend en considération les étapes différentes de u^0 = u^n.
        if eq != True :
          scheme_stages.append(eq)
      un1 = vs[0] + sum([ dt*self.c[i]*Nt(tn+self.b[i]*dt,vs[i]) for i in r ])
      scheme_stages.append(sp.Eq(sp.symbols("v^{n+1}"),un1))
      
      # substitut all occurences of $v$ by $u$
      vs_us = dict([ (v,u) for (v,u) in zip(vs,us) ])
      vs_usexp = dict([ (v,u*sp.exp(-L*(tn+bs*dt))) for (v,u,bs) in zip(vs,us,self.b) ])
      vs_usexp[vs[-1]] = us[-1]*sp.exp(-L*(tn+dt))
      
      self._lawson_scheme = [ sp.Eq(us,(sp.solve(eq.subs(vs_usexp),us)[0]).simplify().expand()) for (eq,bs,us) in zip(scheme_stages,self.b,us[1:]) ]
      
      # simplification by substitution to evaluate only s times N
      if shu_osher:
        N_subs = []
        for i,eq in enumerate(self._lawson_scheme[:-1]):
          self._lawson_scheme[i] = eq.subs(N_subs).simplify().expand()
          N_subs.append((N(us[i]),sp.solve(self._lawson_scheme[i],N(us[i]))[0]))
        self._lawson_scheme[-1] = self._lawson_scheme[-1].subs(N_subs).expand()

    return self._lawson_scheme

class scheme (runge_kutta) :
  r"""
    `scheme` class represent a Runge-Kutta method from its scheme. This
    is not the main class in `ponio` to represent a Runge-Kutta method
    because it is hard to manipulate a general propose of a scheme
    reprensetation. An example of use of this class:

    ```py
      un,u1,u2 = scheme.ui[0:3]
      L = scheme.L
      rk33_so = scheme(3,label="RK(3,3) Shu-Osher")
      rk33_so[0] = un + dt*L(un)
      rk33_so[1] = sp.Rational(3,4)*un + sp.Rational(1,4)*u1 + sp.Rational(1,4)*dt*L(u1)
      rk33_so[2] = sp.Rational(1,3)*un + sp.Rational(2,3)*u2 + sp.Rational(2,3)*dt*L(u2)

      [sp.pprint(eq) for eq in rk3so.scheme()]
    ```
  """
  ui = _name_stages("u",42)
  L = sp.Function("L")

  def __init__ (self,s,label="") :
    self.s = s
    self._scheme = [0]*s
    self._ui = runge_kutta_scheme.ui[:s]
    self._ui.append(runge_kutta_scheme.ui[-1])
    self.label = label
    self._R = None
  
  def __getitem__ (self,s):
    return self._scheme[s]
  def __setitem__ (self,key,value):
    self._scheme[key] = sp.Eq(self._ui[key],value)

  def scheme (self):
    return self._scheme

  def stability_function (self):
    if self._R == None:
      z = sp.symbols("z",complex=True)
      lamb = sp.symbols("\\lambda",complex=True)
      expr = self._scheme[-1].rhs
      for i,us in enumerate(reversed(self._scheme[:-1])):
        expr = expr.subs(self._ui[self.s-1-i],us.rhs)
      self._R = expr.replace(runge_kutta_scheme.L,lambda x:lamb*x).subs(lamb*sp.symbols("\\Delta\\ t",real=True),z).expand().subs(self._ui[0],1).collect(z)
    return self._R
