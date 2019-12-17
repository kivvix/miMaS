#! /usr/bin/env python3

import numpy as np
import sympy as sp
from scipy.optimize import newton
from scipy import special

def maxwellian ( rho , u , T ):
  return lambda v : rho/(np.sqrt(2.*np.pi))*np.exp((v-u)**2/(2.*T))

def Z (ksi):
  return np.sqrt(np.pi)*np.exp(-ksi**2)*(1j-erfi(ksi))

def Z_sp (ksi):
  return sp.sqrt(sp.pi)*sp.exp(-ksi**2)*(sp.I - sp.erfi(ksi)) # définition d'Anaïs et Éric
  #return sp.I*sp.sqrt(sp.pi/2)*sp.exp(-ksi**2/2)*(1-sp.erf(-sp.I/sp.sqrt(2)*ksi)) # définition code Maple

def find_zeros ( D , dD ):
  """
    applique juste la méthode de Newton sur tout une zone sur D avec pour dérivée dD
  """
  zeros = set() # un `set` pour n'avoir qu'une seule fois les valeurs de sortie de la méthode de Newton (on ne stocke pas les valeurs initiales)
  zeros_ = set() # puisqu'on travaille sur des flottants avec une méthode itérative, ce `set` stocke les valeurs approchées qui seront toujours les mêmes
  xx = np.linspace(-2,2,11)
  yy = np.linspace(-1,1,11)
  X,Y = np.meshgrid( xx,yy )
  Z = X + 1j*Y
  for zinit in Z.flatten():
    try:
      z0 = newton( D , zinit , tol=1e-8 )
      z0_ = np.complex(float(f"{np.real(z0):.4f}"),float(f"{np.imag(z0):.4f}"))
      # si la valeur approchée n'est pas déjà dans notre `set` des valeurs approchées
      if z0_ not in zeros_ :
        zeros.add(z0)
        zeros_.add(z0_) # on sotcke le nouveau zéro
    except: # si la méthode de newton renvoie une erreur
      pass
  return zeros


# helpers pour écrire la relation de dispersion d'une somme de Dirac et de Maxwelliennes
def diff_dirac(u):
  r"""
    calcul de $\int_R \frac{d_vf}{v-\frac{\omega}{k}}\,\mathrm{d}v$ pour $f(v) = \delta_u(v)
  """
  return lambda k,w:-1/(u-(w/k))**2

def diff_maxwellian(rho,u,T):
  r"""
    calcul de $\int_R \frac{d_vf}{v-\frac{\omega}{k}}\,\mathrm{d}v$ pour $f(v) = M_{\rho,u,T}(v)
  """
  return lambda k,w:-rho/T*( 1 + ((w/k)-u)/(sp.sqrt(2*T))*Z_sp(((w/k)-u)/sp.sqrt(2*T)) )

def numerize(D,param):
  """
    lambdifie l'expression D en une fonction `numpy` avec les substitutions de variables de `param`
  """
  return sp.lambdify( w , D.subs(param) ,  modules=['numpy',{'erf':special.erf,'erfi':special.erfi}] )


def slope (D,param,display=False):
  """
    cherche les zéros de la relation de dispersion `D`, avec les paramètres d'initialisation de `param`
    puis renvoie la plus petite partie imaginaire en valeur absolue
  """
  np.seterr(all='ignore')
  nDkw = numerize(D,param)
  zeros = find_zeros(nDkw,None)
  if display:
    print("\n".join([ f"\033[1m{np.real(z): .4f} {np.imag(z):+.4f}i\033[0m\t\033[37m{np.abs(nDkw(z))}\033[0m" for z in zeros if np.abs(nDkw(z))<0.5]))
  slope = min([ np.abs(np.imag(z)) for z in zeros if np.abs(nDkw(z))<0.5 ]) # `if np.abs(nDkw(z))<1` pour tester si Newton a convergé
  return slope


# calcul avec SymPy de D(k,w), et sa dérivée par rapport à w
a,u = sp.symbols(r"\alpha u",real=True,positive=True)
Tc = sp.symbols(r"T_c",real=True,positive=True)

# relations de dispersions pour différents cas tests
D_d2m = lambda k,w: 1 - 1/(k**2)*( (1-a)*diff_dirac(0)(k,w) + diff_maxwellian(a/2,u,1)(k,w) + diff_maxwellian(a/2,-u,1)(k,w) )
D_3m  = lambda k,w: 1 - 1/(k**2)*( diff_maxwellian(1-a,0,Tc)(k,w) + diff_maxwellian(a/2,u,1)(k,w) + diff_maxwellian(a/2,-u,1)(k,w) )

D_db     = lambda k,w: 1 - 1/(k**2)*( diff_maxwellian(sp.Rational(1,2),-2,1)(k,w) + diff_maxwellian(sp.Rational(1,2),2,1)(k,w) )
D_landau = lambda k,w: 1 - 1/(k**2)*( diff_maxwellian(1,0,1)(k,w) )
D_bot    = lambda k,w: 1 - 1/(k**2)*( diff_maxwellian(1-a,0,Tc)(k,w) + diff_maxwellian(a,u,sp.Rational(1,4))(k,w) )
D_tb     = lambda k,w: 1 - 1/(k**2)*( diff_maxwellian(1-a,0,Tc)(k,w) + diff_maxwellian(a/2,u,1.)(k,w) + diff_maxwellian(a/2,-u,1.)(k,w) )
D_tbvphl = lambda k,w: 1 - 1/(k**2)*( (1-a)*(k/w)**2                 + diff_maxwellian(a/2,u,1.)(k,w) + diff_maxwellian(a/2,-u,1.)(k,w) )
D_2m     = lambda k,w: 1 - 1/(k**2)*( diff_maxwellian(a/2,u,1.)(k,w) + diff_maxwellian(a/2,-u,1.)(k,w) )

D = D_tb

k = sp.symbols(r"k",real=True,positive=True)
w = sp.symbols(r"\omega")
Dkw = D(k,w)
#dwDkw = Dkw.diff(w) # on peut calculer la dérivée directement avec `SymPy` mais 



param = [(k,0.5),(a,0.2),(u,2.),(Tc,0.01)]
nDkw = numerize(Dkw,param)
ndwDkw = None #lambdify(dwDkw,k_num,param)

print("VPHL")
print( slope(D_tbvphl(k,w),param,display=True) )
print("VP-FK")
print( slope(Dkw,param,display=True) )


"""
param = [(k,0.5),(a,0.2),(u,4.)]
temperatures = [0.1,0.01,0.0075,0.005,0.0025,0.00125,0.001,0.0001]
Dkw = D_tb(k,w).expand().simplify()
slopes_kin  = [ slope(Dkw,param+[(Tc,Tc_num)])     for Tc_num in temperatures ]
Dkw = D_tbvphl(k,w).expand().simplify()
slopes_vphl = [ slope(Dkw,param) ]*len(temperatures)

print("\n".join([
  f"\033[37m{Tc_num:.4f}\033[0m : \033[1m{s1:.3f}\033[0m\t\033[1m{s2:.3f}\033[0m"
  for Tc_num,s1,s2 in zip(temperatures,slopes_kin,slopes_vphl)
]))
"""
