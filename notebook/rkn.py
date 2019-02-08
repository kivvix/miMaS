#! /usr/bin/env python3

"""
script pour parcourir une portion de C et calculer le module du
polynome caractéristique de RKn en ce point

usage :
    ./rkn.py [n]
    n est faculatif, par défaut 3, si on veut étudier RK19 : ./rkn.py 19 > rk19.dat
    le fichier de sortie est formaté une surface avec gnuplot :
    gnuplot -e "spl rkn.dat w pm; pause -1"
"""

import sys
import sympy as sp
import math

z = sp.symbols("z",complex=True)
def poly_carac(N):
  return sum([ sp.Rational(1,(math.factorial(n)))*z**n for n in range(N+1) ])

def poly_rk76():
  un = sp.symbols("u_n")
  dt = sp.symbols("\\Delta\\ t",real=True)
  lamb = sp.symbols("\\lambda",complex=True)
  z = sp.symbols("z",complex=True)

  def L(u):
    return lamb*u

  nu = sp.symbols("\\nu")
  s21 = sp.sqrt(21)
  
  k1 = dt*L(un)
  k2 = dt*L(un+nu*k1)
  k3 = dt*L(un+ ((4*nu-1)*k1+k2)/(8*nu) )
  k4 = dt*L(un+ ((10*nu-2)*k1 + 2*k2 + 8*nu*k3)/(27*nu) )
  k5 = dt*L(un+ (-((77*nu-56)+(17*nu-8)*s21)*k1
                -8*(7+s21)*k2 + 48*(7+s21)*nu*k3
                -3*(21+s21)*nu*k4)/(392*nu) )
  k6 = dt*L(un+ (-5*((287*nu-56)-(59*nu-8)*s21)*k1
                - 40*(7-s21)*k2 + 320*s21*nu*k3 + 3*(21-121*s21)*nu*k4
                + 392*(6-s21)*nu*k5)/(1960*nu) )
  k7 = dt*L(un+ ( 15*((30*nu-8)-(7*nu*s21))*k1 + 120*k2
                - 40*(5+7*s21)*nu*k3 + 63*(2+3*s21)*nu*k4
                - 14*(49-9*s21)*nu*k5 + 70*(7+s21)*nu*k6)/(180*nu) )
  
  un1_rk6 = un + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180
  return un1_rk6.subs(lamb*dt,z).expand().collect(un)/un.collect(z)

def main(N):
  #p = poly_carac(N)
  p = poly_rk76()
  a,b = sp.symbols("a b",real=True)
  p=p.subs(z,a+sp.I*b)

  #x_max = 1.0 ; x_min = -5.0
  #y_max = 5.0 ; y_min = -y_max
  x_max = 2.0 ; x_min = -4.0
  y_max = 4.0 ; y_min = -y_max
  dx = (x_max-x_min)/100.
  dy = (y_max-y_min)/100.

  for i in range(100):
    for j in range(100):
      ax = dx*i+x_min
      by = dy*j+y_min
      print("%f %f %f"%( ax , by , sp.Abs(p.subs({a:ax,b:by}).evalf()).evalf() ))
    print(" ")

if __name__ == '__main__':
  n = 3
  if len(sys.argv) > 1:
    n = int(sys.argv[1])
  main(n)

