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

def main(N):
  p = poly_carac(N)
  a,b = sp.symbols("a b",real=True)
  p=p.subs(z,a+sp.I*b)

  x_max = 1.0 ; x_min = -5.0
  y_max = 5.0 ; y_min = -y_max
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

