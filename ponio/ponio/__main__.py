#! /usr/bin/env python3
# -*- coding: utf-8 -*-

#import ponio
import os,sys
import argparse

import jinja2 as jinja
import pickle

import numpy as np
import sympy as sp
from sympy import I

import ponio.runge_kutta as rk
import ponio.utils as utils

def runge_kutta_list () :
  """
    Function return list of Runge Kutta methods
    All Butcher tableaus are write as :
    c | A
    ------
      | b
  """
  rk_list = []

  """
    # Runge-Kutta methods order 1
  """

  A = [[0]]
  b = [1]
  c = [0]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="Euler method") )
  del A,b,c

  alpha = sp.Rational(3,4)
  A = [[0,0],
       [alpha,0]]
  b = [0,1]
  c = [0,alpha]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK NSSP(2,1)") )
  del alpha
  del A,b,c

  A = [[0,0],
       [1,0]]
  b = [sp.Rational(1,2),sp.Rational(1,2)]
  c = [0,sp.Rational(1,2)]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK (2,1)") )
  del A,b,c

  """
    # Runge-Kutta methods order 2
  """

  alpha = sp.Rational(1,2)
  A = [[0,0],
       [alpha,0]]
  b = [1-sp.Rational(1,2*alpha),sp.Rational(1,2*alpha)]
  c = [0,alpha]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK (2,2) mid-point rule") )
  del alpha
  del A,b,c

  alpha = 1
  A = [[0,0],
       [alpha,0]]
  b = [1-sp.Rational(1,2*alpha),sp.Rational(1,2*alpha)]
  c = [0,alpha]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK (2,2) trapezoidal rule") )
  del alpha
  del A,b,c

  A = [[ 0                , 0                , 0 ],
       [ sp.Rational(1,2) , 0                , 0 ],
       [ sp.Rational(1,2) , sp.Rational(1,2) , 0 ]]
  b =  [ sp.Rational(1,3) , sp.Rational(1,3) , sp.Rational(1,3) ]
  c =  [ 0 , sp.Rational(1,2) , 1 ]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK SSP(3,2)") )
  del A,b,c
  
  A = [[ 0                , 0 , 0 ],
       [ sp.Rational(1,3) , 0 , 0 ],
       [ 0                , 1 , 0 ]]
  b =  [ sp.Rational(1,2) , 0 , sp.Rational(1,2) ]
  c =  [ 0 , sp.Rational(1,3) , 1 ]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK NSSP(3,2)") )
  del A,b,c
  
  A = [[0,0,0],
       [sp.Rational(1,2),0,0],
       [0,sp.Rational(1,2),0]]
  b = [0,0,1]
  c = [0,sp.Rational(1,2),sp.Rational(1,2)]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK (3,2) best") )
  del A,b,c

  """
    # Runge-Kutta methods order 3
  """

  A = [[0,0,0],
       [1,0,0],
       [sp.Rational(1,4),sp.Rational(1,4),0]]
  b = [sp.Rational(1,6),sp.Rational(4,6),sp.Rational(1,6)]
  c = [0,sp.Rational(1,2),1]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK SSP (3,3)") )
  del A,b,c

  A = [[0,0,0],
       [sp.Rational(1,2),0,0],
       [-1,2,0]]
  b = [sp.Rational(1,6),sp.Rational(4,6),sp.Rational(1,6)]
  c = [0,sp.Rational(1,2),1]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK (3,3)") )
  del A,b,c

  A = [[0,0,0],
       [-sp.Rational(4,9),0,0],
       [sp.Rational(7,6),-sp.Rational(1,2),0]]
  b = [sp.Rational(1,4),0,sp.Rational(3,4)]
  c = [0,-sp.Rational(4,9),sp.Rational(2,3)]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK NSSP (3,3)") )
  del A,b,c

  A = [[ 0 , 0 , 0 , 0 , 0 ]                ,
       [ sp.Rational(1,7)  , 0 , 0 , 0 , 0 ],
       [ 0 , sp.Rational(3,16) , 0 , 0 , 0 ],
       [ 0 , 0 , sp.Rational(1,3)  , 0 , 0 ],
       [ 0 , 0 , 0 , sp.Rational(2,3)  , 0 ]]
  b = [sp.Rational(1,4),0,0,0,sp.Rational(3,4)]
  c = [0,sp.Rational(1,7),sp.Rational(3,16),sp.Rational(1,3),sp.Rational(2,3)]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK NSSP (5,3)") )
  del A,b,c

  A = [[0,0,0],
       [sp.Rational(2,3),0,0],
       [sp.Rational(1,3),sp.Rational(1,3),0]]
  b = [sp.Rational(1,4),0,sp.Rational(3,4)]
  c = [0,sp.Rational(2,3),sp.Rational(2,3)]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK (3,3) (233e)") )
  del A,b,c

  """
    # Runge-Kutta methods order 4
  """
  
  A = [[ 0               , 0               , 0 , 0 ],
       [ sp.Rational(1,2), 0               , 0 , 0 ],
       [ 0               , sp.Rational(1,2), 0 , 0 ],
       [ 0               , 0               , 1 , 0 ]]
  b = [sp.Rational(1,6),sp.Rational(1,3),sp.Rational(1,3),sp.Rational(1,6)]
  c = [0,sp.Rational(1,3),sp.Rational(2,3),1]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK (4,4)") )
  del A,b,c

  A = [[ 0               , 0 , 0 , 0 ],
       [ sp.Rational(1,3), 0 , 0 , 0 ],
       [-sp.Rational(1,3), 1 , 0 , 0 ],
       [ 1               ,-1 , 1 , 0 ]]
  b = [sp.Rational(1,8),sp.Rational(3,8),sp.Rational(3,8),sp.Rational(1,8)]
  c = [0,sp.Rational(1,3),sp.Rational(2,3),1]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK (4,4) 3/8-rule") )
  del A,b,c

  A = [[0,0,0,0],
       [sp.Rational(1,4),0,0,0],
       [0,sp.Rational(1,2),0,0],
       [1,-2,2,0]]
  b = [sp.Rational(1,6),0,sp.Rational(2,3),sp.Rational(1,6)]
  c = [0,sp.Rational(1,4),sp.Rational(1,2),1]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK (4,4) (235j)") )
  del A,b,c

  """
    # Runge-Kutta methods order 5
  """
  
  A = [[ 0                 , 0                , 0                , 0                 , 0                , 0 ],
       [ sp.Rational(1,4)  , 0                , 0                , 0                 , 0                , 0 ],
       [ sp.Rational(1,8)  , sp.Rational(1,8) , 0                , 0                 , 0                , 0 ],
       [ 0                 , 0                , sp.Rational(1,2) , 0                 , 0                , 0 ],
       [ sp.Rational(3,16) ,-sp.Rational(3,8) , sp.Rational(3,8) , sp.Rational(9,16) , 0                , 0 ],
       [-sp.Rational(3,7)  , sp.Rational(8,7) , sp.Rational(6,7) ,-sp.Rational(12,7) , sp.Rational(8,7) , 0 ]]
  b =  [ sp.Rational(7,90) , 0                , sp.Rational(32,90) , sp.Rational(12,90) , sp.Rational(32,90) , sp.Rational(7,90) ]
  c = [0,sp.Rational(1,4),sp.Rational(1,4),sp.Rational(1,2),sp.Rational(3,4),1]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK (6,5) (236a)") )
  del A,b,c

  A = [[0                      , 0                      ,0                      , 0                   , 0                      ,0                 ,0],
       [sp.Rational(1,5)       , 0                      ,0                      , 0                   , 0                      ,0                 ,0],
       [sp.Rational(3,40)      , sp.Rational(9,40)      ,0                      , 0                   , 0                      ,0                 ,0],
       [sp.Rational(44,45)     ,-sp.Rational(56,15)     ,sp.Rational(32,9)      , 0                   , 0                      ,0                 ,0],
       [sp.Rational(19372,6561),-sp.Rational(25360,2187),sp.Rational(64448,6561),-sp.Rational(212,729), 0                      ,0                 ,0],
       [sp.Rational(9017,3168) ,-sp.Rational(355,33)    ,sp.Rational(46732,5247), sp.Rational(49,176) ,-sp.Rational(5103,18656),0                 ,0],
       [sp.Rational(35,384)    , 0                      ,sp.Rational(500,1113)  , sp.Rational(125,192),-sp.Rational(2187,6784) ,sp.Rational(11,84),0]]
  b =  [sp.Rational(35,384)    , 0                      ,sp.Rational(500,1113)  , sp.Rational(125,192),-sp.Rational(2187,6784) ,sp.Rational(11,84),0]
  c = [0,sp.Rational(1,5),sp.Rational(3,10),sp.Rational(4,5),sp.Rational(8,9),1,1]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="DP 5") )
  del A,b,c

  """
    # Runge-Kutta methods order 6
  """
  
  A = [[ 0                           , 0                        , 0                         , 0                          , 0                     ,0                    ,0                 ,0],
       [ sp.Rational(1,9)            , 0                        , 0                         , 0                          , 0                     ,0                    ,0                 ,0],
       [ sp.Rational(1,24)           , sp.Rational(1,8)         , 0                         , 0                          , 0                     ,0                    ,0                 ,0],
       [ sp.Rational(1,6)            ,-sp.Rational(1,2)         , sp.Rational(2,3)          , 0                          , 0                     ,0                    ,0                 ,0],
       [ sp.Rational(935,2536)       ,-sp.Rational(2781,2536)   , sp.Rational(309,317)      , sp.Rational(321,1268)      , 0                     ,0                    ,0                 ,0],
       [-sp.Rational(12710,951)      , sp.Rational(8287,317)    ,-sp.Rational(40,317)       ,-sp.Rational(6335,317)      , 8                     ,0                    ,0                 ,0],
       [ sp.Rational(5840285,3104064),-sp.Rational(7019,2536)   ,-sp.Rational(52213,86224)  , sp.Rational(1278709,517344),-sp.Rational(433,2448) ,sp.Rational(33,1088) ,0                 ,0], 
       [-sp.Rational(5101675,1767592), sp.Rational(112077,25994), sp.Rational(334875,441898),-sp.Rational(973617,883796) ,-sp.Rational(1421,1394),sp.Rational(333,5576),sp.Rational(36,41),0]]
  b =  [sp.Rational(41,840),0,sp.Rational(9,35),sp.Rational(9,280),sp.Rational(34,105),sp.Rational(9,280),sp.Rational(9,35),sp.Rational(41,840)]
  c = [0,sp.Rational(1,9),sp.Rational(1,6),sp.Rational(1,3),sp.Rational(1,2),sp.Rational(2,3),sp.Rational(5,6),1]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK (8,6)") )
  del A,b,c

  A = [[  0.0                        ,  0.0                        ,   0.0                        ,  0.0                        ,  0.0                        ,  0.0                        , 0.0 ],
       [  0.202276644898140634933337 ,  0.0                        ,   0.0                        ,  0.0                        ,  0.0                        ,  0.0                        , 0.0 ],
       [  0.075853741836802738100001 ,  0.227561225510408214300004 ,   0.0                        ,  0.0                        ,  0.0                        ,  0.0                        , 0.0 ],
       [  1.359282217283300317252891 , -5.237885702628806615657060 ,   4.753603485345506298404170 ,  0.0                        ,  0.0                        ,  0.0                        , 0.0 ],
       [ -0.321092002258021684715280 ,  1.651353127922382381290896 ,  -0.905286676763720493279991 ,  0.075025551099359796704375 ,  0.0                        ,  0.0                        , 0.0 ],
       [  0.292321839349363565719798 , -0.748269386089829516522437 ,   0.592470844966485039986419 , -0.039554538849143620302490 ,  0.028031240623124531118711 ,  0.0                        , 0.0 ],
       [-20.662761894904085188637368 , 63.852320946332118743247958 , -74.151750947688834248615863 ,  0.864117644373384395219349 , 14.505481659294823706193336 , 16.592592592592592592592588 , 0.0 ]]
  b =  [  0.014285714285714285714286 ,  0.0                        ,   0.0                        ,  0.270899470899470899470899 ,  0.429629629629629629629630 ,  0.270899470899470899470899 , 0.014285714285714285714286 ]
  c = [ 0.0 , 0.202276644898140634933337 , 0.303414967347210952400006 , 0.875 , 0.5 , 0.125 , 1.0 ]
  rk_list.append( rk.butcher(A=A,b=b,c=c,label="RK6ES") )
  del A,b,c

  # end function
  return rk_list


def output_pickle (filename) :
    """
      Export list of Runge-Kutta methods as a pickle python file
    """
    if filename == None :
      filename = "rk_list.pickle"

    rk_list = runge_kutta_list() 
    with open(filename,'wb') as f :
      pickle.dump(rk_list,f,protocol=pickle.HIGHEST_PROTOCOL)

def output_json (filename) :
    """
      Export data on Runge-Kutta methods into json file
    """
    if filename == None:
      filename = "rk.json"
    
    for i,rk in enumerate(runge_kutta_list()) :
      print(rk.label)
      # 1. calculer les données sur rk
      rk_dic = {
        'id': i ,
        'label': rk.label ,
        'order': rk.order() ,
        'stages': rk.stages() ,
        'butcher_tableau': rk.butcher_tableau() ,
        'stability_function': sp.latex(rk.stability_function()) ,
        'scheme': utils.latex_aligned(rk.scheme(shu_osher=True)) ,
        'scheme_lawson': utils.latex_aligned(rk.lawson_scheme(shu_osher=True)) ,
        'stability_domain': [[0,0]] ,
        'order_star': "data:image/png;base64," + "foo", #rk.order_star(base64=True) ,
        'cfl': {
          'weno5': sp.latex(5) ,
          'weno3': sp.latex(3) ,
          'y_max': sp.latex(0)
        }
      }
      print(rk_dic)

      # 2. utiliser jinja pour convertir ceci en json (convertion individuelle)
    # 3. tout écrire dans un fichier json
    
    print("json "+filename)

def output_cpp (filename) :
    """
      Generate C++ functions to test every Runge-Kutta methods
    """
    if filename == None:
      filename = "rk.h"

    for rk in runge_kutta_list() :
      pass
      # 1. calculer le schéma
      # 2. ne conserver qu'un label (nom de fonction, donc supprimer les espaces, parenthèses, etc. et remplacer par des `-` ou `_`) tableau des étapes (format c++) et opération finale (format c++)

    print("c++ "+filename)



if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Ponio scripts for code generation.')
  parser.add_argument('--pickles',
          help="Export Runge-Kutta methods list as pickles without ponio classes",
          action='store_const',const=output_pickle,dest='action')
  parser.add_argument('--json',
          help="Export Runge-Kutta methods data as json format for Web application",
          action='store_const',const=output_json,dest='action')
  parser.add_argument('--c++',
          help="Generates C++ code for each Runge-Kutta methods as function",
          action='store_const',const=output_cpp,dest='action')

  parser.add_argument('--output','-o',
          help="Output file name",
          action='store',dest='filename',default=None)

  arg = parser.parse_args(sys.argv[1:])
  arg.action(arg.filename)

