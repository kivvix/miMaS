#!/usr/bin/env python
import subprocess

Tcrhoc = 0.08

class simu_config:
  def __init__ ( self , Nx , Nv , Tc , Tf , alpha , output_dir ):
    self.Nx = Nx
    self.Nv = Nv
    self.Tc = Tc
    self.Tf = Tf
    self.alpha = alpha
    self.output_dir = output_dir

  def write ( self , config_name="config.init" ):
    with open(config_name,'w') as f :
      f.write("\n".join([ "{} {}".format(k,v) for (k,v) in self.__dict__.items() ]))

def compute_Nv(Tc):
  import math
  return math.floor(16.0/(math.sqrt(Tc)/20.0))
def compute_alpha(Tc):
  import math
  #return 1. - Tcrhoc/Tc
  return 0.2

def dic_init(Tc,d):
  return {'Tc':Tc,'Nv':compute_Nv(Tc),'Nx':135,'alpha':compute_alpha(Tc),'Tf':10.0,'output_dir':d}

configs = [
  simu_config(**dic_init(0.08 ,"compare/alpha2/tm0p08"))  ,
  simu_config(**dic_init(0.1  ,"compare/alpha2/tm0p1"))   ,
  simu_config(**dic_init(0.125,"compare/alpha2/tm0p125")) ,
  simu_config(**dic_init(0.15 ,"compare/alpha2/tm0p15"))  ,
  simu_config(**dic_init(0.175,"compare/alpha2/tm0p175")) ,
  simu_config(**dic_init(0.2  ,"compare/alpha2/tm0p2"))
]

"""
configs.append(simu_config(Tc=1e-3,Nv=10120,Nx=135,Tf=10.0,output_dir="compare/tm3"))
configs.append(simu_config(Tc=1e-2,Nv=3200 ,Nx=135,Tf=10.0,output_dir="compare/tm2"))
configs.append(simu_config(Tc=1e-1,Nv=1012 ,Nx=135,Tf=10.0,output_dir="compare/tm1"))
configs.append(simu_config(Tc=1e-4,Nv=32000,Nx=135,Tf=10.0,output_dir="compare/tm4"))
"""

for c in reversed(configs) :
  print(">--")
  print(" ".join([ "{} {}".format(k,v) for (k,v) in c.__dict__.items() ]))
  c.write("config.init")
  subprocess.run("./cmp_tb.out   config.init".split(),shell=True,check=True)
  subprocess.run("./cmp_vhll.out config.init".split(),shell=True,check=True)

