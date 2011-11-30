#!/usr/bin/env python

import subprocess as sp

Ds=range(0,6)
Ls=range(10,210,10)
Ss=[10000]
#Ts=['','-T 0','-T 1']
Ts=['-T 0','-T 1']
Rf="./lsa_sim.py"
Cf="%s %s -D %d -L %d -S %d %s"

for D in Ds:
  for L in Ls:
    for S in Ss:
      for T in Ts:
        if T == '':
          Of="../test/lsa_sim_D%d_L%d_S%d.txt"
          Of_tmp = Of % (D, L, S)
        else:
          Of="../test/lsa_sim_D%d_L%d_S%d_T%s.txt"
          Of_tmp = Of % (D, L, S, T[-1])
        Cf_tmp = Cf % (Rf, Of_tmp, D, L, S, T)
        print Cf_tmp, "..."
        sp.Popen( Cf_tmp, shell=True).communicate()
