#!/usr/bin/env python
from lsa import lsalib
t_set = [0, 0.5, 1, 2]
v_set = []
nv_set = []
bootNum = 100
for t in t_set:
  if t == 0:
    v_set.append(1.25)
    nv_set.append(1.25)
  else:
    P = lsalib.calc_tmatrix(bootNum, t)
    w, vl, vr = lsalib.calc_eigen(P) 
    sigma_square = lsalib.calc_sigma_square(w, vl, vr)
    v_set.append(sigma_square)
    nv_set.append(lsalib.calc_markov_var(P))
print "old method"
for i in range(0,4):
  print "t=", t_set[i], "v=", v_set[i]
print "new method"
for i in range(0,4):
  print "t=", t_set[i], "v=", nv_set[i]
