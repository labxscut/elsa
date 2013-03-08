#!/usr/bin/env python
from lsa import lsalib
x_range = [200, 220, 240, 260, 280, 300, 320, 340,\
    360, 380, 400, 420, 440, 460, 480, 500]
d_set = [0, 1, 2, 3]
p_theo = [ [], [], [], [] ]
x_theo = [ [], [], [], [] ]
for d in d_set:
  p_theo[d] = lsalib.theoPvalue(Rmax=10, Dmax=d)
  for x in x_range:
    x_theo[d].append(p_theo[d][x]) 
for d in d_set:
  print "d=", d, "\t", ",".join(["%.4f" % float(v) for v in x_theo[d]])
