#!/usr/bin/env python3

import numpy as np
import scipy as sp

# Constants
disp_decimal = 8
kcut_min = 100
Rmax_min = 10
Rmax_max = 50
my_decimal = 2
pipi = np.pi**2
pipi_inv = 1/pipi
Q_lam_step = 0.05
Q_lam_max = 0.95

# Global variables
P_table = dict()

try:
    from lsa import compcore
    np.seterr(all='ignore')
except ImportError:
    from . import compcore
    np.seterr(all='warn')

# Add any other necessary imports or global variables here