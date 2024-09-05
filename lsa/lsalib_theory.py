import numpy as np
import sys
from .lsalib_core import my_decimal, Rmax_min, Rmax_max, kcut_min, pipi, pipi_inv

def readPvalue(P_table, R, N, x_sd=1., M=1., alpha=1., beta=1., x_decimal=my_decimal):
    try:
        xi = int(np.around(R*M/(x_sd*np.sqrt(alpha*beta*N))*(10**x_decimal)))
    except OverflowError as ValueError:
        return np.nan
    if xi in P_table:
        return P_table[xi]
    elif xi > max(P_table.keys()):
        return 0.
    else:
        return np.nan

def theoPvalue(Rmax, Dmax=0, precision=.001, x_decimal=my_decimal):
    Rmax = np.max((Rmax, Rmax_min))
    Rmax = np.min((Rmax, Rmax_max))
    print("computing p_table with Rmax=", Rmax, file=sys.stderr)
    P_table = dict()
    for xi in range(0, Rmax*10**(x_decimal)+1):
        if xi == 0:
            P_table[xi] = 1
            continue
        x = xi/float(10**(x_decimal))
        xx = x**2
        pipi_over_xx = pipi/xx
        alpha = precision
        B = 2*Dmax+1
        Kcut = np.max((kcut_min, int(np.ceil(.5 - np.log((alpha/(2**B-1))**(1/B)*xx*(1-np.exp(-pipi_over_xx))/8/2) / pipi_over_xx))))
        A = 1/xx
        Rcdf = 0
        P_two_tail = 1.
        for k in range(1, Kcut+1):
            C = (2*k-1)**2
            Rcdf = Rcdf + (A+pipi_inv/C)*np.exp(-C*pipi_over_xx/2)
            P_current = 1 - (8**B)*(Rcdf**B)
            if P_current < 0:
                break
            else:
                P_two_tail = P_current
        P_table[xi] = P_two_tail
    return P_table