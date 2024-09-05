import numpy as np
from scipy import signal  # Import signal module from scipy

def FIRuncFilter(y, sigma_y, theta, Utheta=None, shift=0):
    # ... existing code ...
    if Utheta is None:
        Utheta = np.zeros_like(theta)  # Use np.zeros_like instead of np.zeros
    # ... existing code ...
    Uy = signal.lfilter(Utheta, np.array([1.0]), uy)  # Use scipy.signal.lfilter
    # ... existing code ...