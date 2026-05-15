import numpy as np

def flux2mag(flux):
    zeropoint = 31.4 # AB zero-point"
    return -2.5*np.log10(flux) + zeropoint
