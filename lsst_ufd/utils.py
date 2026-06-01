import numpy as np

def flux2mag(flux):
    zeropoint = 31.4 # AB zero-point"
    return -2.5*np.log10(flux) + zeropoint

def quality_mask(data, snr):
    mask = (data['detect_isIsolated'] == True)
    mask &= (data['r_psfFlux']/data['r_psfFluxErr']) > snr
    for band in 'griz':
        mask &= (data[f'{band}_psfFlux_flag'] == 0)
        mask &= (data[f'{band}_sizeExtendedness_flag'] == 0)
    return data[mask]
