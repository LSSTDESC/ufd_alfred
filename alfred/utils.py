import numpy as np

def flux2mag(flux):
    zeropoint = 31.4 # AB zero-point"
    flux = np.asarray(flux, dtype=float)  # handles masked arrays and lists
    with np.errstate(divide='ignore', invalid='ignore'):
        mag = -2.5 * np.log10(flux) + zeropoint
    mag[~np.isfinite(mag)] = np.nan  # catches inf and nan from log10(<=0)
    return mag

def fluxerr2magerr(flux, flux_err):
    flux = np.ma.filled(np.ma.asarray(flux, dtype=float), fill_value=np.nan)
    flux_err = np.ma.filled(np.ma.asarray(flux_err, dtype=float), fill_value=np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        magerr = (2.5 / np.log(10)) * (flux_err / flux)
    magerr[~np.isfinite(magerr)] = np.nan
    return magerr

def quality_mask(data, snr):
    mask = (data['detect_isIsolated'] == True)
    mask &= (data['r_psfFlux']/data['r_psfFluxErr'] > snr)
    for band in 'griz':
        mask &= (data[f'{band}_psfFlux_flag'] == 0)
        mask &= (data[f'{band}_sizeExtendedness_flag'] == 0)
    return data[mask]
