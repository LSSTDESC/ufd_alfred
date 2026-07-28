import numpy as np
import pandas as pd

def flux2mag(flux):
    zeropoint = 31.4
    index = flux.index if hasattr(flux, 'index') else None
    flux = np.asarray(flux, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        mag = -2.5 * np.log10(flux) + zeropoint
    mag[~np.isfinite(mag)] = np.nan
    if index is not None:
        return pd.Series(mag, index=index)
    return mag

# def flux2mag(flux):
#    return (flux*u.nJy).to(u.ABmag).value <-- this instead ??

def fluxerr2magerr(flux, flux_err):
    flux = np.ma.filled(np.ma.asarray(flux, dtype=float), fill_value=np.nan)
    flux_err = np.ma.filled(np.ma.asarray(flux_err, dtype=float), fill_value=np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        magerr = (2.5 / np.log(10)) * (flux_err / flux)
    magerr[~np.isfinite(magerr)] = np.nan
    return magerr

tract_dict={'47 Tuc' : [531, 532, 453, 454],
            'ECDFS' : [4848, 4849, 5063, 4636, 4637, 4638, 4847, 4850, 5061, 5062, 5064,
       5065, 5279, 5280, 5281, 5282], #Euclid's "Euclid Deep Field Fornax" encompasses Extended Chandra Deep Field South
            'EDFS' : [2078, 2079, 2080, 2232, 2233, 2234, 2235, 2236, 2237, 2392, 2393,
       2394, 2395, 2396, 2397, 2557, 2558, 2559, 2560, 2561, 2562, 2728,
       2729, 2730, 2731],
            'Fornax' : [4016, 4017, 4218, 4217],
            'FDSG' : [4016, 4017, 4218, 4217], #same as Fornax idk how it'll be referred to
            'Rubin_SV_095-25' : [5525, 5526],
            'Rubin_SV_38_7' : [10463, 10464, 10704],
            'LELF': [10464, 10221, 10222, 10704, 10705, 10463], # Low Ecliptic Latitude Field / Rubin_SV_38_7
            'Seagull' : [7850, 7849, 7610, 7611],
       }

def get_tract(field):
    '''
    Input: field -- str, name of field in LSST data (case sensitive)
    Output: tracts -- list of ints, all the tracts that lie in the definition of that field
                            (depends on data release,
                            this function will always (try to) be the most up-to-date data release)
    '''
    #for DP1, from https://portal.nersc.gov/cfs/lsst/dp1/contributed-notebooks/DP1_Detector_Visits_NB1.html
    try:
        return tract_dict[field]
    except:
        raise Exception('Field not found. Check capitalization or spelling')

def get_field(tract):
    '''
    Input: tract -- int, single tract number
    Output: field -- str, field in which that tract is defined
    '''
    field_dict = {}
    for f in tract_dict.keys():
        for t in tract_dict[f]:
            field_dict[t] = f
    try:
        return field_dict[tract]
    except:
        raise Exception('Tract not found, check scope of data release')