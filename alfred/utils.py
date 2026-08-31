import numpy as np
import pandas as pd
import yaml
import os

# function to check if the data doesn't exist already and if I want to rewrite it
def check_if_query(path, preload):
    '''
    preload = True means that I want to use the preloaded / saved data instead of querying again
    '''
    if not os.path.exists(path):
        #merged data file doesn't exist yet
        return True
    else:
        #merged data file DOES exist
        if preload == True:
            #I want to use the saved data, so don't remerge them
            return False
        else:
            #I want to overwrite it for whatever reason, so remerge/save them
            return True

def columns_to_query(COLS, bands, output_type='list'):
    INCOLS = []
    for band in bands:
        INCOLS += [col.replace('{band}', band) if '{band}' in col else col for col in COLS]
    INCOLS = list(dict.fromkeys(INCOLS))
    if output_type == 'string' or output_type == 'str':
        INCOLS = ", ".join(INCOLS)
    return INCOLS

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
    #flux = np.ma.filled(np.ma.asarray(flux, dtype=float), fill_value=np.nan)
    index = flux.index if hasattr(flux, 'index') else None
    flux = np.asarray(flux, dtype=float)
    flux_err = np.asarray(flux_err, dtype=float)
    #flux_err = np.ma.filled(np.ma.asarray(flux_err, dtype=float), fill_value=np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        magerr = (2.5 / np.log(10)) * (flux_err / flux)
    magerr[~np.isfinite(magerr)] = np.nan
    if index is not None:
        return pd.Series(mag, index=index)
    return magerr
    
with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    survey = cfg['survey']
    tract_dict = cfg[survey]['field2tract_dict']

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