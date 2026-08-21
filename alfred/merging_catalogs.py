import healpy as hp
import numpy as np
import gc
import sys
import os
import yaml
from alfred import plotting_functions
from astropy.table import Table, vstack, join
#Goes up a directories to get the updated astroquery
#I really need to fix this, maybe github submodules or enforcing a version of astroquery
#I think it's version 0.4.11 ?
#sys.path.append(os.path.abspath('../'))
from ugali.utils.projector import match
## function to create euclid + rubin datasets and register to the data registry on NERSC
## don't know where I want this to live quite yet

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    #assuming that it's cool that the whole github repo is considered "home"
    where = cfg['setup']['where']
    home_dir = os.path.expandvars(cfg['setup']['home_dir'][where])
    pckg_dir = os.path.join(home_dir, cfg['setup']['pckg_dir'])
    #external data is gonna be in a directory above - subject to change
    data_dir = os.path.join(home_dir, cfg['setup']['data_dir'])
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
    results_dir = os.path.join(home_dir, cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    survey = cfg['survey']
    euclid_survey = cfg['euclid_survey']

# function to check if the data doesn't exist already and if I want to rewrite it
def check_merge_data(tract, preload = True):
    '''
    preload = True means that I want to use the preloaded / saved data instead of querying again
    '''
    if not os.path.exists(data_dir + f'/merged/{tract}_{survey}_{euclid_survey}_merged.parquet'):
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

# function to add it to the data registry

# then function to merge catalogs, starting and ending with above
def merge_catalogs(lsst_table, euclid_table, tract, preload = True, validation_needed = False):
    if not check_merge_data(tract, preload):
        print("Check tells me data exists and you don't want to remerge. Opening existing file now")
        return Table.read(data_dir + f'/merged/{tract}_{survey}_{euclid_survey}_merged.parquet')
    print('Check tells me to start the merge, starting now')

    lsst_ra, lsst_dec = lsst_table['coord_ra'], lsst_table['coord_dec']

    NSIDE=4096
    ## get the unique pixels of LSST data
    lsst_upix4096 = np.unique(hp.ang2pix(NSIDE, lsst_ra, lsst_dec, lonlat=True), return_counts=False)
    ## then get the pixels of Euclid data
    euclid_pix4096 = hp.ang2pix(NSIDE, euclid_table['right_ascension'], euclid_table['declination'], lonlat=True)
    ## Euclid has more coverage right now. We only keep the sources that lie in the LSST coverage
    spatial_mask = np.isin(euclid_pix4096, lsst_upix4096) #[lsst_cts > 8])
    euclid_field = euclid_table[spatial_mask]
    euclid_ra, euclid_dec = euclid_field['right_ascension'], euclid_field['declination']
    
    del NSIDE, lsst_upix4096, euclid_pix4096, spatial_mask, euclid_table
    gc.collect()
    
    ## match() is from ugali tools -- matching LSST and Euclid sources
    if len(euclid_ra) == 0:
        return 0
    indexlsst, indexeuclid, ds = match(lsst_ra, lsst_dec, euclid_ra, euclid_dec, tol = 0.0003)
    #print('index lsst:', '\n', indexlsst[0:20])
    #print('index euclid:', '\n', indexeuclid[0:20])
    matches_lsst = lsst_table[indexlsst]
    unmatched_lsst = lsst_table[~indexlsst]
    #print(matches_lsst.columns)
    matches_euclid = euclid_field[indexeuclid]
    unmatched_euclid = euclid_field[~indexeuclid]
    if len(matches_lsst) != len(matches_euclid):
        print("Something isn't right: those lengths don't match")
    del indexlsst, indexeuclid, lsst_ra, lsst_dec, euclid_ra, euclid_dec
    gc.collect()

    ## now merging our matches into one catalog with all LSST and Euclid columns
    matches_lsst['_match_id'] = np.arange(len(matches_lsst))
    matches_euclid['_match_id'] = np.arange(len(matches_euclid))
    merged_table = join(matches_lsst, matches_euclid, keys='_match_id')
    if not os.path.exists(data_dir + f'/merged'):
        os.mkdir(data_dir + f'/merged')
    merged_table.write(data_dir + f'/merged/{tract}_{survey}_{euclid_survey}_merged.parquet',
                       format='parquet', overwrite = True)

    if validation_needed==True:
        plotting_functions.match_validation_plots(tract, survey, euclid_survey, 
                                                  merged_table, matches_lsst, matches_euclid,
                                                  unmatched_lsst, unmatched_euclid,
                                                  lsst_table, euclid_field, ds)
    del matches_lsst, matches_euclid, unmatched_lsst, unmatched_euclid, lsst_table, euclid_field, ds
    gc.collect()
    return merged_table
