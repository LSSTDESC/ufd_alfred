import healpy as hp
import numpy as np
import gc
import sys
import os
import yaml
from alfred import utils, plotting_functions
from astropy.table import Table, vstack, join
#Goes up a directories to get the updated astroquery
#I really need to fix this, maybe github submodules or enforcing a version of astroquery
#I think it's version 0.4.11 ?
#sys.path.append(os.path.abspath('../'))
from external.ugali.utils.projector import match
## function to create euclid + rubin datasets and register to the data registry on NERSC
## don't know where I want this to live quite yet

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    #assuming that it's cool that the whole github repo is considered "home"
    where = cfg['setup']['where']
    home_dir = os.path.expandvars(cfg['setup']['home_dir'][where])
    #external data is gonna be in a directory above - subject to change
    data_dir = os.path.join(home_dir, cfg['setup']['data_dir'])
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
    if not os.path.exists(data_dir + f'/merged'):
        os.mkdir(data_dir + f'/merged')
    results_dir = os.path.join(home_dir, cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)


## INSERT function to add it to the data registry

# function to merge catalogs and return the appropriate merged object
def merge_catalogs(PrimaryData, SecondaryData, SearchRegion, preload = True, validation_needed = False):
    # I'm thinking that Primary is the one you want to match to (maybe for reasons of less coverage)
    # don't want to assert yet that one is optical and one is ir
    
    file_path = data_dir + f'/merged/{SearchRegion.nside}_{SearchRegion.pixel}_{PrimaryData.release}_{SecondaryData.release}_merged.parquet'
    # function to check if the data doesn't exist already and if I want to rewrite it
    if not utils.check_if_query(file_path, preload):
        print("Check tells me data exists and you don't want to remerge. Opening existing file now")
        return Table.read(file_path)
    print('Check tells me to start the merge, starting now')

    prim_ra, prim_dec = PrimaryData.ra, PrimaryData.dec

    NSIDE=4096
    if SearchRegion.nside == NSIDE:
        print('Warning: might cause some issues that NSIDE is already smallest resolution possible. Make sure to check merge')
    ## get the unique pixels of primary dataset
    prim_upix4096 = np.unique(hp.ang2pix(NSIDE, prim_ra, prim_dec, lonlat=True), return_counts=False)
    ## then get the pixels of secondary data
    secun_pix4096 = hp.ang2pix(NSIDE, SecondaryData.ra, SecondaryData.dec, lonlat=True)
    ## We only keep the sources that lie in the Primary survey coverage
    spatial_mask = np.isin(secun_pix4096, prim_upix4096) #enforce [prim_cts > 8]) ?
    SecondaryData_masked = SecondaryData.apply_mask(spatial_mask)
    secun_ra, secun_dec = SecondaryData_masked.ra, SecondaryData_masked.dec
    
    del NSIDE, prim_upix4096, secun_pix4096, spatial_mask, SecondaryData
    gc.collect()
    
    ## match() is a spatial match from ugali tools, tol controls how generous you are in saying the sources overlap
    if len(secun_ra) == 0:
        print('uh oh, no overlap detected')
        return 0
    indexprim, indexsecun, ds = match(prim_ra, prim_dec, secun_ra, secun_dec, tol = 0.0003)
    matches_prim = PrimaryData.apply_mask(indexprim)
    unmatched_prim = PrimaryData.apply_mask(~indexprim)
    matches_secun = SecondaryData_masked.apply_mask(indexsecun)
    unmatched_secun = SecondaryData_masked.apply_mask(~indexsecun)
    if len(matches_prim) != len(matches_secun):
        print("Something isn't right: those lengths don't match")
    del indexprim, indexsecun, prim_ra, prim_dec, secun_ra, secun_dec
    gc.collect()

    ## now merging our matches into one catalog with all LSST and Euclid columns
    matches_prim['_match_id'] = np.arange(len(matches_prim))
    matches_secun['_match_id'] = np.arange(len(matches_secun))
    merged_table = join(matches_prim, matches_secun, keys='_match_id')
    merged_table.write(file_path, format='parquet', overwrite = True)

    if validation_needed==True:
        plotting_functions.match_validation_plots(tract, survey, euclid_survey, 
                                                  merged_table, matches_lsst, matches_euclid,
                                                  unmatched_lsst, unmatched_euclid,
                                                  lsst_table, euclid_field, ds)
    del matches_prim, matches_secun, unmatched_prim, unmatched_secun, PrimaryData, SecondaryData_masked, ds
    gc.collect()

    # I don't know how else to do this logic
    if 'lsst' in PrimaryData.release or 'lsst' in SecondaryData.release:
        if 'euclid' in PrimaryData.release or 'euclid' in SecondaryData.release:
            mergedData = LSSTnEuclidData(merged_table, PrimaryData.release, SecondaryData.release, coord_choice='LSST')
    # then more if statements for the other surveys... 
    del merged_table
    gc.collect()
    
    return mergedData
