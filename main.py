import sys
import yaml
import os
import gc
from astropy.coordinates import SkyCoord
from lsst.daf.butler import Butler
from alfred import utils, DataObjects, RegionObjects, merging_catalogs

#-----------------------------------

#this is a bit hard coded too but idk another work around
#if main.py stays same folder as the config then this should work
with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    #assuming that it's cool that the whole github repo is considered "home"
    home_dir = os.path.expandvars(cfg['setup']['home_dir'])
    pckg_dir = os.path.join(home_dir, cfg['setup']['pckg_dir'])
    #external data is gonna be in a directory above - subject to change
    data_dir = os.path.expandvars(cfg['setup']['data_dir'])
    plots_dir = os.path.join(home_dir, cfg['output']['plots_dir'])
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)
    results_dir = os.path.join(home_dir, cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    
    survey = cfg['survey']
    euclid_survey = cfg['euclid_survey']
    repo_config = cfg[survey]['repo_config']
    collection = cfg[survey]['collection']
    tract_list = cfg[survey]['tract_list']
    INCOLS_addition = cfg[survey]['INCOLS_addition']

#initiate the butler instance
butler = Butler(repo_config, collections=collection)
#SkyMap =  butler.get('skyMap', skymap=skymap, collections=collection)

#define which columns to pull up from butler
lsst_INCOLS = [
    'coord_ra',
    'coord_dec',
    'detect_isIsolated',
]
bands='griz'
for band in bands:
    lsst_INCOLS += [f'{band}_psfFlux',
                    f'{band}_cModelFlux',
                    f'{band}_cModelFluxErr',
                    f'{band}_psfFluxErr',
                    f'{band}_extendedness',
                    f'{band}_psfFlux_flag'
    ]
    #depending on which survey, our preferred star-galaxy separating columns will change
    lsst_INCOLS += [col.replace('{band}',band) if '{band}' in col else col for col in INCOLS_addition]

# call up just one tract from EDFS for now
tract_num = 2394
tract = RegionObjects.Tract(tract_num, butler)
field = tract.field
lsst_table = tract.rubin_query(lsst_INCOLS)
      
#insert Euclid data load here
euclid_INCOLS = 'right_ascension, declination, point_like_prob, point_like_flag, ellipticity, mumax_minus_mag, flux_vis_psf, fluxerr_vis_psf, spurious_flag, det_quality_flag, fwhm, segmentation_map_id'
num = 2
for band in ['VIS', 'Y', 'J', 'H']:
    euclid_INCOLS += f", FLAG_{band}, FLUX_{band}_{num}FWHM_APER, FLUXERR_{band}_{num}FWHM_APER".lower()

euclid_table = tract.euclid_query(euclid_INCOLS, preload = False)

#merge catalogs and clean up memory
merged_table = merging_catalogs.merge_catalogs(lsst_table, euclid_table, tract.tract, 
                                               preload = False, validation_needed = False)
merged_data = DataObjects.LSSTnEuclidData(merged_table, survey, euclid_survey, field)
del lsst_table, euclid_table, merged_table
gc.collect()

#clean up quality

#select for stars -- Zerjal + colorcolor cuts

#isochrone_search

#get maps ready -- need fracdet

#compute_char_density