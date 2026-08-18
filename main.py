import sys
import yaml
import os
import gc
import numpy as np
from astropy.coordinates import SkyCoord
from lsst.daf.butler import Butler
from alfred import utils, DataObjects, RegionObjects, merging_catalogs, masks_and_filters, search_tools, plotting_functions

#-----------------------------------

## Set Environment Variables
with open('config.yaml', 'r') as ymlfile:
# this is a bit hard coded too but idk another work around
# if main.py stays same folder as the config then this should work
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    # assuming that it's cool that the whole github repo is considered "home"
    where = cfg['setup']['where']
    home_dir = os.path.expandvars(cfg['setup']['home_dir'][where])
    pckg_dir = os.path.join(home_dir, cfg['setup']['pckg_dir'])
    # external data is gonna be in a directory above - subject to change
    results_dir = os.path.join(home_dir, cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    survey = cfg['survey']
    euclid_survey = cfg['euclid_survey']
    repo_config = cfg[survey]['repo_config']
    collection = cfg[survey]['collection']
    #tract_list = cfg[survey]['tract_list']
    INCOLS_addition = cfg[survey]['INCOLS_addition']

## Initiate the Butler Instance
butler = Butler(repo_config, collections=collection)
#SkyMap =  butler.get('skyMap', skymap=skymap, collections=collection)

## Define Which Columns to Pull Up From Butler
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
    # depending on which survey, our preferred star-galaxy separating columns will change
    lsst_INCOLS += [col.replace('{band}',band) if '{band}' in col else col for col in INCOLS_addition]

# Load Rubin Data Per Tract (just one from EDFS for now)
tract_num = 2394
tract = RegionObjects.Tract(tract_num, butler)
#field = tract.field
lsst_table = tract.rubin_query(lsst_INCOLS)
      
## Load Euclid Data
euclid_INCOLS = 'right_ascension, declination, point_like_prob, point_like_flag, ellipticity, mumax_minus_mag, flux_vis_psf, fluxerr_vis_psf, spurious_flag, det_quality_flag, fwhm, segmentation_map_id'
num = 2
for band in ['VIS', 'Y', 'J', 'H']:
    euclid_INCOLS += f", FLAG_{band}, FLUX_{band}_{num}FWHM_APER, FLUXERR_{band}_{num}FWHM_APER".lower()
euclid_table = tract.euclid_query(euclid_INCOLS, preload = True)

## Merge Catalogs and Clean Up Memory
merged_table = merging_catalogs.merge_catalogs(lsst_table, euclid_table, tract.tract, 
                                               preload = True, validation_needed = False)
merged_data_raw = DataObjects.LSSTnEuclidData(merged_table, survey, euclid_survey, tract.tract)
del lsst_table, euclid_table, merged_table
gc.collect()

## Clean Up Quality

# Q: which band snr should I enforce? - right now doing really lax snr > 3 cut
snr_mask = masks_and_filters.clean_snr(merged_data_raw.g_mag, merged_data_raw.g_magerr, 3)
snr_mask &= masks_and_filters.clean_snr(merged_data_raw.z_mag, merged_data_raw.z_magerr, 3)
snr_mask &= masks_and_filters.clean_snr(merged_data_raw.VIS_mag, merged_data_raw.VIS_magerr, 3)

# this enforces no per band flux flags
lsst_flag_mask = masks_and_filters.clean_lsst(merged_data_raw.data, 'griz')

# Q: which euclid flags to enforce?
# 0=no flags, 8=source close to a border, 512=source within an extended object area
euclid_flag_mask = masks_and_filters.clean_euclid(merged_data_raw.data, [0,8,512])

# mix em all together
total_mask = snr_mask & lsst_flag_mask & euclid_flag_mask
# cleaned up data
merged_data = merged_data_raw.apply_mask(total_mask)

## Select for Stars -- Zerjal + colorcolor Cuts
colorcolor_mask = masks_and_filters.niroptical_color_stars(merged_data)
morphology_mask = masks_and_filters.Zerjal_stars(merged_data)
morphncolor_mask = colorcolor_mask & morphology_mask
# no one cared who I was til I put on the mask
stellar_catalog = merged_data.apply_mask(morphncolor_mask)

## Isochrone Search
'''
nonans_mask = masks_and_filters.clean_nans(stellar_catalog.g_mag, stellar_catalog.g_magerr) & masks_and_filters.clean_nans(stellar_catalog.r_mag, stellar_catalog.r_magerr)
nonans_stellar_catalog = stellar_catalog.apply_mask(nonans_mask)
''' # I don't think applying this mask did much of anything
#distances = np.arange(200, 2000, 100)
distances = np.linspace(200,200,1)
isocut_stars_eachdistance = []
for distance in distances:
    isochrone_stars = search_tools.isochrone_search(stellar_catalog, distance, graph=True, save=True)
    isocut_stars_eachdistance.append(isochrone_stars)
isocut_stars_eachdistance_arr = np.array(isocut_stars_eachdistance)

## get maps ready -- need fracdet
## compute_char_density