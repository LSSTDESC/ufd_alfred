import sys
import yaml
import os
import gc
import numpy as np
from astropy.coordinates import SkyCoord
from lsst.daf.butler import Butler
from simple_adl import coordinate_tools
from alfred import utils, DataObjects, RegionObjects, merging_catalogs, masks_and_filters, search_tools, plotting_functions, mapmaking

#-----------------------------------

# Note: trying to keep everything generalized as Optical and IR survey as this evolves from DES+Euclid -> Rubin+Roman

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

    nside = cfg['nside']

    opt_survey = cfg['opt_survey']
    ir_survey = cfg['ir_survey']
    repo_config = cfg[opt_survey]['repo_config'][where]
    collection = cfg[opt_survey]['collection'][where]
    #tract_list = cfg[survey]['tract_list']
    INCOLS1 = cfg[opt_survey]['INCOLS']
    opt_bands = cfg[opt_survey]['bands']
    INCOLS2 = cfg[ir_survey]['INCOLS']
    ir_bands = cfg[ir_survey]['bands']

## Define Which Area We're Looking At -- by nside and by coordinate
coord = (53.16, -28.10) #just using ECDFS center for now
SearchRegion = RegionObjects.Region(nside, coord)

## Define Which Columns to Pull Up From Data
opt_INCOLS = utils.columns_to_query(INCOLS1, opt_bands, output_type='list')
ir_INCOLS = utils.columns_to_query(INCOLS2, ir_bands, output_type = 'string')


## Load Rubin Data -- means repo_config and collection are defined and we can use butler
if repo_config != 0 and collection != 0:
    ## Initiate the Butler Instance
    butler = Butler(repo_config, collections=collection)
    #SkyMap =  butler.get('skyMap', skymap=skymap, collections=collection)
    tract_arr = SearchRegion.get_rubin_tracts(butler)
    OptData = SearchRegion.rubin_query(butler, tract_arr, opt_INCOLS)
## Load Euclid Data
if 'euclid' in ir_survey:
    IRData = SearchRegion.euclid_query(ir_INCOLS, preload = True)
# the queries return the data as their respective objects

## Merge Catalogs and Clean Up Memory
mergedData_raw = merging_catalogs.merge_catalogs(OptData, IRData, SearchRegion,
                                                 preload = True, validation_needed = False)
print('Merging catalogs completed')
del OptData, IRData #have to think if I'll need these again, can perhaps save them in Region obj
gc.collect()

## Clean Up Quality -- this will depend on which surveys are being used
if 'lsst' in opt_survey:
    # Q: which band snr should I enforce? - right now doing really lax snr > 3 cut
    snr_mask = masks_and_filters.clean_snr(mergedData_raw.g.mag, mergedData_raw.g.magerr, 3)
    snr_mask &= masks_and_filters.clean_snr(mergedData_raw.z.mag, mergedData_raw.z.magerr, 3)
    # this enforces that there are no per band flux flags
    opt_flag_mask = masks_and_filters.clean_lsst(mergedData_raw.data, 'griz')
if 'euclid' in ir_survey:
    snr_mask &= masks_and_filters.clean_snr(mergedData_raw.VIS.mag, mergedData_raw.VIS.magerr, 3)
    # Q: which euclid flags to enforce?
    # 0=no flags, 8=source close to a border, 512=source within an extended object area
    ir_flag_mask = masks_and_filters.clean_euclid(mergedData_raw.data, [0,8,512])

## mix em all together
total_mask = snr_mask & opt_flag_mask & ir_flag_mask
## clean up data
mergedData = mergedData_raw.apply_mask(total_mask)
SearchRegion.data_dict[opt_survey+'-'+ir_survey] = mergedData
print('Data cleaned and stored')

## Select for Stars -- Zerjal + colorcolor Cuts
colorcolor_mask = masks_and_filters.niroptical_color_stars(mergedData)
morphology_mask = masks_and_filters.Zerjal_stars(mergedData)
morphncolor_mask = colorcolor_mask & morphology_mask
# no one cared who I was til I put on the mask
stars = mergedData.apply_mask(morphncolor_mask)
SearchRegion.data_dict['stellar catalog'] = stars
SearchRegion.data = stars
print('Stellar catalog made')

## Some S-G validation plots - NEEDS TO BE UPDATED
"""
plotting_functions.color_magnitude(stars.g_mag, 'g', stars.r_mag, 'r', 
                                   'c',
                                   f'''g vs g-r of {stars.lsst_survey} and {stars.euclid_survey} stars 
                                   in tract {stars.tract}''',
                                   histogram = False,
                                   selection_label = 'Zerjal morphology + colorcolor cut',
                                   save = True, filename = f'{stars.tract}_{stars.lsst_survey}_{stars.euclid_survey}')
plotting_functions.color_color([('g', stars.g_mag),('r', stars.r_mag),('r', stars.r_mag),('i', stars.i_mag)],
                               None, 
                               f'''g-r vs r-i of {stars.lsst_survey} and {stars.euclid_survey} stars in tract
                               {stars.tract}''',
                               histogram = True, 
                               x_lim = (0,2), y_lim = (-1,2),
                               selection_label = 'Zerjal morphology + colorcolor cut',
                               save = True, filename = f'{stars.tract}_{stars.lsst_survey}_{stars.euclid_survey}')
plotting_functions.star_gal_sep(merged_data.i_mag, merged_data.mumax_minus_mag, 'Euclid mumax_minus_mag',
                                merged_data.pointlikeprob, 
                                f'''morphology separation of of {stars.lsst_survey} and {stars.euclid_survey} stars
                                in tract {stars.tract}''', 
                                histogram = True,
                                save = True, filename = f'{stars.tract}_{stars.lsst_survey}_{stars.euclid_survey}')
print('S-G plots ran and saved')
"""
## hotspot search
#distance_array=np.arange(50,1000,50) #distance is given in kpc
distance_array = [400]
for distance in distance_array:
    distance_modulus = coordinate_tools.distanceToDistanceModulus(distance)
    iso_sel, iso_stars = search_tools.isochrone_search(stars.g, stars.r, 
                                                       distance_modulus, stars,
                                                       age=12.0, Z=0.0002, 
                                                       save_graph=False)
    results = search_tools.search_by_distance(stars.survey, SearchRegion, distance_modulus, iso_sel)
            #survey isn't actually used? so maybe just str?
            #region is an object, I think I've added all the attributes and methods necessary to use their functions
    print(results)
'''
    ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array = np.asarray(results)
    if len(results[3]) == 0:
        return
    best_ra_peak, best_dec_peak, best_r_peak, best_distance_modulus, n_obs_peak, n_obs_half_peak, n_model_peak, best_sig_peak = 0, 0, 0, 0, 0, 0, 0, 0
    for i in range(len(results[0])): #this is how long the array is
'''        

## need fracdet eventually, but not prioritizing for now
# maybe put these functions as methods of region object
#full_map = mapmaking.euclid_fullmap('q1.vmpz_healpix_coverage', 'vis', 'coverage', preload=True)
#masked_map = mapmaking.match_map_polygon(full_map, tract.corners)
#tract_map,fracdet_map = mapmaking.rubin_maps(butler, tract.tract, 
#                                             map_name = 'deepCoadd_psf_maglim_map_weighted_mean', band = 'i', 
#                                             nside=2048, 
#                                             save_plot=True, map_title = f'Tract {tract.tract} Rubin i MagLim Map')
#plotting_functions.map_plot(full_map, 'Full Euclid VIS Coverage Map', color_lims = (24,26), 
#         save = True, filename = 'full_coverage_vis_q1')
#plotting_functions.map_plot(masked_map, f'Tract {tract.tract} Euclid VIS Coverage Map', color_lims = (24,26), 
#         save = True, filename = f'{tract.tract}masked_coverage_vis_{euclid_survey}')
#plotting_functions.map_plot(tract_map, f'Tract {tract.tract} Rubin i MagLim Map', color_lims = (24,26), 
#         save = True, filename = f'{tract.tract}_maglim_i_{survey}')
