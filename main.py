import yaml
import os
import gc
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
import pandas as pd
from ugali.utils import projector
from alfred import utils, DataObjects, RegionObjects, merging_catalogs, masks_and_filters, search_tools, plotting_functions, mapmaking
import itertools
import pdb
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
    if 'lsst' in opt_survey:
        repo_config = cfg[opt_survey]['repo_config'][where]
        collection = cfg[opt_survey]['collection'][where]
    #tract_list = cfg[survey]['tract_list']
    opt_INCOLS = cfg[opt_survey]['INCOLS']
    opt_bands = cfg[opt_survey]['bands']
    ir_INCOLS = cfg[ir_survey]['INCOLS']
    ir_bands = cfg[ir_survey]['bands']

## Define Which Area We're Looking At -- by nside and by coordinate
coord = (53.16, -28.10) #just using ECDFS center for now
SearchRegion = RegionObjects.Region(nside, coord)
corners_ra, corners_dec = SearchRegion.corners
buffer_ra, buffer_dec = SearchRegion.buffer_region(corners_ra, corners_dec, SearchRegion.coord_center)

## Load Rubin Data
if 'lsst' in opt_survey:
    from lsst.daf.butler import Butler
    INCOLS = utils.columns_to_query(opt_INCOLS, opt_bands, output_type='list')
    ## Initiate the Butler Instance
    butler = Butler(repo_config, collections=collection)
    #SkyMap =  butler.get('skyMap', skymap=skymap, collections=collection)
    tract_arr = SearchRegion.get_rubin_tracts(butler)
    OptData = SearchRegion.rubin_query(butler, tract_arr, INCOLS)
## or Load DES Data
elif 'des' in opt_survey:
    INCOLS = utils.columns_to_query(opt_INCOLS, opt_bands, output_type='string')
    OptData = SearchRegion.des_query(INCOLS, preload=True)
## Load Euclid Data
if 'euclid' in ir_survey:
    INCOLS = utils.columns_to_query(ir_INCOLS, ir_bands, output_type = 'string')
    IRData = SearchRegion.euclid_query(INCOLS, preload = True)
# the queries return the data as their respective objects

## Merge Catalogs and Clean Up Memory
mergedData_raw = merging_catalogs.merge_catalogs(OptData, IRData, SearchRegion,
                                                 preload = True, validation_needed = False)
print('Merging catalogs completed')

#was having difficulty with masked arrays
'''
results = mergedData_raw.data
print(type(results))
for col in results.colnames:
    print(col, type(results[col]))
'''   

del OptData, IRData #have to think if I'll need these again, can perhaps save them in Region obj
gc.collect()

## Clean Up Quality -- this will depend on which surveys are being used
    # Q: which band snr should I enforce? - right now doing really lax snr > 3 cut
snr_mask = masks_and_filters.clean_snr(mergedData_raw.g.mag, mergedData_raw.g.magerr, 3)
snr_mask &= masks_and_filters.clean_snr(mergedData_raw.z.mag, mergedData_raw.z.magerr, 3)
snr_mask &= masks_and_filters.clean_snr(mergedData_raw.VIS.mag, mergedData_raw.VIS.magerr, 3)
    # Q: which flags should I enforce?
    # 0=no flags, 8=source close to a border, 512=source within an extended object area
flag_mask = mergedData_raw.clean(lsst_bands='griz', euclid_flags=[0,8,512])
## mix em together
total_mask = snr_mask & flag_mask
## clean up data
mergedData = mergedData_raw.apply_mask(total_mask)
SearchRegion.data_dict[opt_survey+'-'+ir_survey] = mergedData
print('Data cleaned and stored')

## Select for Stars -- Zerjal + colorcolor Cuts
colorcolor_mask = masks_and_filters.niroptical_color_stars(mergedData)
morphology_mask = masks_and_filters.Zerjal_stars(mergedData)
morphncolor_mask = colorcolor_mask & morphology_mask

# the masked arrays need to be filled to use simple's cut_isochrone_path
stars_table = mergedData.data[morphncolor_mask]
stars_table_filled = utils.handle_ma_arr(stars_table, solution='fill')
stars = type(mergedData)(stars_table_filled, survey=mergedData.survey, coord_choice=mergedData.coord_choice)
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
distance_array = np.arange(50,1000,50) #distance is given in kpc

Peaks = []

for distance in distance_array:
    print(f'Searching at {distance} kpc')
    distance_modulus = projector.distanceToDistanceModulus(distance)
    ## NOTE: need new des x euclid isochrone
    iso_sel, iso_stars, iso = search_tools.isochrone_search(stars.g, stars.r, 
                                                       distance_modulus, stars,
                                                       SearchRegion,
                                                       age=12.0, Z=0.0002, 
                                                       save_graph=False)
    save=False
    ## need fracdet eventually, but not prioritizing for now
    
    results = np.asarray(search_tools.search_by_distance(stars.survey, SearchRegion, distance_modulus, iso_sel, verbose = False)) #survey isn't actually used in this function it seems? so just putting in a str...?
    
    one_peak_per_row = results.T
    peak_number = np.shape(one_peak_per_row)[0]
    if peak_number==0:
        continue
    for i in range(peak_number):
        #ra_peak, dec_peak, r_peak, sig_peak, distance_modulus, n_obs_peak, n_obs_half_peak, n_model_peak = results_transpose[i]        
        Peaks.append(DataObjects.Peak(one_peak_per_row[i], iso_stars, iso, SearchRegion))
        
if len(Peaks)==0:
    print('No significant hotspots found.')
    utils.write_peak_result(Peaks, results_dir+f'/{SearchRegion.nside}_{SearchRegion.pixel}_{stars.survey}')
        
# for overlapping peaks, save the one with higher sig
'''
moresig_Peaks = []
for Peak1,Peak2 in itertools.combinations(Peaks,2):
    angsep = projector.angsep(Peak1.ra, Peak1.dec, Peak2.ra, Peak2.dec)
    if angsep<Peak1.r:
        two_peaks_sig = np.array([Peak1.sig, Peak2.sig])
        two_peaks = [Peak1,Peak2]
        moresig_Peak = two_peaks[np.argmax(two_peaks_sig)]
        moresig_Peaks.append(moresig_Peak)
        print('peaks 1 and 2: ', Peak1.id, ' | ', Peak2.id)
        print('peak being appended: ', moresig_Peak.id)
        print(' ')
        moresig_Peak.overlapping_peaks.append(two_peaks[np.argmin(two_peaks_sig)])
'''
moresig_Peaks = []
for i in range(len(Peaks)):
    overlaps = []
    for j in range(len(Peaks)):
        if i==j:
            continue
        else:
            angsep = projector.angsep(Peaks[i].ra, Peaks[i].dec, Peaks[j].ra, Peaks[j].dec)
            if angsep<Peaks[i].r:
                overlaps.append(Peaks[j])
    try:
        mostsig_Peak_i = max(overlaps, key=lambda x: x.sig)
        moresig_Peaks.append(mostsig_Peak_i)
    except:
        moresig_Peaks.append(Peaks[i])
#the combination iteration gets duplicates, so drop duplicates (can probably make this logic better...)
moresig_Peaks = list(set(moresig_Peaks))
# sort in order of significance
moresig_Peaks.sort(key=lambda x: x.sig, reverse=True)
            
for Peak in moresig_Peaks:
    print(f'{Peak.sig} sigma; (RA, Dec, d) = ({Peak.ra} deg, {Peak.dec} deg, {Peak.distance} kpc); r = {Peak.r} deg; mu = {Peak.distance_modulus} mag')
utils.write_peak_result(moresig_Peaks, results_dir+f'/{SearchRegion.nside}_{SearchRegion.pixel}_{stars.survey}', save_format='csv')

    
'''
ra_peak_list = []
dec_peak_list = [] 
r_peak_list = []
sig_peak_list = []
distance_modulus_list = []
mc_source_id_list = []
n_obs_peak_list = []
n_obs_half_peak_list = []
n_model_peak_list = []
        ra_peak_list.append(ra_peak)
        dec_peak_list.append(dec_peak)
        r_peak_list.append(r_peak)
        sig_peak_list.append(sig_peak)
        distance_modulus_list.append(distance_modulus)
        n_obs_peak_list.append(n_obs_peak)
        n_obs_half_peak_list.append(n_obs_half_peak)
        n_model_peak_list.append(n_model_peak)
        #mc_source_id_list.append(np.tile(0, len(sig_peaks))) <- what is this? I made the logic moot by appending 1 peak at a time rather than per distance
=
    ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array = np.asarray(results)
    if len(results[3]) == 0:
        return
    best_ra_peak, best_dec_peak, best_r_peak, best_distance_modulus, n_obs_peak, n_obs_half_peak, n_model_peak, best_sig_peak = 0, 0, 0, 0, 0, 0, 0, 0
    for i in range(len(results[0])): #this is how long the array is
'''        

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
