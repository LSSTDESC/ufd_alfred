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

import matplotlib.pyplot as plt
from matplotlib.patches import Annulus
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
    plots_dir = os.path.join(home_dir, cfg['output']['plots_dir'])
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)
    data_dir = os.path.join(os.path.expandvars(cfg['setup']['home_dir'][where]), cfg['setup']['data_dir'])
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

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
    OptData = SearchRegion.rubin_query(butler, INCOLS, tract_arr=tract_arr)
## or Load DES Data
elif 'des' in opt_survey:
    INCOLS = utils.columns_to_query(opt_INCOLS, opt_bands, output_type='string')
    OptData = SearchRegion.des_query(INCOLS, preload = True)
## Load Euclid Data
if 'euclid' in ir_survey:
    INCOLS = utils.columns_to_query(ir_INCOLS, ir_bands, output_type = 'string')
    IRData = SearchRegion.euclid_query(INCOLS, preload = True)
# the queries return the data as their respective objects

## Merge Catalogs and Clean Up Memory
mergedData_raw = merging_catalogs.merge_catalogs(OptData, IRData, SearchRegion,
                                                 preload = True, validation_needed = False) #can i plot buffer and inner outlines?
print('Merging catalogs completed')

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
print('Data cleaned and stored', len(mergedData.data), ' objects found')

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
print('Stellar catalog made, ', len(stars.data), ' stars found')


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
    iso_sel, iso_stars, iso = search_tools.isochrone_search(stars.g, stars.r, 
                                                       distance_modulus, stars,
                                                       SearchRegion,
                                                       age=12.0, Z=0.0002, 
                                                       save_graph=False)

    ## need fracdet eventually, but not prioritizing for now
    
    results, hist = search_tools.search_by_distance(stars.survey, SearchRegion, distance_modulus, iso_sel, 
                                                    verbose = False, plot = False) 
    
    one_peak_per_row = results.T
    peak_number = np.shape(one_peak_per_row)[0]
    if peak_number==0:
        continue
    for i in range(peak_number):
        #ra_peak, dec_peak, r_peak, sig_peak, distance_modulus, n_obs_peak, n_obs_half_peak, n_model_peak = results_transpose[i]
        Peak_i = DataObjects.Peak(one_peak_per_row[i], iso_stars, iso, SearchRegion)
        Peaks.append(Peak_i)
        
if len(Peaks)==0:
    print('No significant hotspots found.')
    utils.write_peak_result(Peaks, results_dir+f'/{SearchRegion.nside}_{SearchRegion.pixel}_{stars.survey}')

        
# for overlapping peaks, save the one with higher sig
print('Number of peaks before: ', len(Peaks))
# first, sort in order of significance
Peaks.sort(key=lambda x: x.sig, reverse=True)
utils.write_peak_result(Peaks, results_dir+f'/BEFOREsimple_{SearchRegion.nside}_{SearchRegion.pixel}_{stars.survey}', save_format='csv')


fig, ax = plt.subplots(1,2, figsize=(10,5))
overlaps = []
for i in range(len(Peaks)):
    radius = Annulus((Peaks[i].ra, Peaks[i].dec), Peaks[i].r, 0.0002, color='r')
    ax[0].add_artist(radius)
    ax[0].scatter(Peaks[i].ra, Peaks[i].dec, s=1, c='r')
    ax[0].set_title('before')
    for j in range(i+1, len(Peaks)):
        angsep = projector.angsep(Peaks[i].ra, Peaks[i].dec, Peaks[j].ra, Peaks[j].dec)
        if (angsep<=Peaks[i].r) | (angsep<=Peaks[j].r):
            overlaps.append((Peaks[i],Peaks[j]))
            
for pair in overlaps:
    more_sig = pair[0] #because in the pair, i always comes first. i < j and this list is sorted by significance
    less_sig = pair[1]
    if more_sig.sig==less_sig.sig: #I don't know what to do if they have the same significance yet, this should save both
        more_sig.distances.append(less_sig.distance)
    try:
        del Peaks[Peaks.index(less_sig)] #this (hopefully) works because they're all different objects 
                                #therefore none are exactly alike and they won't be deleted 2x
    #    print('successfully deleted above')
    except:
        continue
# as far as I can tell the above code works for duplicates (commented out code is for debugging)

print('Number of peaks after: ', len(Peaks))

for Peak in Peaks:
    radius = Annulus((Peak.ra, Peak.dec), Peak.r, 0.0002, color='r')
    ax[1].add_artist(radius)
    ax[1].scatter(Peak.ra, Peak.dec, s=1, c='r')
    ax[1].set_title('after')
#save duplicates diagnostic plot
plt.savefig(plots_dir + f'/{opt_survey}_{ir_survey}_nside{SearchRegion.nside}_pixel{SearchRegion.pixel}/peaks_duplicates_beforeafter.png')
plt.close()


final_Peaks = []
for Peak in Peaks:
    print(f'{Peak.sig} sigma; (RA, Dec, d) = ({Peak.ra} deg, {Peak.dec} deg, {Peak.distances} kpc); r = {Peak.r} deg; mu = {Peak.distance_modulus} mag')
    
    stars_within_radius = Peak.stars_within_the_radius(scale=1.2, set_attr = True)
    iso_background_stars = Peak.stars_within_the_radius(starsData=Peak.iso_stars, r=0.1, set_attr = False)
    all_background_stars = Peak.stars_within_the_radius(starsData=stars, r=0.1, set_attr = False)
    Peak.diagnostic_plots(iso_background_stars, all_background_stars,
                          plots_dir, data_dir, 
                          preload=True, save=True, cutouts=False)  
    
    if len(stars_within_radius.data)<=1:
        if type(stars_within_radius)==DataObjects.Data:
            print('Has type Data')
        elif type(stars_within_radius)==DataObjects.DESnEuclidData:
            print('Has type DESnEuclidData')
        print(' ')
        continue
    else:
        if type(stars_within_radius)==DataObjects.DESnEuclidData:
            print('Has type DESnEuclidData')
        final_Peaks.append(Peak)
        print(' ')
#save all the peaks that are more than 1 star...
utils.write_peak_result(final_Peaks, results_dir+f'/simple_{SearchRegion.nside}_{SearchRegion.pixel}_{stars.survey}', save_format='csv')


# now put those peaks into ugali...


## when I need maps:
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
