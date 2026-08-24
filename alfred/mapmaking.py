from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import yaml
import os

from external.ugali.utils import healpix
import healsparse as hsp
import healpy as hp
from hpgeom import hpgeom
import skyproj
from lsst.daf.butler import Butler

from alfred import utils, plotting_functions

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    #assuming that it's cool that the whole github repo is considered "home"
    where = cfg['setup']['where']
    home_dir = os.path.expandvars(cfg['setup']['home_dir'][where])
    #external data is gonna be in a directory above - subject to change
    data_dir = os.path.join(home_dir, cfg['setup']['data_dir'])
    if not os.path.exists(data_dir+'/maps'):
        os.mkdir(data_dir+'/maps')
    
    survey = cfg['survey']
    collection = cfg[survey]['collection'][where]
    skymap = cfg[survey]['skymap']
    
    euclid_survey = cfg['euclid_survey']

def rubin_maps(butler, tract, 
               map_name = 'deepCoadd_psf_maglim_map_weighted_mean', band = 'i', nside=2048, 
               save_plot=True, map_title = '', fracdet_title = ''):
    #learning this lives in scratch/healsparse.ipynb -- go back there if needing to change something
    tract_map =  butler.get(map_name, band = band,
                            collections = collection, skymap = skymap,
                            tract = tract,
                            parameters={"degrade_nside": nside},)
    fracdet = tract_map.fracdet_map(nside) #coverage_mask ??

    if save_plot == True:
        if map_title == '':
            map_title = f'Tract {tract}, {survey} \n {map_name}'
        if fracdet_title = '':
            fracdet_title = f'Tract {tract}, {survey} FracDet of \n {map_name}'
        map_filename = f'{tract}_{survey}_{map_name}'
        fracdet_filename = f'{tract}_{survey}_FRACDET_{map_name}'
        plotting_functions.map_plot(tract_map, map_title, save = True, filename = map_filename)
        plotting_functions.map_plot(fracdet, fracdet_title, save = True, filename = filename)
        
    return tract_map, fracdet

def map_query(map_name):
    query = f"SELECT * FROM {map_name}"
    results = Euclid.launch_job_async(query).get_results()
    results.sort('file_path')
    return results

def make_euclid_bigmap(map_name, band, preload=True):
    '''
    Saves the whole sky Euclid map as a Healsparse map
    
    map_name suggestions: 'q1.vmpz_healpix_coverage', 'q1.vmpz_healpix_footprint_mask'
        code to get more Euclid map names is in scratch/euclid_astroquery.ipynb
    preload = True means that I want to use the preloaded / saved data instead of querying again
    '''
    query_check = check_if_query(data_dir + f'/maps/{map_name}_{euclid_survey}.fits',preload):
    if not query_check:
        print("Check tells me map exists and you don't want to overwrite. Doing nothing")
    if query_check:
        print("Check tells me map doesn't exists or you want to overwrite. Making map now")
        query = f"SELECT * FROM {map_name}"
        results = Euclid.launch_job_async(query).get_results()
        covmap = results.sort('file_path')

        tiles = [tile_list[0] for tile_list in np.unique(covmap['tile_index_list'])]
        map_list = []
        for tile in tiles:
            output_path = data_path+ f"/coverage_{tile}_vis.fits"
            i = np.where((covmap['tile_index_list']==[tile])&(covmap['filter_name']==band))[0][0]
            path = Euclid.get_product(file_name=covmap['file_name_list'][i], product_id = covmap['product_id'][i], 
                                          output_file= output_path, verbose=False)
            # ! need to figure out if this nside_coverage is telling hsp which coverage to read in 
            # or if it makes an assumption about the properties of this map
            map_list.append(hsp.HealSparseMap.read(output_path, nside_coverage=32))
            
        combined_map = hsp.operations.sum_union(map_list)
        combined_map.write(data_path+ f"/coverage_vis_combined.fits", clobber=False)

def map_tile_queryNsave(tile_id, map_name, band):
    '''
    Saves just the single tile of Euclid map
    '''
    query_check = check_if_query(data_dir + f'/maps/{tile_id}_{euclid_survey}_{map_name}.fits',preload):
    if not query_check:
        print("Check tells me map exists and you don't want to overwrite. Doing nothing")
    if query_check:
        print("Check tells me map doesn't exists or you want to overwrite. Making map now")
        query = f"SELECT * FROM {map_name}"
        results = Euclid.launch_job_async(query).get_results()
        covmap = results.sort('file_path')

        i_list = np.where((covmap['tile_index_list']==[tile_id])&(covmap['filter_name']==band))[0]
        n=1
        map_list = []
        for i in i_list:
            output_path = data_dir + f'/maps/{tile_id}_{euclid_survey}_{map_name}.fits'
            path = Euclid.get_product(file_name=covmap['file_name_list'][i], 
                                      product_id = covmap['product_id'][i], 
                                      output_file=output_path, verbose=False)
            # ! need to figure out if this nside_coverage is telling hsp which coverage to read in 
            # or if it makes an assumption about the properties of this map
            m = hp.read_map(output_path)
            nside = hp.get_nside(m)
            map_list.append(hsp.HealSparseMap.read(output_path, nside_coverage=32))
            n+=1
    return map_list