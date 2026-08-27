from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import yaml
import os

from external.ugali.utils import healpix
from external.astroquery_updated.esa.euclid import Euclid
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

## WHICH Rubin map do I want??

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
        if fracdet_title == '':
            fracdet_title = f'Tract {tract}, {survey} FracDet of \n {map_name}'
        map_filename = f'{tract}_{survey}_{map_name}'
        fracdet_filename = f'{tract}_{survey}_FRACDET_{map_name}'
        plotting_functions.map_plot(tract_map, map_title, save = True, filename = map_filename)
        plotting_functions.map_plot(fracdet, fracdet_title, save = True, filename = filename)
        
    return tract_map, fracdet

def euclid_map_query(map_name):
    '''
    map_name : string, the name of the table that I pulled from this code:
        tables = Euclid.load_tables(only_names=True, include_shared_tables=True)
        tile = [t.name for t in tables if 'q1' in t.name.lower()]
        ex. q1.vmpz_healpix_coverage

    returns table of the maps of that type, sorted by file_path 
        (file_path is where you'd eventually download them from)
    *THIS DOESN'T SAVE THE MAP*
    '''
    query = f"SELECT * FROM {map_name}"
    results = Euclid.launch_job_async(query).get_results()
    results.sort('file_path')
    return results

def euclid_fullmap(map_name, band, simple_name, preload=True):
    '''
    Opens or creates/saves the whole sky Euclid map as a Healsparse map
    
    map_name : str, ex. 'q1.vmpz_healpix_coverage', 'q1.vmpz_healpix_footprint_mask'
        code to get more Euclid map names is in scratch/euclid_astroquery.ipynb
    band : str, the filter of which map you want
    simple_name : str, name of the map for the file (i.e. coverage or footprint)
    preload = True means that I want to use the preloaded / saved data instead of querying again
    '''
    combined_map_path = data_dir + f'/maps/combined_{simple_name}_{band.lower()}_{euclid_survey}.fits'
    query_check = utils.check_if_query(combined_map_path, preload)
    if not query_check:
        print("Check tells me map exists and you don't want to overwrite. Opening map")
        combined_map = hsp.HealSparseMap.read(combined_map_path)
    if query_check:
        print("Check tells me map doesn't exists or you want to overwrite. Making map now")
        map_table = euclid_map_query(map_name)

        tiles = [tile_list[0] for tile_list in np.unique(map_table['tile_index_list'])]
        map_list = []
        for tile in tiles:
            output_path = data_path + f'/{tile}_{simple_name}_{band}_{euclid_survey}.fits'
            i = np.where((map_table['tile_index_list']==[tile]) & (map_table['filter_name']==band))[0][0]
            path = Euclid.get_product(file_name=map_table['file_name_list'][i], product_id = map_table['product_id'][i], 
                                          output_file= output_path, verbose=False)
            # ! need to figure out if this nside_coverage is telling hsp which coverage to read in 
            # or if it makes an assumption about the properties of this map
            map_list.append(hsp.HealSparseMap.read(output_path, nside_coverage=32))
            # am I going to want to delete the individual tile map file?
            
        combined_map = hsp.operations.sum_union(map_list)
        combined_map.write(combined_map_path, clobber=False)
        
    return combined_map

def euclid_tilemap(tile_id, map_name, band, simple_name, preload=True):
    '''
    Opens/makes and saves just the single tile of Euclid map
    
    tile_id : int, Euclid tile index (not tract)
    map_name : str, ex. 'q1.vmpz_healpix_coverage', 'q1.vmpz_healpix_footprint_mask'
        code to get more Euclid map names is in scratch/euclid_astroquery.ipynb
    band : str, filter of which map you want
    simple_name : str, name of the map for the file (i.e. coverage or footprint)
    preload = True means that I want to use the preloaded / saved data instead of querying again
    '''
    
    tile_map_path = data_path + f'/{tile}_{simple_name}_{band}_{euclid_survey}.fits'
    query_check = utils.check_if_query(tile_map_path,preload)
    if not query_check:
        print("Check tells me map exists and you don't want to overwrite. Opening map")
        tile_map = hsp.HealSparseMap.read(tile_map_path)
    if query_check:
        print("Check tells me map doesn't exists or you want to overwrite. Making map now")
        map_table = euclid_map_query(map_name)

        i_list = np.where((map_table['tile_index_list']==[tile_id])&(map_table['filter_name']==band))[0]
        if len(i_list) > 1:
            print('more than one map for this tile/filter choice')
        map_list = []
        # just in case there are multiple files for a tile/band (wouldn't think there would be?)
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
        tile_map = hsp.operations.sum_union(map_list)
        tile_map.write(tile_map_path, clobber=False)

    return tile_map

def match_map_polygon(fullmap,corners,nside=2048):
    '''
    takes in a healsparse map (I'm assuming it's going to be a Euclid map but generalizable) 
    and restricts it to be just the polygon area defined by corners (which I'm assuming is a Rubin tract)

    fullmap : Healsparse map
    corners : list of SkyCoord objects, should already be in order to create convex polygon
    '''
    ras = []
    decs = []
    for coord in corners:
        ras.append(coord.ra.value)
        decs.append(coord.dec.value)
    
    maskedmap = fullmap[hpgeom.query_polygon(nside, ras, decs, inclusive=True,)]

    return maskedmap