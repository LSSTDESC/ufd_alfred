import yaml
import os
from astropy.table import Table
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy import units as u
import lsst.geom as geom
import numpy as np
import gc

from alfred import utils, DataObjects

from external.astroquery_updated.esa.euclid import Euclid
from external.ugali.utils import healpix
#I really need to fix this, maybe github submodules or enforcing a version of astroquery
#I think it's version 0.4.11 or 10?

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    where = cfg['setup']['where']
    survey = cfg['survey']
    skymap = cfg[survey]['skymap']
    repo_config = cfg[survey]['repo_config'][where]
    collection = cfg[survey]['collection'][where]
    data_dir = os.path.join(os.path.expandvars(cfg['setup']['home_dir'][where]), cfg['setup']['data_dir'])
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
    euclid_survey = cfg['euclid_survey']

class Region():
    #for sharding data into large healpy pixel regions
    def __init__(self, nside, pixel):
        self.nside = nside
        self.pixel = pixel
        ra,dec = hp.pix2ang(nside, pixel, lonlat=True)
        self.ra = ra
        self.dec = dec
        phi = healpix.lon2phi(ra)
        theta = healpix.lat2theta(dec)
        self.pixel_neighbors = hp.pixelfunc.get_all_neighbours(nside, theta, phi=phi) 
            #returns 8 nearest pixel indices 
        #self.borders = hp.vec2ang(hp.boundaries(nside, pixel, step=1, nest=False).T, lonlat=True)
        #self.borders_str = 

        #getting these overlapping regions for data querying purposes
        self.rubin_tracts = -1
        self.euclid_tiles = -1
        self.des_tiles = -1

        self.rubin_data = 0
        self.euclid_data = 0
        self.des_data = 0
        self.data = {}

    def region_cut(self, Data, finer_nside = 4096):
        '''
        input a Data Object that has an ra and dec attribute, needs to be in degrees, and apply_mask method
        '''
        subpix = healpix.subpixel(self.pixel, self.nside, finer_nside)
        datapix = healpix.ang2pix(finer_nside, Data.ra, Data.dec)
        mask = np.isin(datapix, subpix)
        region_data = Data.apply_mask(mask)
        return region_data

    def region_borders(self, return_type = 'string', step = 1):
        ra_arr, dec_arr = hp.vec2ang(hp.boundaries(self.nside, self.pixel, step=step, nest=False).T, lonlat=True)
        if return_type == 'string':
            border_str = ''
            for ra, dec in zip(ra_arr, dec_arr):
                border_str += f'{ra}, {dec}, '
            border_str = border_str.removesuffix(', ')
            return border_str
        elif return_type == 'list':
            border_list = []
            for ra, dec in zip(ra_arr, dec_arr):
                border_list.append((ra, dec))
            return border_list
            
        
    def get_rubin_tracts(self, butler):
        '''
        SkyMap = skyMap object, generated from butler
        '''
        
        SkyMap =  butler.get('skyMap', skymap=skymap, collections=collection)
        ra_arr, dec_arr = self.borders
        tract_ids = []
        for ra,dec in zip(ra_arr,dec_arr):
            TractsInfo = SkyMap.findTract(geom.SpherePoint(ra*geom.degrees,dec*geom.degrees))
            tract_ids.append(TractsInfo.tract_id)
        tract_ids = np.unique(tract_ids)
        
        self.rubin_tracts = tract_ids #maybe have them as Tract objects?

        return tract_ids

    def rubin_query(self, butler, tract_arr, INCOLS):
        '''
        I don't have a preload because I figured that we'd just be using the butler and not saving
        
        queries by tract but for a whole array of tracts, then restricts based on a healpix mask 
        of what is actually in that region
        '''
        rubin_data_list = [] #these are a bunch of LSSTData objects -- need to think how to concatenate
        for tract in tract_arr:
            full_tract = butler.get('object', 
                                    dataId={'skymap': skymap, 'tract': tract}, 
                                    collections=collection, parameters={"columns":INCOLS})
            rubin_data_list.append(self.region_cut(LSSTData(full_tract, survey)).data)
                #region_cut returns an LSSTData object, so append just the data of that object
            del full_tract
            gc.collect()
        rubin_data = LSSTData(np.vstack(rubin_data_list), survey)
            #shove all the tract data together, make it an LSST object
        self.rubin_data = rubin_data
        self.data['LSST '+ survey.upper()] = rubin_data
        return rubin_data

    def euclid_query(self, INCOLS, preload = True):
        '''
        checks if the data already exists for that survey and nside/pixel. 
        if it does, we won't query again, if we query, this function saves the new queried data
            preload means "yes I want you to load the data if it exists"
            might be False if you want to overwrite the data
        either way, this assigns the data to self.euclid_data and self.data['Euclid + survey'] 
        '''
        if not os.path.exists(data_dir + f'/{euclid_survey}'):
            print("no Euclid data folder, making one now")
            os.mkdir(data_dir + f'/{euclid_survey}')
        
        file_dir = data_dir + f'/{euclid_survey}/{self.nside}_{self.pixel}_euclid.parquet'
        if not utils.check_if_query(file_dir, preload):
            print("Check tells me Euclid data exists and you don't want to overwrite. Opening existing file now")
            results = Table.read(file_dir)
        else:
            print("Check tells me Euclid data doesn't exist or you do want to overwrite, querying now")
    
            query = f'SELECT {INCOLS} FROM mer_catalogue'
            radius = 1.7 #going for bigger than a tract
            #query += f''' WHERE DISTANCE({self.center.ra.value}, {self.center.dec.value},
            #                            right_ascension, declination) < {radius}'''
            query += f''' WHERE CONTAINS(POINT('ICRS', right_ascension, declination),
                         POLYGON('ICRS', {self.region_borders(return_type='string',step=1)})) = 1'''
    
            results_table = Euclid.launch_job_async(query, verbose=False).get_results()
            if not os.path.exists(data_dir + f'/{euclid_survey}'):
                os.mkdir(data_dir + f'/{euclid_survey}')
            results_table.write(data_dir + f'/{euclid_survey}/{self.tract}_euclid.parquet',
                                   format='parquet', overwrite = True)
        
        results = EuclidData(results_table, euclid_survey)
        self.euclid_data = results
        self.data['Euclid '+ euclid_survey.upper()] = results
        
        return results

    def des_query():
        '''
        thinking that this would be really similar to the Euclid function in form
        taking in the columns, using the region borders, returning results/updating attributes

        self.des_data = query_result
        self.data[DES survey] = query_result
        '''



class Tract():
    def __init__(self, tract, butler):
        '''
        tract = int
        SkyMap = skyMap object, generated from butler
        '''
        self.tract = tract
        self.field = utils.get_field(tract)

        self.butler = butler
        SkyMap =  butler.get('skyMap', skymap=skymap, collections=collection)

        self.center = SkyCoord(SkyMap.generateTract(tract).getCtrCoord().getRa().asDegrees()*u.deg, 
                               SkyMap.generateTract(tract).getCtrCoord().getDec().asDegrees()*u.deg, 
                               frame='icrs')
        self.center_SpherePoint = SkyMap.generateTract(tract).getCtrCoord()
        ras = sorted([SkyMap.getRaDecRange(tract)[0].asDegrees()*u.deg, 
                      SkyMap.getRaDecRange(tract)[1].asDegrees()*u.deg])
        self.ra_range = ras
        decs = sorted([SkyMap.getRaDecRange(tract)[2].asDegrees()*u.deg, 
                       SkyMap.getRaDecRange(tract)[3].asDegrees()*u.deg])
        self.dec_range = decs
        self.corners = [SkyCoord(ra,decs[1],frame='icrs') for ra in ras] + [SkyCoord(ra,decs[0],frame='icrs') for ra in ras[::-1]]
        corners_str = ''
        for coord in self.corners:
            corners_str += f'{coord.ra.value}, {coord.dec.value}, '
        self.corners_str = corners_str.removesuffix(', ')
        self.corners_Angle = SkyMap.getRaDecRange(tract)


# I don't know if I'll want to do a patch class to make these smaller than tract
class Patch(Tract):
    def __init__(self, patch_num):
        self.patch = patch_num
        self.patch_data = self.tract_data[self.tract_data['patch']==patch_num]


