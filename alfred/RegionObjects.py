from astropy.coordinates import SkyCoord
from astropy import units as u
import lsst.geom as geom
import yaml
import os
import sys
from astropy.table import Table
from alfred import utils
#first need to load in the module I need
#Goes up a directory to get the updated astroquery
#sys.path.append(os.path.abspath('../..'))
from external.astroquery_updated.esa.euclid import Euclid
#I really need to fix this, maybe github submodules or enforcing a version of astroquery
#I think it's version 0.4.11 or 10?

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    survey = cfg['survey']
    skymap = cfg[survey]['skymap']
    repo_config = cfg[survey]['repo_config']
    collection = cfg[survey]['collection']
    where = cfg['setup']['where']
    data_dir = os.path.join(os.path.expandvars(cfg['setup']['home_dir'][where]), cfg['setup']['data_dir'])
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
    euclid_survey = cfg['euclid_survey']

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
        ras = [SkyMap.getRaDecRange(tract)[0].asDegrees()*u.deg, SkyMap.getRaDecRange(tract)[1].asDegrees()*u.deg]
        self.ra_range = ras
        decs = [SkyMap.getRaDecRange(tract)[2].asDegrees()*u.deg, SkyMap.getRaDecRange(tract)[3].asDegrees()*u.deg]
        self.dec_range = decs
        self.corners = [SkyCoord(ra,dec,frame='icrs') for ra in ras for dec in decs]
        corners_str = ''
        for ra in ras:
            for dec in decs:
                corners_str += f'{ra.value}, {dec.value}, '
        self.corners_str = corners_str.removesuffix(', ')
        self.corners_Angle = SkyMap.getRaDecRange(tract)

        self.tract_data = 0

    def rubin_query(self, INCOLS):
        '''
        queries the tract, can restrict it to the patch in the child class if desired
        '''
        butler = self.butler
        full_tract = butler.get('object', 
                                dataId={'skymap': skymap, 'tract': self.tract}, 
                                collections=collection, parameters={"columns":INCOLS})
        self.tract_data = full_tract
        return full_tract

    def euclid_query(self, INCOLS, preload = True):
        if not os.path.exists(data_dir + f'/{euclid_survey}'):
            print("no Euclid data folder, making one now")
            os.mkdir(data_dir + f'/{euclid_survey}')
        ## if we've already done this query, just load in that data
            ## (unless user wants to override that for overwriting purposes)
        if not os.path.exists(data_dir + f'/{euclid_survey}/{self.tract}_euclid.parquet'):
            print("Euclid data doesn't exist")
            preload = False

        if preload == True:
            print("Euclid data exists and you don't want to overwrite existing data, opening file now")
            results_table = Table.read(data_dir + f'/{euclid_survey}/{self.tract}_euclid.parquet')
        else:
            print("Euclid data doesn't exist or you want to overwrite existing data, querying now")
            query = f'SELECT {INCOLS} FROM mer_catalogue'
            radius = 1.7 #going for bigger than a tract
            #query += f''' WHERE DISTANCE({self.center.ra.value}, {self.center.dec.value},
            #                            right_ascension, declination) < {radius}'''
            query += f''' WHERE CONTAINS(POINT('ICRS', right_ascension, declination),
                         POLYGON('ICRS', {self.corners_str})) = 1'''

            results_table = Euclid.launch_job_async(query, verbose=False).get_results()
            if not os.path.exists(data_dir + f'/{euclid_survey}'):
                os.mkdir(data_dir + f'/{euclid_survey}')
            results_table.write(data_dir + f'/{euclid_survey}/{self.tract}_euclid.parquet',
                                   format='parquet', overwrite = True)
        return results_table

# I don't know if I'll want to do a patch class to make these smaller than tract
class Patch(Tract):
    def __init__(self, patch_num):
        self.patch = patch_num
        self.patch_data = self.tract_data[self.tract_data['patch']==patch_num]


