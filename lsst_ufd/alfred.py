import yaml
import os
from lsst.daf.butler import Butler
from star_gal_sep import *
#-----------------------------------

#this is a bit hard coded too but idk another work around
#if alfred stays in a folder below the config always then this should work
with open('../config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)

    home_dir_str = cfg['setup']['home_dir']
    home_dir = os.path.expandvars(home_dir_str)
    pckg_dir = os.path.join(home_dir, cfg['setup']['pckg_dir'])
    results_dir = os.path.join(home_dir, cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    plots_dir = os.path.join(home_dir, cfg['output']['plots_dir'])
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)
    
    survey = cfg['survey']
    repo_config = cfg[survey]['repo_config']
    collection = cfg[survey]['collection']
    field2tract_dict = cfg[survey]['field2tract_dict']
    skymap = cfg[survey]['skymap']
    #also have instrument, tract_list, tract2field_dict information

#initiate the butler instance
butler = Butler(repo_config, collections=collection)

INCOLS = [
    'coord_ra',
    'coord_dec',
    'detect_isIsolated',
]
bands="griz"
for band in bands:
    INCOLS += [
        f'{band}_psfFlux',
        f'{band}_cModelFlux',
        f'{band}_cModelFluxErr',
        f'{band}_psfFluxErr',
        f'{band}_extendedness',
        f'{band}_psfFlux_flag'
    ]
    #if survey=='dp1':
    #    INCOLS += [f'{band}_SizeExtendedness']
    if survey=='dp2':
        INCOLS += [f'{band}_model_extendedness']
# maybe inelegant but that's a problem for future kayleigh
data = butler.get('object', collections=[collection],
                  dataId={'skymap': skymap, 'tract': 
                          field2tract_dict['EDFS'][0]}, 
                  parameters={"columns":INCOLS})
# just calling up one tract for now
# thank you alfred

#data = Data(survey, data_arr) 

stars = stellar_catalog(data, survey, 'fluxratioerr', 'i', 0.5, c=1.2)
