import yaml
import os
from lsst.daf.butler import Butler
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
    #also have skymap, instrument, tract_list, tract2field_dict information

butler = Butler(repo_config, collections=collection)
if butler is not None:
    print("yay!")
