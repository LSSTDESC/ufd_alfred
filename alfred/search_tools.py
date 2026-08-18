import sys
#sys.path.append("../") #has to be run in same directory rn
import simple_adl.simple_adl.isochrone as isochrone
import simple_adl.simple_adl.coordinate_tools as coordinate_tools
from alfred import plotting_functions, utils
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import scipy
import yaml
import os

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    #assuming that it's cool that the whole github repo is considered "home"
    where = cfg['setup']['where']
    home_dir = os.path.expandvars(cfg['setup']['home_dir'][where])
    pckg_dir = os.path.join(home_dir, cfg['setup']['pckg_dir'])
    #external data is gonna be in a directory above - subject to change
    data_dir = os.path.expandvars(cfg['setup']['data_dir'])
    plots_dir = os.path.join(home_dir, cfg['output']['plots_dir'])
    if not os.path.exists(plots_dir+'/isochrones'):
        os.mkdir(plots_dir+'/isochrones')
    results_dir = os.path.join(home_dir, cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    survey = cfg['survey']
    euclid_survey = cfg['euclid_survey']

def cut_isochrone_path(g, r, g_err, r_err, isochrone, radius=0.01, mag_max = 26, return_all=False):
    """
    Cut to identify objects within isochrone cookie-cutter.
    """
    if np.all(isochrone.stage == 'Main'):
        # Dotter case
        index_transition = len(isochrone.stage)
    else:
        # Other cases
        index_transition = np.nonzero(isochrone.stage >= isochrone.hb_stage)[0][0] + 1    

    mag_1_rgb = isochrone.mag_1[0: index_transition] + isochrone.distance_modulus
    mag_2_rgb = isochrone.mag_2[0: index_transition] + isochrone.distance_modulus
    
    mag_1_rgb = mag_1_rgb[::-1]
    mag_2_rgb = mag_2_rgb[::-1]

    # Cut one way...
    f_isochrone = scipy.interpolate.interp1d(mag_2_rgb, mag_1_rgb - mag_2_rgb, bounds_error=False, fill_value = 999.)
    color_diff = np.fabs((g - r) - f_isochrone(r))
    cut_2 = (color_diff < np.sqrt(radius**2 + r_err**2 + g_err**2))

     # ...and now the other
    f_isochrone = scipy.interpolate.interp1d(mag_1_rgb, mag_1_rgb - mag_2_rgb, bounds_error=False, fill_value = 999.)
    color_diff = np.fabs((g - r) - f_isochrone(g))
    cut_1 = (color_diff < np.sqrt(radius**2 + r_err**2 + g_err**2))

    cut = np.logical_or(cut_1, cut_2)

    #mag_bins = np.arange(17., 24.1, 0.1)
    mag_bins = np.arange(17., mag_max+0.1, 0.1)
    mag_centers = 0.5 * (mag_bins[1:] + mag_bins[0:-1])
    magerr = np.tile(0., len(mag_centers))
    for ii in range(0, len(mag_bins) - 1):
        cut_mag_bin = (g > mag_bins[ii]) & (g < mag_bins[ii + 1])
        magerr[ii] = np.median(np.sqrt(radius**2 + r_err[cut_mag_bin]**2 + g_err[cut_mag_bin]**2))

    if return_all:
        return cut, mag_centers[f_isochrone(mag_centers) < 100], (f_isochrone(mag_centers) + magerr)[f_isochrone(mag_centers) < 100], (f_isochrone(mag_centers) - magerr)[f_isochrone(mag_centers) < 100]
    else:
        return cut


def isochrone_search(star_data, distance, age=12.0, Z=0.0002, graph=True, save=True):
    '''
    I'm assuming stars_data is a LSSTData or LSSTnEuclidData object
    '''
    #distance is given in kpc
    distance_modulus = coordinate_tools.distanceToDistanceModulus(distance)

    #the isochrone with Euclid, Roman, and LSST bands is 'mixed'
    iso = isochrone.Isochrone(
          age=age,
          metallicity=Z,
          distance_modulus=distance_modulus,
          survey= 'mixed',
          band_1= 'g',
          band_2= 'r')
    
    #cut = cut_isochrone_path(star_data.g_mag, star_data.r_mag,
    #                         star_data.g_magerr, star_data.r_magerr,
    #                         iso, radius = 0.1)
    
    cut = cut_isochrone_path(star_data.g_mag, star_data.r_mag,
                             star_data.g_magerr, star_data.r_magerr,
                             iso, radius = 0.1)
    isochrone_stars = star_data.apply_mask(cut)

    if graph==True:
        plotting_functions.isochrone_plot(iso, distance_modulus,
                                          star_data, isochrone_stars,
                                          save=save)
    
    return isochrone_stars
        

  

