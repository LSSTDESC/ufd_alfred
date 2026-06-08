import sys
sys.path.append("../") #has to be run in same directory rn
import simple_adl.simple_adl.isochrone as isochrone
import simple_adl.simple_adl.coordinate_tools as coordinate_tools
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import scipy

def cut_isochrone_path(g, r, g_err, r_err, isochrone, radius=0.1, return_all=False):
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
    cut_2 = (color_diff < np.sqrt(0.1**2 + r_err**2 + g_err**2))

     # ...and now the other
    f_isochrone = scipy.interpolate.interp1d(mag_1_rgb, mag_1_rgb - mag_2_rgb, bounds_error=False, fill_value = 999.)
    color_diff = np.fabs((g - r) - f_isochrone(g))
    cut_1 = (color_diff < np.sqrt(0.1**2 + r_err**2 + g_err**2))

    cut = np.logical_or(cut_1, cut_2)
    
    mag_max = 24 #??? idk
    #mag_bins = np.arange(17., 24.1, 0.1)
    mag_bins = np.arange(17., mag_max+0.1, 0.1)
    mag_centers = 0.5 * (mag_bins[1:] + mag_bins[0:-1])
    magerr = np.tile(0., len(mag_centers))
    
    for ii in range(0, len(mag_bins) - 1):
        cut_mag_bin = (g > mag_bins[ii]) & (g < mag_bins[ii + 1])
        magerr[ii] = np.median(np.sqrt(0.1**2 + r_err[cut_mag_bin]**2 + g_err[cut_mag_bin]**2))
    if return_all:
        return cut, mag_centers[f_isochrone(mag_centers) < 100], (f_isochrone(mag_centers) + magerr)[f_isochrone(mag_centers) < 100], (f_isochrone(mag_centers) - magerr)[f_isochrone(mag_centers) < 100]
    else:
        return cut

def isochrone_search(star_data, distance, age=12.0, Z=0.0002, graph=True, plots_dir=''):
  #distance in kpc
  distance_modulus = coordinate_tools.distanceToDistanceModulus(distance)
  
  iso = isochrone.Isochrone(
          age=age,
          metallicity=Z,
          distance_modulus=distance_modulus,
          survey= 'lsst',
          band_1= 'g',
          band_2= 'r')
  star_data['g mag'] = flux2mag(star_data['g_psfFlux'])
  star_data['r mag'] = flux2mag(star_data['r_psfFlux'])
  star_data['g mag err'] = -2.5/np.log(10)*(star_data['g_psfFluxErr']/star_data['g_psfFlux'])
  star_data['r mag err'] = -2.5/np.log(10)*(star_data['r_psfFluxErr']/star_data['r_psfFlux'])
  star_data['g mag err'][~np.isfinite(star_data['g mag err'])] = np.nan
  star_data['r mag err'][~np.isfinite(star_data['r mag err'])] = np.nan
  
  cut = cut_isochrone_path(star_data['g mag'], star_data['r mag'], star_data['g mag err'], star_data['r mag err'], iso)
  star_data = star_data[cut]
 
  if graph==True: 
    fig, ax = plt.subplots(1,1, figsize=(6,6))
    index = np.min(np.where(iso.stage == iso.hb_stage)[0]) + 1
    ax.set(xlabel = 'g-r', ylabel = 'g', xlim = (-1,4), ylim = (28,18))
    ax.plot(iso.mag_1[0:index] - iso.mag_2[0:index], iso.mag_1[0:index] + distance_modulus, color='k')
    ax.plot(iso.mag_1[index:] - iso.mag_2[index:], iso.mag_1[index:] + distance_modulus, color = 'k')
    ax.scatter(data['g mag'] - data['r mag'], 
             data['g mag'], c='b')
  
    plt.savefig(plots_dir + '/iso_test.png')

  

