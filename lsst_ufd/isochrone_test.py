import sys
sys.path.append("../")
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

'''
distance = 300 # kpc
#distance_modulus = ugali.utils.projector.distanceToDistanceModulus(distance)
distance_modulus = coordinate_tools.distanceToDistanceModulus(distance)

survey = 'LSST'
bands = ['g', 'r']

iso = isochrone.Isochrone(
        age=12.0,
        metallicity=0.0002,
        distance_modulus=distance_modulus,
        survey= survey.lower(),
        band_1= bands[0],
        band_2= bands[1])

g = np.linspace(28,18,50)
g_r = np.linspace(-1,4,50)
G, G_R = np.meshgrid(g,g_r)
R = -1*G_R + G
data = Table()
data['g'] = G
data['r'] = R
data['g_err'] = np.linspace(0,0,50)
data['r_err'] = np.linspace(0,0,50)

cut = cut_isochrone_path(data['g'], data['r'], data['g_err'], data['r_err'], iso, radius=0.1, return_all=False)
data = data[cut]

fig, ax = plt.subplots(1,1, figsize=(6,6))
index = np.min(np.where(iso.stage == iso.hb_stage)[0]) + 1
ax.set(xlabel = 'g-r', ylabel = 'g', xlim = (-1,4), ylim = (28,18))
ax.plot(iso.mag_1[0:index] - iso.mag_2[0:index], iso.mag_1[0:index] + distance_modulus)
ax.plot(iso.mag_1[index:] - iso.mag_2[index:], iso.mag_1[index:] + distance_modulus)
ax.plot(data['g']-data['r'], data['g'], marker = 'o')

plt.savefig('/global/u2/k/kexcell/ultrafaints/plots/iso_test.png')
'''
