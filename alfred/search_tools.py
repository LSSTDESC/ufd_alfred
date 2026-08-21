import sys
#sys.path.append("../") #has to be run in same directory rn
import external.simple_adl.simple_adl.isochrone as isochrone
import external.simple_adl.simple_adl.coordinate_tools as coordinate_tools
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
        
#~~~~~~~~~ still learning this~~~~~~~~~
def get_maps():
    ''' for DP2 maps
    # define collection parameters
    REPO = 'dp2_prep'
    
    butler = Butler(REPO)
    registry = butler.registry
    
    # I don't think this dp2 repo was chained? So each stage was needed to be read at separate collections
    collection = ['LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage1',
                  'LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage2',
                  'LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage3',
                  'LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage4',]
    '''
    
    hsp_map = butler.get('deepCoadd_psf_size_consolidated_map_weighted_mean',
                         collections = collection,
                         band = 'i',
                         skymap = 'lsst_cells_v1',
                         tract = 5063,
                         parameters={"degrade_nside": 2048},)
    nside_frac=2048
    fracdet = hsp_map.fracdet_map(nside_frac).coverage_mask
    basis_1 = 'coord_ra'
    basis_2 = 'coord_dec'
    data = stars
    ra_select = np.median(data['coord_ra'])
    dec_select = np.median(data['coord_dec'])
    print('min and max ra:', data['coord_ra'].min(), data['coord_ra'].max(), data['coord_ra'].max()-data['coord_ra'].min())
    print('min and max dec:', data['coord_dec'].min(), data['coord_dec'].max(), data['coord_dec'].max()-data['coord_dec'].min())
      

def compute_char_density():
    ##COMPUTE CHARACTERISTIC DENSITY

    #only the stars that are reasonably bright
    cut_magnitude_threshold = (data['g mag'] < mag_max)
    #project into a cartesian image rather than spherical coords
    proj = projector.Projector(ra_select, dec_select)
    x, y = proj.sphereToImage(data[basis_1][cut_magnitude_threshold], data[basis_2][cut_magnitude_threshold]) # Trimmed magnitude range for hotspot finding
    #x_full, y_full = proj.sphereToImage(data[basis_1], data[basis_2]) # If we want to use full magnitude range for significance evaluation
    
    delta_x = 0.01
    area = delta_x**2
    smoothing = 2. / 60. # Was 3 arcmin
        #^ smoothing goes into sigma which is standard deviation for Gaussian kernel
        # but i don't see the impact of it on the graph
    bins = np.arange(-1., 1. + 1.e-10, delta_x)
        # why -8 to 8?? Am I missing something about how a region is defined
    centers = 0.5 * (bins[0: -1] + bins[1:]) #averaging the left and right boundaries of bins
    yy, xx = np.meshgrid(centers, centers)
    
    h = np.histogram2d(x, y, bins=[bins, bins])[0]
    
    h_g = scipy.ndimage.gaussian_filter(h, smoothing / delta_x)
    #don't quite understand the coordinates on these just yet
    plt.hist2d(x, y, bins=[bins, bins])
    plt.show()
    plt.imshow(h.T, extent=[bins[0], bins[-1], bins[0], bins[-1]], origin='lower')
    plt.show()
    plt.imshow(h_g.T, extent=[bins[0], bins[-1], bins[0], bins[-1]], origin='lower')
    plt.show()
    
    ## h_goodcoverage is supposed to be a mask where you restrict only to pixels with more than one observation
    #right now it only cuts out pixels with 0 observed stars (so just the empty pixels at the edge)
    delta_x_coverage = 0.1
    area_coverage = (delta_x_coverage)**2
    bins_coverage = np.arange(-5., 5. + 1.e-10, delta_x_coverage)
    h_coverage = np.histogram2d(x, y, bins=[bins_coverage, bins_coverage])[0]
    #h_goodcoverage = np.histogram2d(x[cut_goodcoverage], y[cut_goodcoverage], bins=[bins_coverage, bins_coverage])[0]
    h_goodcoverage = np.histogram2d(x, y, bins=[bins_coverage, bins_coverage])[0]
    
    n_goodcoverage = h_coverage[h_goodcoverage > 0].flatten()
    
    #characteristic_density = np.mean(n_goodcoverage) / area_coverage # per square degree
    characteristic_density = np.median(n_goodcoverage) / area_coverage # per square degree
    print('Characteristic density = {:0.1f} deg^-2'.format(characteristic_density))

    nside = 512 #need an nside < 2048 for splitting up subpix idk what it should be tho

    # Use pixels with fracdet ~1.0 to estimate the characteristic density
    if fracdet is not None:
        ## don't think this code is necessary if fracdet is the coveragemask
        #fracdet_zero = np.tile(0., len(fracdet)) #making new array same as fracdet but of 0s
        #cut = (fracdet != hp.UNSEEN)
        #fracdet_zero[cut] = fracdet[cut]
        
        print(nside_fracdet) #may be a problem that this i thought was supposed to be 2048
                            #but there was something about coarser resolution I didn't understand
        
        subpix_region_array = []
        for pix in np.unique(healpix.angToPix(nside, data[basis_1], data[basis_2])):
            subpix_region_array.append(healpix.subpixel(pix, nside, nside_fracdet))
        subpix_region_array = np.concatenate(subpix_region_array)
        print(subpix_region_array)
    
        # Compute mean fracdet in the region so that this is available as a correction factor
        cut = (fracdet[subpix_region_array] != hp.UNSEEN)
        mean_fracdet = np.mean(fracdet[subpix_region_array[cut]])
    
        # smau: this doesn't seem to be used in the non-local density estimation
        subpix_region_array = subpix_region_array[fracdet[subpix_region_array] > 0.99]
        subpix = healpix.angToPix(nside_fracdet, 
                                  data[basis_1][cut_magnitude_threshold], 
                                  data[basis_2][cut_magnitude_threshold]) # Remember to apply mag threshold to objects
        characteristic_density_fracdet = float(np.sum(np.in1d(subpix, subpix_region_array))) \
                                         / (hp.nside2pixarea(nside_fracdet, degrees=True) * len(subpix_region_array)) # deg^-2
        print('Characteristic density fracdet = {:0.1f} deg^-2'.format(characteristic_density_fracdet))
        
        # Correct the characteristic density by the mean fracdet value
        characteristic_density_raw = 1. * characteristic_density
        characteristic_density /= mean_fracdet 
        print('Characteristic density (fracdet corrected) = {:0.1f} deg^-2'.format(characteristic_density))