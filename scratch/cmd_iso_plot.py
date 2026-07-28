#import ugali.isochrone
import simple_adl.simple_adl.isochrone as isochrone
import simple_adl.simple_adl.coordinate_tools as coordinate_tools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from useful_functions_global import *

my_path =  '/sdf/data/rubin/user/kexcel/'
my_plotspath = my_path + 'plots/'

merged_df = pd.read_csv('/sdf/data/rubin/user/kexcel/dp1_euclid_merged.csv')

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
'''
iso = ugali.isochrone.factory(
    survey='des',
    name='Bressan2012',
    age=12,  # Gyr
    metallicity=0.00020, # Z
    distance_modulus=distance_modulus,
    band_1 = 'g',
    band_2 = 'r'
)
'''
'''
def cut_isochrone_path(g, r, g_err, r_err, isochrone, mag_max, radius=0.1, return_all=False):
    """
    Cut to identify objects within isochrone cookie-cutter.
    """
    if np.all(isochrone.stage == 'Main'):
        # Dotter case
        index_transition = len(isochrone.stage)
    else:
        # Other cases
        index_transition = np.nonzero(isochrone.stage >= isochrone.hb_stage)[0][0] # + 1 

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

iso_selection = cut_isochrone_path(region.data[survey.mag_dered_1], region.data[survey.mag_dered_2], 
                                   region.data[survey.mag_err_1], region.data[survey.mag_err_2], iso, survey.catalog['mag_max'],
                                   radius=0.1, return_all=True)
'''

def color_magnitude(star_df1, plt_title1, plt1colors, file_name, n = 1, cbar1label = None, scnd_plot = [None, None, None, None], isochrone_dict = False, save = False):
    '''
    Plots color-magnitude diagram, g vs g-r
    Expected that star_df would be sorted from a certain stellar classifier
        (i.e. just the sources LSST would classify as stars), but then those
        sources are plotted with both Euclid and LSST classifier colors

    Parameters
    ----------
    star_df : pandas dataframe
        Dataframe of the sources classified as stars
    euc_colors : dataframe column
        Data to be used for subplot 1 colorbar
        Euclid stellar classifier data (e.g. 'POINT_LIKE_PROB')
    euc_label : string
        Subplot 1 colorbar label
        Name of Euclid star classifier
    lsst_colors : dataframe column
        Data to be used for subplot 2 colorbar
        LSST-based stellar classifier data (e.g. 'i_SizeExtendedness')
    lsst_label : string
        Subplot 2 colorbar label
        Name of LSST star classifier
    starselector_name : string
        The survey with which the star_df has been sorted
        (e.g. if sorting based on 'i_SizeExtendedness' from DP1, this argument should be DP1)
    isochrone (optional) : default False
        if 
    save (optional) : default False
        If True the file will be saved
    file_num (optional) : default ''
        You can optionally add a number if you don't want to overwrite
        the file previously saved with same name
        (file titles have form 'colormag_{starselector_name}stars_{lsst_label}_selector_{file_num}')

    Returns
    -------
    Pretty plot
    '''

    fig, ax = plt.subplots(1,n, figsize=(15,6))

    if n==1:
        ax1 = ax
    else:
        ax1 = ax[0]
        star_df2, plt_title2, plt2colors, cbar2label = scnd_plot

    g_mag = flux2mag(star_df1["g_psfFlux"].values)
    r_mag = flux2mag(star_df1["r_psfFlux"].values)

    _ = ax1.scatter(g_mag - r_mag, g_mag,
                c=plt1colors,
                s=2)
    if cbar1label is not None:
        cbar1 = plt.colorbar(_)#, ticks=[0.48, 0.4, 0.3, 0.2, 0.1, 0.02])
        cbar1.ax.invert_yaxis()
        #cbar1.ax.set_yticklabels(['more \n extended', '0.4', '0.3', '0.2', '0.1', 'less \n extended'])
        cbar1.set_label(cbar1label)
    ax1.set(xlabel="g - r (LSST)", ylabel="g (LSST)", 
            xlim = (-1, 4), ylim = (28, 18), 
            title=plt_title1)

    if n>1:
        g_mag = flux2mag(star_df2["g_psfFlux"].values)
        r_mag = flux2mag(star_df2["r_psfFlux"].values)
        _ = ax[1].scatter(g_mag - r_mag, g_mag,
                    c=plt2colors,
                    s= 2)
        ax[1].set(xlabel="g - r (LSST)", ylabel="g (LSST)", 
                  xlim = (-1, 4), ylim = (28, 18), 
                  title=plt_title2)
        if cbar2label is not None:
            cbar2 = plt.colorbar(_) #, ticks=[0.02, 0.2, 0.4, 0.6, 0.8, 0.98])
            #cbar2.ax.set_yticklabels(['more \n extended', '0.2', '0.4', '0.6', '0.8', 'less \n extended'])
            cbar2.set_label(cbar2label)
            
    if isochrone_dict != False:
        iso, distance_dict = isochrone_dict
        for distance in distance_dict.keys():
            if distance is None:
                distance_modulus = iso.distance_modulus
            else:
                distance_modulus = coordinate_tools.distanceToDistanceModulus(distance)
            index = np.min(np.where(iso.stage == iso.hb_stage)[0]) + 1
            label = 'Distance = %i kpc'%(distance)
            linestyle = distance_dict[distance]['ls']
            color = distance_dict[distance]['color']
            ax1.plot(iso.mag_1[0:index] - iso.mag_2[0:index], iso.mag_1[0:index] + distance_modulus, 
                     label=label, ls = linestyle, c = color)
            ax1.plot(iso.mag_1[index:] - iso.mag_2[index:], iso.mag_1[index:] + distance_modulus, 
                     ls = linestyle, c = color)
            ax1.legend(title = f'{survey} {bands[0]}, {bands[1]} bands')
            if n != 1:
                ax[1].plot(iso.mag_1[0:index] - iso.mag_2[0:index], iso.mag_1[0:index] + distance_modulus, 
                           label=label, ls = linestyle, c = color)
                ax[1].plot(iso.mag_1[index:] - iso.mag_2[index:], iso.mag_1[index:] + distance_modulus, 
                           ls = linestyle, c = color)
                ax[1].legend(title = f'{survey} {bands[0]}, {bands[1]} bands')
            
    plt.tight_layout()
    if save == True:
        plt.savefig(my_plotspath + file_name)



euclid_stars = merged_df[merged_df['POINT_LIKE_PROB']>0.5]
lsst_stars = merged_df[merged_df['i_SizeExtendedness']<0.5]

color_magnitude(lsst_stars, "LSST Only", 'C9', "cmd_presentation.png", 
                n = 2, cbar1label = None, scnd_plot = [euclid_stars, "LSST + Euclid", 'C9', None], 
                isochrone_dict = (iso, {300 : {'color' : 'black', 'ls' : '-'}, 2000 : {'color' : 'black', 'ls' : '--'}}), 
                save = True)

#fig, ax = plt.subplots(1,1, figsize=(7,5))
#simple_plt.plot_cmd_sep(flux2mag(lsst_stars["g_psfFlux"].values), flux2mag(lsst_stars["r_psfFlux"].values), ax, sdata=None, cbar=True, show_iso=True, iso_selection=iso, fname="cmd_presentation2.png", save=True)