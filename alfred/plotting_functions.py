import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import seaborn as sns
from alfred import utils
import yaml
import os
## ~~~~~~~~~ PLOTS ~~~~~~~~~~~~

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
pad=20
size=1

with open('config.yaml', 'r') as ymlfile:
# this is a bit hard coded too but idk another work around
# if main.py stays same folder as the config then this should work
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    # assuming that it's cool that the whole github repo is considered "home"
    where = cfg['setup']['where']
    home_dir = os.path.expandvars(cfg['setup']['home_dir'][where])
    plots_dir = os.path.join(home_dir, cfg['output']['plots_dir'])
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)


def isochrone_plot(iso, distance_modulus, uncut_data, cut_data, save = True):
    fig, ax = plt.subplots(1,1, figsize=(6,6))
    index = np.min(np.where(iso.stage == iso.hb_stage)[0]) + 1
    ax.set(xlabel = 'g-r', ylabel = 'g', xlim = (-1,4), ylim = (28,18), title = f'Tract: {uncut_data.tract}, {uncut_data.lsst_survey} and {uncut_data.euclid_survey} Data')
    ax.plot(iso.mag_1[0:index] - iso.mag_2[0:index], iso.mag_1[0:index] + distance_modulus, color='k')
    ax.plot(iso.mag_1[index:] - iso.mag_2[index:], iso.mag_1[index:] + distance_modulus, color = 'k')

    #ax.scatter(uncut_data.g_mag - uncut_data.r_mag, uncut_data.g_mag, c='r', alpha = 0.3, label = 'Before cut')
    #ax.scatter(cut_data.g_mag - cut_data.r_mag, cut_data.g_mag, c='b', alpha = 0.5, label = 'After cut')
    ax.scatter(uncut_data.g_mag - uncut_data.r_mag, 
           uncut_data.g_mag,
           s=10, c = 'r', alpha =0.3,
           label = 'Before cut')
    ax.scatter(cut_data.g_mag - cut_data.r_mag, 
               cut_data.g_mag, 
               s=10, c = 'b', alpha =0.5, 
               label = 'After cut')
    ax.legend()

    plt.show()
    if save == True:
        plt.savefig(plots_dir + f'/isochrones/{uncut_data.tract}_{uncut_data.lsst_survey}_{uncut_data.euclid_survey}.png')


#to do: fix this
def match_validation_plots(merged_df, matches_lsst, matches_euclid, 
                           unmatched_lsst, unmatched_euclid, 
                           lsst_table, euclid_field, ds):
    ## Match Verification
        b = 60
        #1D histogram of matches and not matches
        fig, ax = plt.subplots(1,1, figsize=(13,5))
        match_vis_mag = flux2mag(matches_euclid['FLUX_VIS_2FWHM_APER']*(10**3))
        unmatch_vis_mag = flux2mag(unmatched_euclid['FLUX_VIS_2FWHM_APER']*(10**3))
        total_vis_mag = flux2mag(euclid_field['FLUX_VIS_2FWHM_APER']*(10**3))
        plt.hist(match_vis_mag, bins = b, histtype = 'step', color='b', label = 'Matched Euclid Sources')
        plt.hist(unmatch_vis_mag, bins = b, histtype = 'step', color='r', label = 'Unmatched Euclid Sources')
        plt.hist(total_vis_mag, bins = b, histtype = 'step', color='k', label = 'Total Euclid Sources')
        plt.xlabel('FLUX_VIS_2FWHM_APER mag')
        plt.xlim(16,36)
        plt.ylabel('Number counts')
        plt.yscale('log')
        plt.title(f'Tract {tract}: Euclid Source Match/Unmatch')
        plt.legend()
        plt.show()
        
        match_i_mag = flux2mag(matches_lsst['i_psfFlux'])
        unmatch_i_mag = flux2mag(unmatched_lsst['i_psfFlux'])
        total_i_mag = flux2mag(lsst_datafile['i_psfFlux'])
        fig, ax = plt.subplots(1,1, figsize=(13,5))
        plt.hist(match_i_mag, bins = b, histtype = 'step', color='c', label = 'Matched LSST Sources')
        plt.hist(unmatch_i_mag, bins = b, histtype = 'step', color='r', label = 'Unmatched LSST Sources')
        plt.hist(total_i_mag, bins = b, histtype = 'step', color='k', label = 'Total LSST Sources')
        plt.xlabel('i_psfFlux mag')
        plt.xlim(16,36)
        plt.ylabel('Number counts')
        plt.yscale('log')
        plt.title(f'Tract {tract}: LSST Source Match/Unmatch')
        plt.legend()
        plt.show()

        #2D histogram of matches in Euclid and LSST, does it look the same?
        fig, ax = plt.subplots(1,2, figsize=(13,5))
        _, _, _, im = ax[0].hist2d(matches_euclid['RIGHT_ASCENSION'], matches_euclid['DECLINATION'], bins=100)
        plt.colorbar(im, ax=ax[0])
        ax[0].set(title = f'Tract {tract}: Matches in Euclid', ylabel = "Dec (deg)", xlabel = "RA (deg)")
        ax[0].invert_xaxis()
        _, _, _, im = ax[1].hist2d(matches_lsst['coord_ra'], matches_lsst['coord_dec'], bins=100)
        plt.colorbar(im, ax=ax[1])
        ax[1].set(title = f'Tract {tract}: Matches in LSST', ylabel = "Dec (deg)", xlabel = "RA (deg)")
        ax[1].invert_xaxis()
        plt.show()

        '''
        #histogram of separation
        ds = ds * 3600 #ds is in degrees, want to plot in arcsecs
        #I forced in the function that matches would be <1"
        plt.hist(ds, histtype='step', range=(0,1))
        plt.xlabel('separation [arcsec]')
        plt.title(f'Tract {tract}: Matched Source Separation')
        plt.tight_layout()
        plt.show()

        #checking that dec and DECLINATION relation is slope of 1
        plt.scatter(merged_df['coord_dec'],merged_df['DECLINATION'],)
        '''
        i_mag = flux2mag(lsst_datafile['i_psfFlux'])
        vis_mag = flux2mag(euclid_field['FLUX_VIS_PSF']*10**3)
        
        print(len(i_mag))
        
        lsst_datafile1 = lsst_datafile #[(i_mag < 22)]
        euclid_field1 = euclid_field #[(vis_mag < 22)]
        
        plt.scatter(euclid_field1['RIGHT_ASCENSION'], euclid_field1['DECLINATION'], 
                    marker = '+', label = 'Euclid Sources', #c = flux2mag(euclid_field1['FLUX_VIS_PSF']*10**3), 
                   )
        plt.scatter(lsst_datafile1['coord_ra'], lsst_datafile1['coord_dec'], 
                    marker = 'x', label = 'LSST Sources', #c = flux2mag(lsst_datafile1['i_psfFlux']),
                   )
        plt.xlim(59.85, 59.84)
        plt.xlabel('RA (deg)')
        plt.ylim(-48.55,-48.54)
        plt.ylabel('DEC (deg)')
        #plt.colorbar()
        plt.legend()
        plt.title('Euclid and LSST before matching, \n restricted i_mag & vis_mag < 22')
        plt.show()
        
        ra_diff = (merged_df['RIGHT_ASCENSION'] - merged_df['coord_ra'])*3600
        dec_diff = (merged_df['DECLINATION'] - merged_df['coord_dec'])*3600
        plt.hist(ra_diff, bins = 100)
        plt.title('Merged Catalog RA difference')
        plt.xlabel('Euclid RA - LSST RA (arcsec)')
        plt.xlim(-1,1)
        plt.show()
        
        plt.hist(dec_diff, bins = 100)
        plt.xlim(-0.5,0.5)
        plt.xlabel('Euclid DEC - LSST DEC (arcsec)')
        plt.title('Merged Catalog DEC difference')
        plt.show()
    #save to merged_verification

def color_magnitude(df, band1, band2,
                    color_data, color_label,
                    selection_label, title,
                    colors = 'viridis',
                    histogram = False, ax = None,
                    colorbar_limits = (0,1),
                    x_lim = (-1, 4), y_lim = (30, 18),
                    save = False, filename = None):

    # NOTE: not generalized yet, ONLY plots LSST bands
    """
    Plots color-magnitude diagram, band1 vs band1-band2

    Parameters
    ----------
    df : Astropy table, Pandas dataframe, or structured Numpy array
        Data structure of objects to plot
    band1 : string
        LSST photometry band to be plotted on y-axis
        function wil plot {band}_psfFlux mag
    band2 : string
        LSST photometry band for the x-axis
        function wil plot {band}_psfFlux mag
    color_data : data column OR None
        If plot is to be scatterplot or contour, then should be data for the mapping
        If plot is to be histogram, then None
    color_label : string
        Colorbar label
    title : string
        Title of plot/subplot
    colors (optional) : string, default 'viridis'
        Color scheme to use
    histogram (optional) : boolean, default False
        If True, will plot either histogram or contour (depending on colors type)
        If False, will plot scatterplot
    ax (optional) : default None, else plt axes object
        Axes to plot subplots on
        If none, will set subplot axes (1,1)
    x_lim (optional) : tuple, default (-0.1, 0.6)
        Limits of x-axis
    y_lim (optional) : tuple, default (30, 18)
        Limits of y-axis
        Note: should be backwards because of weirdness of magnitude system
    save (optional) : default False
        File names have form 'colormag_{title}.png'
    plots_path (optional) : default None, else string
        Where to save
    filename (optional) : default None, else str
        Name of file to save as

    Returns
    -------
    Pretty plot
    """

    if ax == None:
        fig, axes = plt.subplots(1,1, figsize=(7,5))
        ax = axes

    band1_mag = flux2mag(df[f'{band1}_psfFlux'])
    band2_mag = flux2mag(df[f'{band2}_psfFlux'])
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='.*colormapping.*')
        warnings.filterwarnings("ignore", message='.*labels.*')
        if histogram == True:
            if color_data is None:
                ax.hist2d(band1_mag - band2_mag, band1_mag, bins=200,
                          range = [[x_lim[0],x_lim[1]],[y_lim[1], y_lim[0]]],
                          cmin=1, cmap = colors, label = selection_label, )
            else:
                _ = ax.scatter(band1_mag - band2_mag, band1_mag, s = size,
                               cmap = colors, c = color_data,
                               vmin=colorbar_limits[0], vmax=colorbar_limits[1],
                               label = selection_label)
                x = band1_mag - band2_mag
                y = band1_mag
                sns.kdeplot(x=x, y=y, fill=False, color="k")
        else:
            _ = ax.scatter(band1_mag - band2_mag, band1_mag, s = size,
                           cmap = colors, c = color_data,
                           vmin=colorbar_limits[0], vmax=colorbar_limits[1],
                           label = selection_label)

        ax.set_title(title, pad=pad)
        ax.set(xlabel = f"{band1} - {band2}", ylabel = f"{band1}", xlim = x_lim, ylim = y_lim)
        # cleaning up, not displaying everything if not needed
        if selection_label[0] != '_':
            ax.legend()
        if color_label is not None:
            plt.colorbar(_,label=color_label)
        plt.tight_layout()

        if save == True:
            if filename is None:
                title = title.replace(' ', '').replace('-', '_').lower()
                filename = title
            plt.savefig(plots_path + f'/colormag/{filename}.png')


def color_color(band_list,
                color_data, color_label,
                selection_label, title,
                colors = 'viridis',
                histogram = False, ax = None,
                y_lim = None, x_lim = None,
                save = False, filename = None):
    """
    Plots color-color for arbitrary bands

    Parameters
    ----------
    band_list : list of tuples
        1st entry of tuple is string labeling band
        2nd entries are the band data
    colors : data column OR None
        If plot is to be scatterplot or contour, then should be data for the mapping
        If plot is to be histogram, then None
    color_label : str
        Colorbar label
    selection_label : str
        Labels how the data was cut, to be displayed in legend
    title : str
        Plot title
    colors (optional) : str, default 'viridis'
        Color scheme to use
    histogram (optional) : boolean, default False
        If True, will plot either histogram or contour (depending on colors type)
        If False, will plot scatterplot
    ax (optional) : default None, else plt axes object
        Axes to plot subplots on
        If none, will set subplot axes (1,1)
    y_lim (optional) : default None, else tuple
        Limits of y-axis
        Found it necessary sometimes to zoom and enhance on LSST selector
    x_lim (optional) : default None, else tuple
        Limits of x-axis
    save (optional) : default False
        File names have form 'colorcolor_{title}.png'
    plots_path : default None, else string
        Where to save
    filename (optional) : default None, else str
        Name of file to save as

    Returns
    -------
    Pretty plot
    """

    if ax == None:
        fig, axes = plt.subplots(1,1, figsize=(7,6))
        ax = axes

    band_x1, band_x2, band_y1, band_y2 = band_list
    x1_mag = band_x1[1]
    x2_mag = band_x2[1]
    y1_mag = band_y1[1]
    y2_mag = band_y2[1]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='.*colormapping.*')
        warnings.filterwarnings("ignore", message='.*labels.*')
        if histogram == True:
            if color_data is None:
            # 2D histograms require colormap which is str type
                ax.hist2d(x1_mag - x2_mag, y1_mag - y2_mag, bins=200,
                          range = [[x_lim[0],x_lim[1]],[y_lim[0], y_lim[1]]],
                          cmin=1, cmap = colors, label = selection_label)
            # otherwise if color dimension is needed for data, it'll plot contours over scatter
            else:
                x=x1_mag - x2_mag
                y=y1_mag - y2_mag
                _ = ax.scatter(x, y, s = size,
                               c = color_data, cmap = colors,
                               label = selection_label)
                sns.kdeplot(x=x, y=y, fill=False, color="k")
        else:
            _ = ax.scatter(x1_mag - x2_mag, y1_mag - y2_mag,
                           c = color_data, vmin=0, vmax=1,
                           label = selection_label,
                           s = size, cmap = colors)

        ax.set_title(title, pad=pad)
        ax.set(ylabel = f'{band_y1[0]} - {band_y2[0]}',
               xlabel = f'{band_x1[0]} - {band_x2[0]}')
        # cleaning up the place, not displayin/changing everything if not needed
        if selection_label[0] != '_':
            ax.legend()
        if color_label is not None:
            plt.colorbar(_,label=color_label)
        if y_lim is not None:
            ax.set(ylim=y_lim)
        if x_lim is not None:
            ax.set(xlim=x_lim)
        plt.tight_layout()

        if save == True:
            if filename is None:
                title = title.replace(' ', '').replace('-', '_').lower()
                filename = title
            plt.savefig(plots_path + f'/colorcolor/{filename}.png')

#----

def star_gal_sep(df, separator,
                 color_data, color_label,
                 selection_label, title, colors = 'viridis_r',
                 histogram = False, ax = None, line_plt = None,
                 colorbar_limits = (0,1), y_lim = (-0.1, 0.6), x_lim = (18,28),
                 save = False, filename = None):
    '''
    Plots the morphology across magnitudes to show star-galaxy classifiers
    Shows magnitude where things start getting confused

    Parameters
    ----------
    df : Astropy table, Pandas dataframe, or structured Numpy array
        Data structure of objects to plot
    separator : str
        Separation parameter to plot on the y-axis
        Generally assumed to be a Euclid column title
        LSST i psf - cmodel or i psf / cmodel are specially-handled cases
    colors : data column OR None
        If plot is to be scatterplot or contour, then should be data for the mapping
        If plot is to be histogram, then None
    color_label : str
        Colorbar label
        If None, won't display a colorbar
    selection_label : str
        Description of how the data was cut to be displayed in legend
        If _nolegend, won't display a legend
    title : str
        Plot title
    colors (optional) : str, default 'viridis_r'
        Color scheme to use
    histogram (optional) : boolean, default False
        If True, will plot either histogram or contour (depending on colors type)
        If False, will plot scatterplot
    ax (optional) : default None, else plt axes object
        Axes to plot subplots on
        If none, will set subplot axes (1,1)
    colorbar_limits (optional) : tuple, default (0, 1)
        Limits for the colorbar (changes vmin, vmax)
    y_lim (optional) : tuple, default (-0.1, 0.6)
        Limits of y-axis
        Found it necessary sometimes to zoom and enhance on LSST selector
    x_lim (optional) : tuple, default (18,28)
        Limits of x-axis
    save (optional) : default False
        If True the file will be saved
    plots_path (optional) : str, default empty
        Where to save
    filename (optional) : default None, else str
        Name of file to save as

    Returns
    -------
    Pretty plot
    '''
    if ax == None:
        fig, axes = plt.subplots(1,1, figsize=(7,5))
        ax = axes

    # x-axis data
    i_mag = flux2mag(df['i_psfFlux'])
    # y-axis data depends on which morphology you want to demonstrate
    #  rubin morphologies, special because they requires operations:
    if separator == 'i psf - cmodel':
        i_mag_cmodel = flux2mag(df["i_cModelFlux"])
        y = i_mag - i_mag_cmodel
        y_label = f'LSST i_psfFlux - i_cModelFlux'
    elif separator == 'i psf / cmodel':
        i_mag_cmodel = flux2mag(df["i_cModelFlux"])
        y = i_mag / i_mag_cmodel
        y_label = f'LSST i_psfFlux / i_cModelFlux'
    #  otherwise I'm assuming it's just using a Euclid separator
    else:
        y = df[separator]
        y_label = 'Euclid ' + separator

    # the code will start yelling about the colorbar and labels depending on how you plot it
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='.*colormapping.*')
        warnings.filterwarnings("ignore", message='.*labels.*')
        # if you want to display by histogram/contour plot bc scatter gets oversaturated
        if histogram == True:
            # 2D histograms require colormap which is str type
            if color_data is None:
                ax.hist2d(i_mag, y, bins=200,
                          range = [[x_lim[0],x_lim[1]],[y_lim[0], y_lim[1]]],
                          cmin=1, cmap = colors, label = selection_label)
            # otherwise if color dimension is needed for data, it'll plot contours over scatter
            else:
                _ = ax.scatter(i_mag, y,
                               c = color_data, vmin=colorbar_limits[0], vmax=colorbar_limits[1],
                               label = selection_label,
                               s = size, cmap = colors)
                sns.kdeplot(x=i_mag, y=y, fill=False, color="k")
        # normal scatter plot
        else:
            _ = ax.scatter(i_mag, y,
                           c = color_data, vmin=colorbar_limits[0], vmax=colorbar_limits[1],
                           label = selection_label,
                           s = size, cmap = colors)

        ax.set(xlabel = 'LSST i_psfFlux mag', ylabel = y_label, ylim = y_lim, xlim = x_lim)
        ax.set_title(title, pad=pad) #title is separate so I can have the pad
        # cleaning up the place, not displaying everything if not needed
        if selection_label[0] != '_':
            ax.legend()
        if line_plt is not None:
            ax.plot(line_plt[0],line_plt[1],line_plt[2])
        if color_label is not None:
            plt.colorbar(_,label=color_label)
        plt.tight_layout()

        if save == True:
            if filename is None:
                title.replace(' ', '').replace('-','_').lower()
                filename = title
            plt.savefig(plots_path + f'/stargalsep/{filename}.png')




## ~~~~~ probably very outdated function ~~~~
def color_magnitude2(star_df, pltL_dict, pltR_dict, title = None, save = False, filename = None):
    """
    Plots 2 subplots of color-magnitude diagrams, band1 vs band1-band2

    Parameters
    ----------
    star_df : Astropy Table
        Stellar catalog
    pltL/R_dict : dictionary
        Params for left/right subplot
        {band1 : , band2 : , colors : , color_label : , title : }
        See color_magnitude for definitions of above
    title (optional) : default None, else string,
        Super title for plot
    save (optional) : default False
        File names have form 'colormag_{title}.png'
    plots_path (optional) : default None, else string
        Where to save

    Returns
    -------
    Pretty plot
    """

    fig, axes = plt.subplots(1,2, figsize = (14,5))
    axes = axes.flatten()
    dict_list = [pltL_dict, pltR_dict]

    for i, ax in enumerate(axes):
        color_magnitude(star_df,
                        dict_list[i]['band1'], dict_list[i]['band2'],
                        dict_list[i]['colors'], dict_list[i]['color_label'], dict_list[i]['title'],
                        ax = ax)

    if title is not None:
        fig.suptitle(title)
    plt.tight_layout()
    if save == True:
        if filename is None:
            title = title.replace(' ', '').replace('-', '_').lower()
            filename = title
        plt.savefig(plots_path + f'/colormag/{filename}.png')
