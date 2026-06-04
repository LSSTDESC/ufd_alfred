import numpy as np
import matplotlib.pyplot as plt
from star_gal_sep import *
from utils import *
#-----------------------------------

def flux_ratio(data, band, c):
    psfFlux = data[f'{band}_psfFlux']
    psfFluxErr = data[f'{band}_psfFluxErr']
    cModelFlux = data[f'{band}_cModelFlux']
    cModelFluxErr = data[f'{band}_cModelFluxErr']

    flux_ratio = psfFlux / cModelFlux
    flux_ratio_err = np.sqrt((psfFluxErr / cModelFlux)**2
                         + ((psfFlux / cModelFlux**2)*cModelFluxErr)**2)
    return (1 - flux_ratio) + c*flux_ratio_err


def stellar_catalog(data, survey, func, band, lt_threshold, c=5/2):
    '''
    accepts which classifier to use, which threshold you want the stars to be less than
    returns selected stars'''
    if survey == 'dp1':
        #i'm sure there's a better way to do this
        if func == 'fluxratioerr':
            data['starClassifier'] = flux_ratio(data, band, c)
    if survey == 'dp2':
        if func == 'model_extendedness':
            data['starClassifier'] = data[f'{band}_model_extendedness']

    stars = data[data['starClassifier'] < lt_threshold]
    return stars
    
#diagnostic plots
def color_magnitude(star_df, band1, band2, colors, color_label, title, ax = None, save = False, plots_path = None):
    """
    Plots color-magnitude diagram, band1 vs band1-band2

    Parameters
    ----------
    star_df : Astropy Table
        Stellar catalog
    band1 : string
        Photometry band to be plotted on y-axis
    band2 : string
        Photometry band for the x-axis
        function wil plot {band}_psfFlux
    colors : Table column
        Colorbar data
    color_label : string
        Colorbar label
    title : string
        Title of plot/subplot
    ax (optional) : default None, else plt axes object 
        Axes to plot subplots on
        If none, will set subplot axes (1,1) 
    save (optional) : default False
        File names have form 'colormag_{title}.png'
    plots_path (optional) : default None, else string
        Where to save

    Returns
    -------
    Pretty plot
    """

    if ax == None:
        fig, axes = plt.subplots(1,1, figsize=(9,7))
        ax = axes

    band1_mag = flux2mag(star_df[f'{band1}_psfFlux'])
    band2_mag = flux2mag(star_df[f'{band2}_psfFlux'])

    _ = ax.scatter(band1_mag - band2_mag, band1_mag,
                c = colors,
                s = 10, cmap = 'jet')
    cbar = plt.colorbar(_)
    cbar.ax.invert_yaxis()
    cbar.set_label(color_label)
    ax.set(xlabel = f"{band1} - {band2}", ylabel = f"{band1}", xlim = (-1, 4),
              ylim = (30, 18), title = title)

    plt.tight_layout()
    if save == True:
        title = title.replace(' ', '').lower()
        plt.savefig(plots_path + f'/colormag_{title}.png')


def color_magnitude2(star_df, pltL_dict, pltR_dict, title = None, save = False, plots_path = None):
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

    fig, axes = plt.subplots(1,2, figsize = (18,7))
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
        title = title.replace(' ', '').lower()
        plt.savefig(plots_path + f'/colormag_{title}.png')


def color_color(star_df, band_x1, band_x2, band_y1, band_y2, colors, color_label, title, ax = None, save = False, plots_path = None):
    """
    Plots color-color, r-i vs g-r for different magnitudes

    Parameters
    ----------
    star_df : Astropy Table
        Stellar catalog
    band_x1, band_x2 : string
        Bands to subtract for the x-axis
    band_y1, band_y2 : string
        Bands to subtract for the y-axis
        Function plots band1_psfFlux - band2_psfFlux
    colors : Table column
        Data to be used for colorbar
    color_label : string
        Colorbar label
    title : string
        Plot title
    ax (optional) : default None, else plt axes object
        Axes to plot subplots on
        If none, will set subplot axes (1,1)
    save (optional) : default False
        File names have form 'colorcolor_{title}.png'
    plots_path : default None, else string
        Where to save

    Returns
    -------
    Pretty plot
    """

    if ax == None:
        fig, axes = plt.subplots(1,1, figsize=(9,7))
        ax = axes

    x1_mag = flux2mag(star_df[f"{band_x1}_psfFlux"])
    x2_mag = flux2mag(star_df[f"{band_x2}_psfFlux"])
    y1_mag = flux2mag(star_df[f"{band_y1}_psfFlux"])
    y2_mag = flux2mag(star_df[f"{band_y2}_psfFlux"])

    _ = ax.scatter(x1_mag - x2_mag, y1_mag - y2_mag,
                c = colors,
                s = 15, cmap = "jet_r")
    plt.colorbar(_, label = color_label)
    ax.set(ylabel = f'{band_y1} - {band_y2}', ylim = (-0.6,2),
              xlabel = f'{band_x1} - {band_x2}', xlim = (0.6,2.5),
              title = title)
    
    plt.tight_layout()
    if save == True:
        plt.savefig(plots_path + f'/colorcolor_{title.replace(' ', '').lower()}.png')
        print(plots_path + f'/colorcolor_{title.replace(' ', '').lower()}.png')


def star_gal_sep(df, survey_sep, colors, color_label, title, ax = None, y_bounds=(-0.2, 1.0), save = False, plots_path = ''):
    '''
    Plots the difference between PSF and cModel flux across magnitudes
    color-coded by star-galaxy classifiers
    Shows magnitude where things start getting confused
    
    Parameters
    ----------
    df : astropy table
        Table of all sources
    survey_sep : string
        Survey used for the separation, determines the y-axis
    colors : table column, potentially a 
        Data to be used for colorbar
    color_label : string
        Colorbar label
    title : string
        Plot title
    ax (optional) : default None, else plt axes object
        Axes to plot subplots on
        If none, will set subplot axes (1,1)
    y_bounds (optional) : tuple of floats
        Argument for y limits
        Found it necessary sometimes to zoom and enhance on LSST selector
    save (optional) : default False
        If True the file will be saved
    plots_path (optional) : default empty string
        Where to save

    Returns
    -------
    Pretty plot
    '''
    if ax == None:
        fig, axes = plt.subplots(1,1, figsize=(9,7))
        ax = axes

    i_mag = flux2mag(df['i_psfFlux'])
    if survey_sep == 'Euclid':
        y = df['MUMAX_MINUS_MAG']
        y_label = 'Euclid mu_max - mag'
    else:
        i_mag_cmodel = flux2mag(df["i_cModelFlux"])
        y = i_mag - i_mag_cmodel
        y_label = f'{survey_sep} i_psfFlux - i_cModelFlux'

    _ = ax.scatter(i_mag, y,
                   c = colors,
                   s = 10, cmap = "viridis_r")
    ax.set(title = title, xlabel = 'LSST i_psfFlux', ylabel = y_label, ylim = y_bounds)
    plt.colorbar(_,label = color_label)
    plt.tight_layout()
    
    if save == True:
        plt.savefig(plots_path + f'/stargalsep_{title.replace(' ', '').lower()}.png')


"""
def color_colorN(star_df, plt_dict_list, title = None, save = False, plots_path = None):
    


    
    N = len(plt_dict_list)
    fig, axes = plt.subplots(1,N, figsize = (9*N,7))
    axes = axes.flatten()

    for i, ax in enumerate(axes):
        plt_dict = plt_dict_list[i]
        color_color(star_df,
                    plt_dict['band_x1'], plt_dict['band_x2'],
                    plt_dict['band_y1'], plt_dict['band_y2'],
                    plt_dict['colors'], plt_dict['color_label'],
                    plt_dict['title'], ax = ax)

    if title is not None:
        fig.suptitle(title)

    plt.tight_layout()
    if save == True:
        plt.savefig(plots_path + f'/colormag_{title.replace(' ', '').lower()}.png')
"""
"""
def color_color_gradient(star_df, plt_dict, title = None, save = False, plots_path = None):
    ###
#plt_dict['band_x1'], plt_dict['band_x2'],
#                    plt_dict['band_y1'], plt_dict['band_y2'],
#                    plt_dict['colors'], plt_dict['color_label'],
#                    plt_dict['title']

    band_x1 = plt_dict['band_x1']
    band_x2 = plt_dict['band_x2']
    band_y1 = plt_dict['band_y1']
    band_y2 = plt_dict['band_y2']
### NEEDS TO BE FINISHED!!!!!! 
    bright = {'band_x1' : [], 'band_x2' : [], 'band_y1' : [], 'band_y2' : []
              'colors' : [], 'color_label' : plt_dict['color_label'], 'title' }
    middle = {'band_x1' : [], 'band_x2' : [], 'band_y1' : [], 'band_y2' : []
              'colors' : [],}
    dim = {'band_x1' : [], 'band_x2' : [], 'band_y1' : [], 'band_y2' : []
              'colors' : [],}
    for i in range(len(band_x1)):
        mag_i = band_x1[i]
        if (mag_i > 18) and (mag_i < 22):
            bright['band_x1'].append(mag_i)
            bright['band_x2'].append(band_x2[i])
            bright['band_y1'].append(band_y1[i])
            bright['band_y2'].append(band_y2[i])
            bright['colors'].append(plt_dict['colors'][i])
        elif (mag_i > 22) and (mag_i < 24):
            middle['band_x1'].append(mag_i)
            middle['band_x2'].append(band_x2[i])
            middle['band_y1'].append(band_y1[i])
            middle['band_y2'].append(band_y2[i])
            middle['colors'].append(plt_dict['colors'][i])
        else:
            dim['band_x1'].append(mag_i)
            dim['band_x2'].append(band_x2[i])
            dim['band_y1'].append(band_y1[i])
            dim['band_y2'].append(band_y2[i])
            dim['colors'].append(plt_dict['colors'][i])

    color_colorN(star_df, [bright, middle, dim], title = None, save = False, plots_path = None)
"""
