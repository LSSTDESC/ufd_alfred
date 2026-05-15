import numpy as np
import matplotlib.pyplot as plt
from star_gal_sep import *
from utils import *
#-----------------------------------


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
    dict_list = [pltA_dict, pltB_dict]
    
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
        plt.savefig(plots_path + f'colorcolor_{title.replace(' ', '').lower()}.png')

def color_colorN(star_df, plt_dict_list, title = None, save = False, plots_path = None):
    """

    """
    
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

def color_color_gradient(star_df, plt_dict, title = None, save = False, plots_path = None):
    """
plt_dict['band_x1'], plt_dict['band_x2'],
                    plt_dict['band_y1'], plt_dict['band_y2'],
                    plt_dict['colors'], plt_dict['color_label'],
                    plt_dict['title']
    """

    band_x1 = plt_dict['band_x1']
    band_x2 = plt_dict['band_x2']
    band_y1 = plt_dict['band_y1']
    band_y2 = plt_dict['band_y2']
### NEEDS TO BE FINISHED!!!!!! 
    bright = {'band_x1' : [], 'band_x2' : [], 'band_y1' : [], 'band_y2' : []
              'colors' : [],}
    middle = {"r mag" : [], "g mag" : [], "i mag" : [], "ext" : []}
    dim = {"r mag" : [], "g mag" : [], "i mag" : [], "ext" : []}

    for i in range(len(r_mag)):
        r_i = r_mag[i]
        if (r_i > 18) and (r_i < 22):
            bright["r mag"].append(r_i)
            bright["g mag"].append(g_mag[i])
            bright["i mag"].append(i_mag[i])
            bright["ext"].append(euc_colors.to_list()[i])
        elif (r_i > 22) and (r_i < 24):
            middle["r mag"].append(r_i)
            middle["g mag"].append(g_mag[i])
            middle["i mag"].append(i_mag[i])
            middle["ext"].append(euc_colors.to_list()[i])
        else:
            dim["r mag"].append(r_i)
            dim["g mag"].append(g_mag[i])
            dim["i mag"].append(i_mag[i])
            dim["ext"].append(euc_colors.to_list()[i])

    color_colorN(star_df, plt_dict_list, title = None, save = False, plots_path = None)


def star_gal_sep(merged_df, colors1, label1, colors2, label2, surveyname, y_bounds, save = False, add_lines = False, file_num = ''):
    '''
    Plots color-color, r-i vs g-r for different magnitudes
    Expected that star_df would be sorted from an LSST-based stellar classifier
    Point of this plot is to show how LSST's classification changes as sources get fainter

    Parameters
    ----------
    merged_df : pandas dataframe
        Dataframe of all sources
    colors1 : dataframe column
        Data to be used for subplot 1 colorbar
    label1 : string
        Subplot 1 colorbar label
        If label1 and label2 are to be different,
        label 1 should be the one for which you want the file named
    colors2 : dataframe column
        Data to be used for subplot 2 colorbar
    label2 : string
        Subplot 2 colorbar label
    surveyname : string
        Survey from which i band photometry is being pulled
        (where merged_df is getting its i_psfFlux and i_cModelFlux)
    y_bounds : tuple of floats
        Argument for subplot 1's ylim
        Found it necessary to zoom and enhance on LSST selector
    save (optional) : default False
        If True the file will be saved
    file_num (optional) : default ''
        You can optionally add a number if you don't want to overwrite
        the file previously saved with same name
        (file titles have form 'star-gal-sep_{label1}_selector_{file_num}')

    Returns
    -------
    Pretty plot
    '''

    fig, ax = plt.subplots(1,2, figsize=(18,6))

    i_mag = flux2mag(merged_df["i_psfFlux"].values)
    i_mag_cmodel = flux2mag(merged_df["i_cModelFlux"].values)

    _=ax[0].scatter(i_mag,
                i_mag - i_mag_cmodel,
                c = colors1,
                s = 10, cmap = "viridis_r")
    ax[0].set(xlabel = f'{surveyname} i_psfFlux', ylabel = f'{surveyname} i psf - cmodel mag', ylim = y_bounds)
    plt.colorbar(_,label = label1)
    _=ax[1].scatter(i_mag,
                merged_df['MUMAX_MINUS_MAG'],
                c = colors2,
                s = 10, cmap = 'viridis_r')
    plt.colorbar(_,label = label2)
    ax[1].set(xlabel = f'{surveyname} i_psfFlux', ylabel = "Euclid mu_max - mag")
    if add_lines == True:
        ax[1].plot(i_mag, (-0.07*i_mag)-1.2, label = "Cutoff Attempt")
    plt.tight_layout()
    if save == True:
        plt.savefig(my_plotspath + f'star-gal-sep_{label1}_selector_{file_num}.png')
