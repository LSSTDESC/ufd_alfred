import numpy as np
import matplotlib.pyplot as plt
from star_gal_sep import *
from utils import *
#-----------------------------------


def color_magnitude(star_df, band1, band2, colors, color_label, title, ax = None, save = False, plots_path = None):
    '''
    Plots color-magnitude diagram, band1 vs band1-band2

    Parameters
    ----------
    star_df : astropy table, stellar catalog
    band1 : string, band to be plotted on y-axis
    band2 : string, band for the x-axis 
        (function plots {band}_psfFlux)
    colors : table column, colorbar data
    color_label : string, colorbar label
    title : string, title of plot
    ax (optional) : default None (if none will set subplot axes (1,1)) 
    save (optional) : default False
        (file titles have form 'colormag_{title}.png')
    plots_path (optional) : default None, else string, where to save

    Returns
    -------
    Pretty plot
    '''
    if ax = None:
        fig, axes = plt.subplots(1,1, figsize=(9,7))
        ax = axes[0]

    band1_mag = flux2mag(star_df[f"{band1}_psfFlux"].values)
    band2_mag = flux2mag(star_df[f"{band2}_psfFlux"].values)

    _ = ax.scatter(band1_mag - band2_mag, band1_mag,
                c = colors,
                s = 10, cmap = "jet")
    cbar = plt.colorbar(_)#, ticks=[0.48, 0.4, 0.3, 0.2, 0.1, 0.02])
    cbar.ax.invert_yaxis()
    #cbar1.ax.set_yticklabels(['more \n extended', '0.4', '0.3', '0.2', '0.1', 'less \n extended'])
    cbar.set_label(color_label)
    ax.set(xlabel = f"{band1} - {band2}", ylabel = f"{band1}", xlim = (-1, 4),
              ylim = (30, 18), title = title)

    plt.tight_layout()
    if save == True:
        plt.savefig(plots_path + f'colormag_{title}.png')

def color_magnitude2(star_df, pltA_dict, pltB_dict, title = None, save = False, plots_path = None)
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
        plt.savefig(plots_path + f'colormag_{title}.png')


def color_color(star_df, euc_colors, euc_label, lsst_label, starselector_name, save = False, file_num = ''):
    '''
    Plots color-color, r-i vs g-r for different magnitudes
    Expected that star_df would be sorted from an LSST-based stellar classifier
        Point of this plot is to show how LSST's classification changes as sources get fainter

    Parameters
    ----------
    star_df : pandas dataframe
        Dataframe of the sources LSST classifies as stars
    euc_colors : dataframe column
        Data to be used for colorbar
        Euclid stellar classifier data (e.g. 'POINT_LIKE_PROB')
    euc_label : string
        Colorbar label
        Name of Euclid star classifier
    lsst_label : string
        Name of LSST star classifier used for star_df
        Only affects file name
    starselector_name : string
        The survey with which the star_df has been sorted
        (for now will just be DP1, just allowing for flexibility as more LSST data released)
    save (optional) : default False
        If True the file will be saved
    file_num (optional) : default ''
        You can optionally add a number if you don't want to overwrite
        the file previously saved with same name
        (file titles have form 'colorcolor_{starselector_name}stars_{lsst_label}_selector_{file_num}')

    Returns
    -------
    Pretty plot
    '''

    fig, ax = plt.subplots(1,3, figsize=(20,6))
    plt.suptitle(f'{starselector_name} \'Stars\' Photometry', fontsize = 18)

    g_mag = flux2mag(star_df["g_psfFlux"].values)
    r_mag = flux2mag(star_df["r_psfFlux"].values)
    i_mag = flux2mag(star_df["i_psfFlux"].values)

    bright = {"r mag" : [], "g mag" : [], "i mag" : [], "ext" : []}
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
    for dic in [bright, middle, dim]:
        for key in dic.keys():
            dic[key] = np.array(dic[key])

    _ = ax[0].scatter(bright["g mag"] - bright["r mag"], bright["r mag"] - bright["i mag"],
                c=bright["ext"],
                s=15, cmap="jet_r")
    _ = ax[1].scatter(middle["g mag"] - middle["r mag"], middle["r mag"] - middle["i mag"],
                c=middle["ext"],
                s=15, cmap="jet_r")
    _ = ax[2].scatter(dim["g mag"] - dim["r mag"], dim["r mag"] - dim["i mag"],
                c=dim["ext"],
                s=15, cmap="jet_r")
    plt.colorbar(_,label=euc_label)
    ax[0].set(ylabel='r - i', ylim = (-0.6,2),
              xlabel='g - r', xlim = (-0.6,2.5),
              title='r mag 20 - 22')
    ax[1].set(ylabel='r - i', ylim = (-0.6,2),
              xlabel='g - r', xlim = (-0.6,2.5),
              title='r mag 22 - 24')
    ax[2].set(ylabel='r - i', ylim = (-0.6,2),
              xlabel='g - r', xlim = (-0.6,2.5),
              title='r mag 24 - 26')
    plt.tight_layout()
    if save == True:
        plt.savefig(my_plotspath + f"colorcolor_{starselector_name}stars_{lsst_label}_selector_{file_num}.png")


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
