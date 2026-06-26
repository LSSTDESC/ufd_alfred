import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import seaborn as sns

def flux2mag(flux):
    zeropoint = 31.4
    index = flux.index if hasattr(flux, 'index') else None
    flux = np.asarray(flux, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        mag = -2.5 * np.log10(flux) + zeropoint
    mag[~np.isfinite(mag)] = np.nan
    if index is not None:
        return pd.Series(mag, index=index)
    return mag

def fluxerr2magerr(flux, flux_err):
    flux = np.ma.filled(np.ma.asarray(flux, dtype=float), fill_value=np.nan)
    flux_err = np.ma.filled(np.ma.asarray(flux_err, dtype=float), fill_value=np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        magerr = (2.5 / np.log(10)) * (flux_err / flux)
    magerr[~np.isfinite(magerr)] = np.nan
    return magerr

tract_dict={'47 Tuc' : [531, 532, 453, 454],
            'ECDFS' : [4848, 4849, 5063, 4636, 4637, 4638, 4847, 4850, 5061, 5062, 5064,
       5065, 5279, 5280, 5281, 5282], #Euclid's "Euclid Deep Field Fornax" encompasses Extended Chandra Deep Field South
            'EDFS' : [2078, 2079, 2080, 2232, 2233, 2234, 2235, 2236, 2237, 2392, 2393,
       2394, 2395, 2396, 2397, 2557, 2558, 2559, 2560, 2561, 2562, 2728,
       2729, 2730, 2731],
            'Fornax' : [4016, 4017, 4218, 4217],
            'FDSG' : [4016, 4017, 4218, 4217], #same as Fornax idk how it'll be referred to
            'Rubin_SV_095-25' : [5525, 5526],
            'Rubin_SV_38_7' : [10463, 10464, 10704],
            'LELF': [10464, 10221, 10222, 10704, 10705, 10463], # Low Ecliptic Latitude Field / Rubin_SV_38_7
            'Seagull' : [7850, 7849, 7610, 7611],
       }

def get_tracts(field):
    '''
    Input: field -- str, name of field in LSST data (case sensitive)
    Output: tracts -- list of ints, all the tracts that lie in the definition of that field
                            (depends on data release,
                            this function will always (try to) be the most up-to-date data release)
    '''
    #for DP1, from https://portal.nersc.gov/cfs/lsst/dp1/contributed-notebooks/DP1_Detector_Visits_NB1.html
    try:
        return tract_dict[field]
    except:
        raise Exception('Field not found. Check capitalization or spelling')

def get_field(tract):
    '''
    Input: tract -- int, single tract number
    Output: field -- str, field in which that tract is defined
    '''
    field_dict = {}
    for f in tract_dict.keys():
        for t in tract_dict[f]:
            field_dict[t] = f
    try:
        return field_dict[tract]
    except:
        raise Exception('Tract not found, check scope of data release')


## ~~~~~~~~~ PLOTS ~~~~~~~~~~~~

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
pad=20
size=1

def color_magnitude(df, band1, band2,
                    colors, color_label, selection_label,
                    title,
                    histogram = False, ax = None,
                    x_lim = (-1, 4), y_lim = (30, 18),
                    save = False, plots_path = None, filename = None):

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
    colors : data column OR str
        If plot is to be scatterplot or contour, then should be data for the mapping
        If plot is to be histogram, then should be str of the cmap to be used
    color_label : string
        Colorbar label
    title : string
        Title of plot/subplot
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
            if type(colors) == str:
                ax.hist2d(band1_mag - band2_mag, band1_mag, bins=200,
                          range = [[x_lim[0],x_lim[1]],[y_lim[1], y_lim[0]]],
                          cmin=1, cmap = colors, label = selection_label)
            else:
                _ = ax.scatter(band1_mag - band2_mag, band1_mag, s = size,
                               cmap = 'viridis_r',
                               c = colors, vmin=0, vmax=1,
                               label = selection_label)
                x = band1_mag - band2_mag
                y = band1_mag
                sns.kdeplot(x=x, y=y, fill=False, color="k")
        else:
            _ = ax.scatter(band1_mag - band2_mag, band1_mag, s = size,
                           cmap = 'viridis_r',
                           c = colors, vmin=0, vmax=1,
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
                colors, color_label,
                selection_label, title,
                histogram = False, ax = None,
                y_lim = None, x_lim = None,
                save = False, plots_path = None, filename = None):
    """
    Plots color-color for arbitrary bands

    Parameters
    ----------
    band_list : list of tuples
        1st entry of tuple is string labeling band
        2nd entries are the band data
    colors : data column OR str
        If plot is to be scatterplot or contour, then should be data for the mapping
        If plot is to be histogram, then should be str of the cmap to be used
    color_label : str
        Colorbar label
    selection_label : str
        Labels how the data was cut, to be displayed in legend
    title : str
        Plot title
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
            # 2D histograms require colormap which is str type
            if type(colors) == str:
                ax.hist2d(x1_mag - x2_mag, y1_mag - y2_mag, bins=200,
                          range = [[x_lim[0],x_lim[1]],[y_lim[0], y_lim[1]]],
                          cmin=1, cmap = colors, label = selection_label)
            # otherwise if color dimension is needed for data, it'll plot contours over scatter
            else:
                x=x1_mag - x2_mag
                y=y1_mag - y2_mag
                _ = ax.scatter(x, y, s = size,
                               c = colors, cmap = "viridis_r",
                               label = selection_label)
                sns.kdeplot(x=x, y=y, fill=False, color="k")
        else:
            _ = ax.scatter(x1_mag - x2_mag, y1_mag - y2_mag,
                           c = colors, vmin=0, vmax=1,
                           label = selection_label,
                           s = size, cmap = "viridis_r")

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
                 colors, color_label,
                 selection_label, title,
                 histogram = False, ax = None, line_plt = None,
                 colorbar_limits = (0,1), y_lim = (-0.1, 0.6), x_lim = (18,28),
                 save = False, plots_path = '', filename = None):
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
    colors : data column OR str
        If plot is to be scatterplot or contour, then should be data for the mapping
        If plot is to be histogram, then should be str of the cmap to be used
    color_label : str
        Colorbar label
        If None, won't display a colorbar
    selection_label : str
        Description of how the data was cut to be displayed in legend
        If _nolegend, won't display a legend
    title : str
        Plot title
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
            if type(colors) == str:
                ax.hist2d(i_mag, y, bins=200,
                          range = [[x_lim[0],x_lim[1]],[y_lim[0], y_lim[1]]],
                          cmin=1, cmap = colors, label = selection_label)
            # otherwise if color dimension is needed for data, it'll plot contours over scatter
            else:
                _ = ax.scatter(i_mag, y,
                               c = colors, vmin=colorbar_limits[0], vmax=colorbar_limits[1],
                               label = selection_label,
                               s = size, cmap = "viridis_r")
                sns.kdeplot(x=i_mag, y=y, fill=False, color="k")
        # normal scatter plot
        else:
            _ = ax.scatter(i_mag, y,
                           c = colors, vmin=colorbar_limits[0], vmax=colorbar_limits[1],
                           label = selection_label,
                           s = size, cmap = "viridis_r")

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
def color_magnitude2(star_df, pltL_dict, pltR_dict, title = None, save = False, plots_path = None, filename = None):
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