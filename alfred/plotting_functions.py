import numpy as np
import matplotlib.pyplot as plt
import warnings
import seaborn as sns
from alfred import utils
import yaml
import os
from matplotlib.patches import Circle, Annulus
from pyvo.dal import sia
from astropy.utils.data import download_file
from astropy.wcs import WCS
## ~~~~~~~~~ PLOTS ~~~~~~~~~~~~

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
pad=20
size=1

try:
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
except:
    print('No config file, this will mess up saving')

#~~~~~~~~~START CANDIDATE SCATTERPLOT/IMAGE FUNCTIONS~~~~~~~~~~~~~~~~~~~
width=0.0002
def candidate_scatterplot(Peak, ax = None, legend=True):
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    stars = Peak.member_candidates
    ax.scatter(stars.ra,stars.dec,
               s=50, facecolors='none', edgecolors='C0',
               label='Candidate member stars')
    ax.scatter(Peak.ra, Peak.dec,
               marker='X',c='C1',
               label=f'Peak center and radius outline')
    radius = Annulus((Peak.ra, Peak.dec),Peak.r, width, color='C1')
    ax.add_artist(radius)
    ax.set(xlabel='RA (deg)', ylabel='Dec (deg)',
           xlim=(Peak.ra-Peak.r*1.2, Peak.ra+Peak.r*1.2),ylim=(Peak.dec-Peak.r*1.2, Peak.dec+Peak.r*1.2))
    if legend==True:
        ax.legend()
    ax.invert_xaxis()
    if ax is None:
    	plt.close()

def three_cutouts()

def euclid_cutout(Peak, ax=None, annotation=True,legend=True):
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    #insert cutout code here
    
    ##
    stars = Peak.member_candidates
    ax.scatter(stars.ra,stars.dec,
	       s=50, facecolors='none', edgecolors='C0',
	       label='Candidate member stars')
    if annotation==True:
        ax.scatter(Peak.ra, Peak.dec,
    	       marker='X',c='C1',
    	       label=f'Peak center ({round(Peak.ra,2)},{round(Peak.dec,2)} deg) and radius = {round(Peak.r,2)}')
        radius = Annulus((Peak.ra, Peak.dec),Peak.r, width, color='C1')
        ax.add_artist(radius)
    ax.set(title = "Euclid",
           xlabel='RA (deg)', ylabel='Dec (deg)',
           xlim=(Peak.ra-Peak.r*1.2, Peak.ra+Peak.r*1.2),ylim=(Peak.dec-Peak.r*1.2, Peak.dec+Peak.r*1.2))
    if legend==True:
        ax.legend()
    ax.invert_xaxis()
    if ax is None:
    	plt.close()

def des_cutout(Peak, band, ax=None, annotation=True,legend=True, access_url="https://datalab.noirlab.edu/sia/nsc_dr2", save=False, output_path=None):

    ## insert cutout code here
    service = sia.SIAService(access_url) # that default access url is to DR2
    fov = Peak.r*1.5
    imgTable = service.search((Peak.ra,Peak.dec), (fov/np.cos(Peak.dec*np.pi/180), fov), verbosity=2).to_table()
    table = imgTable[(imgTable['prodtype'] == 'image') & (imgTable['obs_bandpass']==band) & \
                    (np.char.find(np.ma.filled(imgTable['object'].astype(str), fill_value='none'),'DES')!=-1)]
    row = table[np.argmin(table['seeing'])] #does this make sense??
    url = row['access_url']
    filename = download_file(url,cache=True,show_progress=False,timeout=120)
    hdu = fits.open(filename)[0]
    image = hdu.data
    hdr = hdu.header
    wcs = WCS(hdr)
    # ax being None means that this is a standalone plot
    if ax is None:
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
    ax.imshow(image, cmap='gray', vmin=image.min(), vmax=image.min()+(image.max()-image.min())/100.)
    ##
    
    stars = Peak.member_candidates
    ax.scatter(stars.ra,stars.dec, transform=ax.get_transform('icrs'), s=1000./(stars.i.mag-12), 
               facecolors='none', edgecolors='C0', linewidths=3,
               label='Candidate member stars')
    if annotation==True:
        ax.scatter(Peak.ra, Peak.dec,
    	       marker='X',c='C1',
    	       label=f'Peak center ({round(Peak.ra,2)},{round(Peak.dec,2)} deg) and radius = {round(Peak.r,2)}')
        radius = Annulus((Peak.ra, Peak.dec),Peak.r, width, color='C1')
        ax.add_artist(radius)
    ax.set(title = "DES Legacy DR9",
           xlabel='RA (deg)', ylabel='Dec (deg)',
           xlim=(Peak.ra-Peak.r*1.2, Peak.ra+Peak.r*1.2),ylim=(Peak.dec-Peak.r*1.2, Peak.dec+Peak.r*1.2))
    if legend==True:
        ax.legend()
    ax.invert_xaxis()
    if save == True:
        if output_path is None:
            print('You need to specify path in order to save)
            plt.close()
            return
        plt.savefig(output_path)
    # ax being None means this is standalone therefore we can close it now
    if ax is None:
    	plt.close()
        
def lsst_cutout(Peak, ax=None, annotation=True,legend=True):
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    #insert cutout code here
    
    ##
    stars = Peak.member_candidates
    ax.scatter(stars.ra,stars.dec,
	       s=50, facecolors='none', edgecolors='C0',
	       label='Candidate member stars')
    if annotation==True:
        ax.scatter(Peak.ra, Peak.dec,
    	       marker='X',c='C1',
    	       label=f'Peak center ({round(Peak.ra,2)},{round(Peak.dec,2)} deg) and radius = {round(Peak.r,2)}')
        radius = Annulus((Peak.ra, Peak.dec),Peak.r, width, color='C1')
        ax.add_artist(radius)
    ax.set(title = "LSST",
           xlabel='RA (deg)', ylabel='Dec (deg)',
           xlim=(Peak.ra-Peak.r*1.2, Peak.ra+Peak.r*1.2),ylim=(Peak.dec-Peak.r*1.2, Peak.dec+Peak.r*1.2))
    if legend==True:
        ax.legend()
    ax.invert_xaxis()
    if ax is None:
    	plt.close()
    

#~~~~~~~~~~START MAPPING FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def map_plot(hsp_map, title, color_lims = (24,26),
	     save = True, filename = '', ax = None):
    if ax is None:
    	fig, ax = plt.subplots(figsize=(12, 8))
    sp = skyproj.MollweideSkyproj(ax=ax)
    sp.draw_hspmap(hsp_map, vmin = color_lims[0], vmax = color_lims[1])

    plt.title(title, pad=25)
    plt.colorbar(shrink=0.5)

    if save == True:
        if not os.path.exists(plots_dir + f'/maps'):
            os.mkdir(plots_dir + f'/maps')
        if filename == '':
            filename = title.lower.replace(' ','').replace('-','_').replace(',','_')
        plt.savefig(plots_dir + f'/maps/{filename}.png')
    if ax is None:
    	plt.close()

#~~~~~~~~~~START ISOCHRONE FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def isochrone_plot(iso, distance_modulus, 
                   isostars_band1, isostars_band2,
                   title,
                   allstars_band1=None, allstars_band2=None,
                   save = True, filename = '', ax=None):
    '''
    Plots a CMD (any bands) with isochrone line on top

    Parameters
    ----------
    iso : Isochrone object
    distance_modulus : float, converted from distance using ugali coordinate tools
    isostars_band1, isostars_band2 : Band class of the stars within isochrone cut (e.g. isostars.g)
    title : string, title for the plot (required)
    allstars_band1, allstars_band2 : default None, else Band class of the stars NOT within isochrone cut (e.g. allstars.g)
                                     if they're not None, the data will be plotted as a scatter
    save : default True, decides whether to save the file or not
    filename : default '', else will be the title with special characters and spaces removed and letters turned lowercase

    Returns
    -------
    Pretty plot, saves to plots_dir/isochrones/{filename}.png
    '''
    if ax is None:
    	fig, ax = plt.subplots(1,1, figsize=(6,6))
        
    index = np.min(np.where(iso.stage == iso.hb_stage)[0]) + 1

    ax.plot(iso.mag_1[0:index] - iso.mag_2[0:index], iso.mag_1[0:index] + distance_modulus, color='k')
    ax.plot(iso.mag_1[index:] - iso.mag_2[index:], iso.mag_1[index:] + distance_modulus, color = 'k')
    #ax.scatter(uncut_data.g_mag - uncut_data.r_mag, uncut_data.g_mag, c='r', alpha = 0.3, label = 'Before cut')
    #ax.scatter(cut_data.g_mag - cut_data.r_mag, cut_data.g_mag, c='b', alpha = 0.5, label = 'After cut')
    
    if allstars_band1 is not None and allstars_band2 is not None:
        ax.scatter(allstars_band1.mag - allstars_band2.mag,
                   allstars_band1.mag, s=10, c = 'r', alpha =0.3, label = 'All stars')
    elif allstars_band1 is None and allstars_band2 is not None or allstars_band1 is not None and allstars_band2 is None:
        print('Only one band of uncut data was inputted. Input both bands to plot this data')
        
    ax.scatter(isostars_band1.mag - isostars_band2.mag,
               isostars_band1.mag,
               s=10, c = 'b', alpha =0.5,
               label = 'Stars within isochrone template')
        
    ax.set(xlabel = f'{isostars_band1.str} - {isostars_band2.str}', ylabel = f'{isostars_band1.str}', xlim = (-1,4), ylim = (28,18), title = title)
    ax.legend()

    if save == True:
        if not os.path.exists(plots_dir + f'/isochrones'):
            os.mkdir(plots_dir + f'/isochrones')
        if filename == '':
            specials = ":,.;&-*()=+/'\n"
            filename = title.replace(' ', '').replace("'",'').lower()
            for char in specials:
                filename = filename.replace(char, '_')
        plt.savefig(plots_dir + f'/isochrones/{filename}.png')
    if ax is None: 	
    	plt.close()

#~~~~~~~~~~START MATCH VERIFICATION FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def oneD_hist_matches(match1Band, unmatch1Band, full1Band, match2Band, unmatch2Band, full2Band, SearchRegion, PrimaryData, SecondaryData):
    '''
    Plots 2 one-dimensional histograms of total sources, matched sources, and unmatched sources for Euclid and LSST

    Parameters - all of which come from the merging_catalogs function
    ----------
    match(1/2)Band, unmatch(1/2)Band : Band objects of the matched and unmatched sources from 1=Primary, 2=Secondary surveys
    full(1/2)Band : Band objects of the full data (matched+unmatched) from, again, 1=Primary, 2=Secondary surveys
                Band has attributes mag, magerr, flux, fluxerr, and str
    SearchRegion : Region object, has all the attributes to define where this data is looking
    Primary/SecondaryData : some Data objects (probably a child class), in this context just tells us which data release this is

    Returns
    --------
    1 pretty plot, 2 subplots
    Saves to plots_dir/merged_verification/1dhistmatches_{nside}_{pixel}_{primary survey}_{secondary survey}.png
    '''

    fig, ax = plt.subplots(1,2, figsize=(20,7))
    b = 60

    #one dataset...
    ax[0].hist(match1Band.mag, bins = b, histtype = 'step', color='b', label = 'Matched')
    ax[0].hist(unmatch1Band.mag, bins = b, histtype = 'step', color='r', label = 'Unmatched')
    ax[0].hist(full1Band.mag, bins = b, histtype = 'step', color='k', label = 'All Sources')
    ax[0].set(xlabel=f'{match1Band.str} mag', xlim=(16, 36), ylabel='Number counts', yscale='log',
              title=f'{PrimaryData.release.replace('_',' ').upper()} Sources Match/Unmatch')

    #... then the other
    ax[1].hist(match2Band.mag, bins = b, histtype = 'step', color='b', label = 'Matched')
    ax[1].hist(unmatch2Band.mag, bins = b, histtype = 'step', color='r', label = 'Unmatched')
    ax[1].hist(full2Band.mag, bins = b, histtype = 'step', color='k', label = 'All Sources')
    ax[1].set(xlabel=f'{match2Band.str} mag', xlim=(16, 36), ylabel='Number counts', yscale='log',
              title=f'{SecondaryData.release.replace('_',' ').upper()} Sources Match/Unmatch')

    plt.legend()

    plt.savefig(plots_dir + '/merged_verification/' + f'''1dhistmatches_
                                                          {SearchRegion.nside}_{SearchRegion.pixel}_
                                                          {PrimaryData.release}_{SecondaryData.release}.png''')
    plt.close()

def twoD_hist_matches(SearchRegion, matchPrimaryData, matchSecondaryData):
    '''
    Plots 2 two-dimensional histograms of matched sources for 2 surveys

    Parameters - all of which come from the merging_catalogs function
    ----------
    SearchRegion : Region object, has all the attributes to define where this data is looking
    Primary/SecondaryData : some Data objects (probably a child class) to map the locations of matched sources

    Returns
    --------
    Pretty plot with 2 subplots (one side is the matches in Euclid, other is the matches in LSST)
    Saves to plots_dir/merged_verification/2dhist_{nside}_{pixel}_{primary survey}_{secondary survey}.png
    '''
    fig, ax = plt.subplots(1,2, figsize=(13,5))
    _, _, _, im = ax[0].hist2d(matchPrimaryData.ra, matchPrimaryData.dec, bins=100)
    plt.colorbar(im, ax=ax[0])
    ax[0].set(title = f'Matches in {matchPrimaryData.release.replace('_',' ').upper()}', ylabel = "Dec (deg)", xlabel = "RA (deg)")
    ax[0].invert_xaxis()
    _, _, _, im = ax[1].hist2d(matchSecondaryData.ra, matchSecondaryData.dec, bins=100)
    plt.colorbar(im, ax=ax[1])
    ax[1].set(title = f'Matches in {matchSecondaryData.release.replace('_',' ').upper()}', ylabel = "Dec (deg)", xlabel = "RA (deg)")
    ax[1].invert_xaxis()
    plt.savefig(plots_dir + '/merged_verification/' + '2dhist_' + f'{SearchRegion.nside}_{SearchRegion.pixel}_{matchPrimaryData.release}_{matchSecondaryData.release}.png')
    plt.close()


def source_scatterplot(PrimaryData, SecondaryData, SearchRegion, mag_cut=None, mag_cut_label=None, ra_limits=None, dec_limits = None):
    '''
    Plots 2-D scatterplot of where *all* the source coordinates of the two surveys are, within 0.01x0.01 deg box

    Parameters - all of which come from the merging_catalogs function
    ----------
    Primary/SecondaryData : some Data objects (probably a child class) to map the locations of all sources, to see if they overlap
    SearchRegion : Region object, has all the attributes to define where this data is looking
    mag_cut (optional) : default None, else mask. A mask of a magnitude cut to be applied, e.g. (i.mag > 24)
    mag_cut_label (optional) : default None, else string. Labels what is the magnitude cut for the title
    ra_limits, dec_limits (optional) : default None, else a tuple of floats. The coordinate range you wish to plot (in degrees) (ra should be given backwards if a flipped x-axis is desired)

    Returns
    --------
    Pretty plot
    Saves to plots_dir/merged_verification/2dscatter_{nside}_{pixel}_{primary survey}_{secondary survey}_{mag_cut_label}.png
    '''
    title = f'{PrimaryData.release.replace('_',' ').upper()} and {SecondaryData.release.replace('_',' ').upper()} before matching'
    filename = '2dscatter_' + f'{SearchRegion.nside}_{SearchRegion.pixel}_{PrimaryData.release}_{SecondaryData.release}'
    if mag_cut is not None:
        primarydata = PrimaryData.apply_mask(mag_cut)
        secondarydata = SecondaryData.apply_mask(mag_cut)
        title += f' with magnitude cut {mag_cut_label}'
        filename += f'_{mag_cut_label.maketrans({'<': 'lt', '>': 'gt', ' ': ''})}.png'
    else:
        primarydata = PrimaryData
        secondarydata = SecondaryData
        filename += '.png'

    plt.scatter(primarydata.ra, primarydata.dec,
                marker = '+', label = f'{primarydata.release.replace('_',' ').upper()}'
               )
    plt.scatter(secondarydata.ra, secondarydata.dec,
                marker = 'x', label = f'{secondarydata.release.replace('_',' ').upper()}'
               )
    if ra_limits is None:
        ra_limits = (np.median(primarydata.ra)+0.01, np.median(primarydata.ra))
        dec_limits = (np.median(primarydata.dec),np.median(primarydata.dec)+0.01)
    plt.xlim(ra_limits[0], ra_limits[1])
    plt.xlabel('RA (deg)')
    plt.ylim(dec_limits[0], dec_limits[1])
    plt.ylabel('DEC (deg)')
    #plt.colorbar()
    plt.legend()
    plt.title(title)
    plt.savefig(plots_dir + '/merged_verification/' + filename)
    plt.close()

def coord_diff_hist(merged_df_coord1, merged_df_coord2, PrimaryData, SecondaryData, SearchRegion, ds=None):
    '''
    Plots 2 one-dimensional histograms of the RA and Dec differences (Euclid - LSST for both)

    Parameters - all of which come from the merging_catalogs function
    ----------
    merged_df : Table or other dataframe type, merged catalog (both Euclid and LSST data)
    tract : int, tract these datasets are querying
    Primary/SecondaryData : some Data objects (probably a child class) just need the releases
    SearchRegion : Region object, has all the attributes to define where this data is looking

    Returns
    --------
    1 pretty plot, 2 subplots
    Saves to plots_dir/merged_verification/{coorddiff}_{nside}_{pixel}_{primary survey}_{secondary survey}.png
    '''
    #commented out version of this below
    '''
    #histogram of separation
    ds = ds * 3600 #ds is in degrees, want to plot in arcsecs
    #I forced in the function that matches would be <1"
    plt.hist(ds, histtype='step', range=(0,1))
    plt.xlabel('separation [arcsec]')
    plt.title(f'Tract {tract}: Matched Source Separation')
    plt.tight_layout()
    plt.show()
    '''
    fig, ax = plt.subplots(1,2, figsize=(20,7))
    ra_diff = (merged_df_coord1[0] - merged_df_coord2[0])*3600
    dec_diff = (merged_df_coord1[1] - merged_df_coord2[1])*3600
    
    ax[0].hist(ra_diff, bins = 100)
    ax[0].set(title='Merged Catalog RA difference', xlabel=f'{PrimaryData.release} RA - {SecondaryData.release} RA (arcsec)', xlim=(-1,1))
    
    ax[1].hist(dec_diff, bins = 100)
    ax[1].set(title='Merged Catalog DEC difference', xlabel=f'{PrimaryData.release} DEC - {SecondaryData.release} DEC (arcsec)', xlim=(-1,1))
    
    plt.savefig(plots_dir + '/merged_verification/' + 'coorddiff_' + f'''{SearchRegion.nside}_{SearchRegion.pixel}
                                                                        _{PrimaryData.release}_{SecondaryData.release}.png''')
    plt.close()


def match_validation_plots(match1Band, unmatch1Band, full1Band,
                           match2Band, unmatch2Band, full2Band,
                           merged_df_coord1, merged_df_coord2,
                           matchPrim, matchSecun,
                           SearchRegion, PrimaryData, SecondaryData):
    '''
    All the outputs from the merging_catalogs function going into different validation plots

    Parameters - all of which come from the merging_catalogs function
    ----------
    tract : int, tract these datasets are querying
    lsst_survey : str, label of the survey release from which the data comes, e.g. 'dp2'
    euclid_survey : str, label of the survey release from which the data comes, e.g. 'q1'
    merged_df : Table or other dataframe type, merged catalog (both Euclid and LSST data)
    matches_euclid, unmatched_euclid : numpy arrays of the matched and unmatched sources from Euclid
    matches_lsst, unmatched_lsst : numpy arrays of the matched and unmatched sources from LSST
    lsst_table : Table or other dataframe type, all LSST data from that area 
    euclid_field : Table or other dataframe type, all Euclid data from that area
    ds : degree separation of each source, generated by ugali match function

    Returns
    --------
    Bunch of pretty plots
    Saves to plots_dir/merged_verification/
    '''
    if not os.path.exists(plots_dir + f'/merged_verification'):
        os.mkdir(plots_dir + f'/merged_verification')

    oneD_hist_matches(match1Band, unmatch1Band, full1Band, match2Band, unmatch2Band, full2Band,
                      SearchRegion, PrimaryData, SecondaryData)
    twoD_hist_matches(SearchRegion, matchPrim, matchSecun)
    source_scatterplot(PrimaryData, SecondaryData, SearchRegion, mag_cut=None, mag_cut_label=None)
    coord_diff_hist(merged_df_coord1, merged_df_coord2, PrimaryData, SecondaryData, SearchRegion, ds=None)
    print('Match validation plots ran and saved')


#~~~~~~~~~~START COLOR-MAG FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def color_magnitude(band1_mag, band1_str,
                    band2_mag, band2_str,
                    color_data, title,
                    color_label = '', selection_label = '_nolegend',
                    colors = 'viridis',
                    histogram = False, ax = None,
                    colorbar_limits = (0,1),
                    x_lim = (-1, 4), y_lim = (30, 18),
                    save = False, filename = None):
    """
    Plots color-magnitude diagram, band1 vs band1-band2

    Parameters
    ----------
    band1_mag : array or column
        Band magnitude to be plotted on y-axis
    band1_str : string
        Label of band1
    band2_mag : array or column
        Band magnitude to be subtracted for the x-axis
    band2_str : string
        Label of band2
    color_data : data column OR None
        If plot is to be scatterplot or contour, then should be data for the mapping
        If plot is to be histogram, then None
    title : string
        Title of plot/subplot
    color_label (optional) : string, default ''
        Colorbar label
        If color_data is None, color_label won't be used
    selection_label (optional) : string, default '_nolegend'
        Label for legend if desired to annotate which selection was chosen for this plot
        If kept _nolegend, won't display a legend
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
    Pretty plot, saves to plots_dir/colormag/{filename}.png
    """

    if ax == None: #to allow me to have this as a subplot
        fig, axes = plt.subplots(1,1, figsize=(7,5))
        ax = axes

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
        ax.set(xlabel = f"{band1_str} - {band2_str}", ylabel = f"{band1_str}", xlim = x_lim, ylim = y_lim)
        # cleaning up, not displaying everything if not needed
        if selection_label[0] != '_':
            ax.legend()
        if color_label is not None:
            plt.colorbar(_,label=color_label)
        plt.tight_layout()

        if save == True:
            if not os.path.exists(plots_dir + f'/colormag'):
                os.mkdir(plots_dir + f'/colormag')
            if filename is None:
                filename = title.replace(' ', '').replace('-', '_').lower()
            plt.savefig(plots_dir + f'/colormag/{filename}.png')
        plt.close()

#~~~~~~~~~~START COLOR-COLOR FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def color_color(band_list,
                color_data, title,
                color_label = '', selection_label = '_nolegend',
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
    title : str
        Plot title
    color_label (optional) : str, default ''
        Colorbar label. If colors is None, color_label won't be used
    selection_label (optional) : str, default '_nolegend'
        Labels how the data was cut, to be displayed in legend
        If kept _nolegend, won't display a legend
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
    Pretty plot, saves to plots_dir/colorcolor/{filename}.png
    """

    if ax == None:
        fig, axes = plt.subplots(1,1, figsize=(7,6))
        ax = axes

    band_x1, band_x2, band_y1, band_y2 = band_list
    #cleaning out nans because the hist2d doesn't like them
    nanmask = (~np.isnan(band_x1[1])) & (~np.isnan(band_x2[1])) & (~np.isnan(band_y1[1])) & (~np.isnan(band_y2[1]))
    x1_mag = band_x1[1][nanmask]
    x2_mag = band_x2[1][nanmask]
    y1_mag = band_y1[1][nanmask]
    y2_mag = band_y2[1][nanmask]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='.*colormapping.*')
        warnings.filterwarnings("ignore", message='.*labels.*')
        if histogram == True:
            if color_data is None:
            # 2D histograms require colormap which is str type
                ax.hist2d(x1_mag - x2_mag, y1_mag - y2_mag, bins=200,
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
        if color_label != '':
            plt.colorbar(_,label=color_label)
        if y_lim is not None:
            ax.set(ylim=y_lim)
        if x_lim is not None:
            ax.set(xlim=x_lim)
        plt.tight_layout()

        if save == True:
            if not os.path.exists(plots_dir + f'/colorcolor'):
                os.mkdir(plots_dir + f'/colorcolor')
            if filename is None:
                title = title.replace(' ', '').replace('-', '_').lower()
                filename = title
            plt.savefig(plots_dir + f'/colorcolor/{filename}.png')

#~~~~~~~~~~START STAR-GAL MORPHOLOGY FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def star_gal_sep(xax_mag, yax_sep, y_label,
                 color_data, title, 
                 colors = 'viridis_r', x_label = 'LSST i mag', 
                 selection_label = '_nolegend', color_label = None,
                 histogram = False, ax = None, line_plt = None,
                 colorbar_limits = (0,1), y_lim = None, x_lim = (18,28),
                 save = False, filename = None):
    '''
    Plots the morphology across magnitudes to show star-galaxy classifiers
    Shows magnitude where things start getting confused

    Parameters
    ----------
    xax_mag : array or column
        Data for the x-axis, typically an i magnitude
        Assumed to be LSST i mag; if not, change x_label
    yax_sep : array or column
        Data for the y-axis, some morphology measurement
    y_label : string
        Descriptor of the morphology measurement to label y-axis
    colors : data column OR None
        If plot is to be scatterplot or contour, then should be data for the mapping
        If plot is to be histogram, then None
    title : str
        Plot title
    colors (optional) : str, default 'viridis_r'
        Color scheme to use
    x_label : str, default 'LSST i mag'
        Descriptor to label x-axis
    selection_label (optional) : str, default '_nolegend'
        Description of how the data was cut to be displayed in legend
        If kept _nolegend, won't display a legend
    color_label (optional) : str, default None
        Colorbar label
        If None, won't display a colorbar
    histogram (optional) : boolean, default False
        If True, will plot either histogram or contour (depending on colors type)
        If False, will plot scatterplot
    ax (optional) : default None, else plt axes object
        Axes to plot subplots on
        If none, will set subplot axes (1,1)
    colorbar_limits (optional) : tuple, default (0, 1)
        Limits for the colorbar (changes vmin, vmax)
    y_lim (optional) : tuple, default None
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
    Pretty plot, saves to plots_dir/stargalsep/{filename}.png
    '''
    if ax == None:
        fig, axes = plt.subplots(1,1, figsize=(7,5))
        ax = axes

    # the code will start yelling about the colorbar and labels depending on how you plot it
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='.*colormapping.*')
        warnings.filterwarnings("ignore", message='.*labels.*')
        # if you want to display by histogram/contour plot bc scatter gets oversaturated
        if histogram == True:
            # 2D histograms require colormap which is str type
            if color_data is None:
                ax.hist2d(xax_mag, yax_sep, bins=200,
                          range = [[x_lim[0],x_lim[1]],[y_lim[0], y_lim[1]]],
                          cmin=1, cmap = colors, label = selection_label)
            # otherwise if color dimension is needed for data, it'll plot contours over scatter
            else:
                _ = ax.scatter(xax_mag, yax_sep,
                               c = color_data, vmin=colorbar_limits[0], vmax=colorbar_limits[1],
                               label = selection_label,
                               s = size, cmap = colors)
                sns.kdeplot(x=xax_mag, y=yax_sep, fill=False, color="k")
        # normal scatter plot
        else:
            _ = ax.scatter(xax_mag, yax_sep,
                           c = color_data, vmin=colorbar_limits[0], vmax=colorbar_limits[1],
                           label = selection_label,
                           s = size, cmap = colors)

        ax.set(xlabel = x_label, ylabel = y_label, xlim = x_lim)
        ax.set_title(title, pad=pad) #title is separate so I can have the pad
        # cleaning up the place, not displaying everything if not needed
        if y_lim is not None:
            ax.set(ylim = y_lim)
        if selection_label[0] != '_':
            ax.legend()
        if line_plt is not None:
            ax.plot(line_plt[0],line_plt[1],line_plt[2])
        if color_label is not None:
            plt.colorbar(_,label=color_label)
        plt.tight_layout()

        if save == True:
            if not os.path.exists(plots_dir + f'/stargalsep'):
                os.mkdir(plots_dir + f'/stargalsep')
            if filename is None:
                title.replace(' ', '').replace('-','_').lower()
                filename = title
            plt.savefig(plots_dir + f'/stargalsep/{filename}.png')




## ~~~~~ very outdated function trying to automate doing CMD subplots ~~~~
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
