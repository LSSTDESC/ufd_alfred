import healpy as hp
import numpy as np
import gc
import sys
import os
import yaml
from astropy.table import Table, vstack, join
#Goes up a directories to get the updated astroquery
#I really need to fix this, maybe github submodules or enforcing a version of astroquery
#I think it's version 0.4.11 ?
sys.path.append(os.path.abspath('../'))
from ugali.utils.projector import match
## function to create euclid + rubin datasets and register to the data registry on NERSC
## don't know where I want this to live quite yet

with open('config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    #assuming that it's cool that the whole github repo is considered "home"
    home_dir = os.path.expandvars(cfg['setup']['home_dir'])
    pckg_dir = os.path.join(home_dir, cfg['setup']['pckg_dir'])
    #external data is gonna be in a directory above - subject to change
    data_dir = os.path.expandvars(cfg['setup']['data_dir'])
    plots_dir = os.path.join(home_dir, cfg['output']['plots_dir'])
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)
    results_dir = os.path.join(home_dir, cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    survey = cfg['survey']
    euclid_survey = cfg['euclid_survey']

# function to check if the data doesn't exist already and if I want to rewrite it
def check_merge_data(tract, preload = True):
    '''
    preload = True means that I want to use the preloaded / saved data instead of querying again
    '''
    if not os.path.exists(data_dir + f'/merged/{tract}_{survey}_{euclid_survey}_merged.parquet'):
        #merged data file doesn't exist yet
        return True
    else:
        #merged data file DOES exist
        if preload == True:
            #I want to use the saved data, so don't remerge them
            return False
        else:
            #I want to overwrite it for whatever reason, so remerge/save them
            return True

# function to add it to the data registry

# then function to merge catalogs, starting and ending with above
def merge_catalogs(lsst_table, euclid_table, tract, preload = True, validation_needed = False):
    if not check_merge_data(tract, preload):
        print("Check tells me data exists and you don't want to remerge. Opening existing file now")
        return Table.read(data_dir + f'/merged/{tract}_{survey}_{euclid_survey}_merged.parquet')
    print('Check tells me to start the merge, starting now')

    lsst_ra, lsst_dec = lsst_table['coord_ra'], lsst_table['coord_dec']

    NSIDE=4096
    ## get the unique pixels of LSST data
    lsst_upix4096 = np.unique(hp.ang2pix(NSIDE, lsst_ra, lsst_dec, lonlat=True), return_counts=False)
    ## then get the pixels of Euclid data
    euclid_pix4096 = hp.ang2pix(NSIDE, euclid_table['right_ascension'], euclid_table['declination'], lonlat=True)
    ## Euclid has more coverage right now. We only keep the sources that lie in the LSST coverage
    spatial_mask = np.isin(euclid_pix4096, lsst_upix4096) #[lsst_cts > 8])
    euclid_field = euclid_table[spatial_mask]
    euclid_ra, euclid_dec = euclid_field['right_ascension'], euclid_field['declination']
    
    del NSIDE, lsst_upix4096, euclid_pix4096, spatial_mask, euclid_table
    gc.collect()
    
    ## match() is from ugali tools -- matching LSST and Euclid sources
    if len(euclid_ra) == 0:
        return 0
    indexlsst, indexeuclid, ds = match(lsst_ra, lsst_dec, euclid_ra, euclid_dec, tol = 0.0003)
    #print('index lsst:', '\n', indexlsst[0:20])
    #print('index euclid:', '\n', indexeuclid[0:20])
    matches_lsst = lsst_table[indexlsst]
    unmatched_lsst = lsst_table[~indexlsst]
    #print(matches_lsst.columns)
    matches_euclid = euclid_field[indexeuclid]
    unmatched_euclid = euclid_field[~indexeuclid]
    if len(matches_lsst) != len(matches_euclid):
        print("Something isn't right: those lengths don't match")
    del indexlsst, indexeuclid, lsst_ra, lsst_dec, euclid_ra, euclid_dec
    gc.collect()

    ## now merging our matches into one catalog with all LSST and Euclid columns
    matches_lsst['_match_id'] = np.arange(len(matches_lsst))
    matches_euclid['_match_id'] = np.arange(len(matches_euclid))
    merged_table = join(matches_lsst, matches_euclid, keys='_match_id')
    merged_table.write(data_dir + f'/merged/{tract}_{survey}_{euclid_survey}_merged.parquet', 
                       format='parquet', overwrite = True)

    if validation_needed==True:
        match_validation_plots(merged_table, matches_lsst, matches_euclid, 
                               unmatched_lsst, unmatched_euclid, 
                               lsst_table, euclid_field, ds)
    del matches_lsst, matches_euclid, unmatched_lsst, unmatched_euclid, lsst_table, euclid_field, ds
    gc.collect()
    return merged_table
    
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