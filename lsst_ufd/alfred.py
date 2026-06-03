import yaml
import os
from lsst.daf.butler import Butler
from star_gal_sep import *
from utils import *
from sg_diagnostic_plots import *
from isochrone_test import *
import simple_adl.simple_adl.isochrone as isochrone
import simple_adl.simple_adl.coordinate_tools as coordinate_tools
#-----------------------------------

#this is a bit hard coded too but idk another work around
#if alfred stays in a folder below the config always then this should work
with open('../config.yaml', 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)

    home_dir_str = cfg['setup']['home_dir']
    home_dir = os.path.expandvars(home_dir_str)
    pckg_dir = os.path.join(home_dir, cfg['setup']['pckg_dir'])
    results_dir = os.path.join(home_dir, cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    plots_dir = os.path.join(home_dir, cfg['output']['plots_dir'])
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)
    
    survey = cfg['survey']
    repo_config = cfg[survey]['repo_config']
    collection = cfg[survey]['collection']
    field2tract_dict = cfg[survey]['field2tract_dict']
    skymap = cfg[survey]['skymap']
    #also have instrument, tract_list, tract2field_dict information

#initiate the butler instance
butler = Butler(repo_config, collections=collection)

INCOLS = [
    'coord_ra',
    'coord_dec',
    'detect_isIsolated',
]
bands='griz'
for band in bands:
    INCOLS += [
        f'{band}_psfFlux',
        f'{band}_cModelFlux',
        f'{band}_cModelFluxErr',
        f'{band}_psfFluxErr',
        f'{band}_extendedness',
        f'{band}_psfFlux_flag'
    ]
    if survey=='dp1':
        INCOLS += [f'{band}_sizeExtendedness']
        INCOLS += [f'{band}_sizeExtendedness_flag']
    if survey=='dp2':
        INCOLS += [f'{band}_model_extendedness']
# maybe inelegant but that's a problem for future kayleigh
data = butler.get('object', collections=[collection],
                  dataId={'skymap': skymap, 'tract': 
                          field2tract_dict['EDFS'][0]}, 
                  parameters={'columns':INCOLS})
# just calling up one tract for now
# thank you alfred

#clean up data
data = quality_mask(data, snr=5)

#separate stars from galaxies
stars = stellar_catalog(data, survey, 'fluxratioerr', 'i', 0.3, c=1.2)

#call up diagnostic plots
color_magnitude(stars,
                'g', 'r', stars['i_sizeExtendedness'], 'i_sizeExtendedness',
                'Test Plot', save = True, plots_path = plots_dir)
color_magnitude2(stars,
                 {'band1' : 'g', 'band2' : 'r', 'colors' : stars['i_sizeExtendedness'], 
                  'color_label' : 'i_sizeExtendedness', 'title' : 'i_sizeExtendedness'},
                 {'band1' : 'g', 'band2' : 'r', 'colors' : stars['g_sizeExtendedness'], 
                  'color_label' : 'g_sizeExtendedness', 'title' : 'g_sizeExtendedness'},
                 title = 'Test Plot 2', save = True, plots_path = plots_dir)
color_color(stars,'g','r','r','i', stars['i_sizeExtendedness'], 'i_sizeExtendedness',
            'ColorColor Test', save = True, plots_path = plots_dir)
star_gal_sep(data, 'LSST', data['i_sizeExtendedness'], 'i_sizeExtendedness',
             'DP1 SG Sep', save = True, plots_path = plots_dir)

#isochrone cut test
distance = 300 # kpc
#distance_modulus = ugali.utils.projector.distanceToDistanceModulus(distance)
distance_modulus = coordinate_tools.distanceToDistanceModulus(distance)

iso = isochrone.Isochrone(
        age=12.0,
        metallicity=0.0002,
        distance_modulus=distance_modulus,
        survey= 'lsst',
        band_1= 'g',
        band_2= 'r')
g = flux2mag(data['g_psfFlux'])
r = flux2mag(data['r_psfFlux'])
g_err = flux2mag(data['g_psfFluxErr'])
r_err = flux2mag(data['r_psfFluxErr'])
cut = cut_isochrone_path(g, r, g_err, r_err, iso)
isocut_data = data[cut]

fig, ax = plt.subplots(1,1, figsize=(6,6))
index = np.min(np.where(iso.stage == iso.hb_stage)[0]) + 1
ax.set(xlabel = 'g-r', ylabel = 'g', xlim = (-1,4), ylim = (28,18))
ax.plot(iso.mag_1[0:index] - iso.mag_2[0:index], iso.mag_1[0:index] + distance_modulus)
ax.plot(iso.mag_1[index:] - iso.mag_2[index:], iso.mag_1[index:] + distance_modulus)
ax.scatter(flux2mag(isocut_data['g_psfFlux']) - flux2mag(isocut_data['r_psfFlux']), 
           flux2mag(isocut_data['g_psfFlux']))

plt.savefig(plots_dir + '/iso_test.png')
print(plots_dir + '/iso_test.png')
