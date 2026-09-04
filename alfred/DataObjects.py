from alfred import utils
from ugali.utils import projector

class Data():
    def __init__(self, data, *args, **kwargs):
        self.data = data

    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return Data(new_data)

class Band():
    def __init__(self, val, valerr, name, input_type='flux'):
        if input_type=='flux':
            self.mag = utils.flux2mag(val)
            self.magerr = utils.fluxerr2magerr(val, valerr)
        elif input_type=='mag':
            self.mag = val
            self.magerr = valerr
        self.str = name

class LSSTData(Data):
    def __init__(self, data, lsst_survey='',**kwargs):
        super().__init__(data,**kwargs)
        self.release = lsst_survey
        self.survey = lsst_survey
        #self.tract = tract
        #self.field = utils.get_field(tract)

        ## coordinates
        self.ra_limits = (data['coord_ra'].min(), data['coord_ra'].max())
        self.dec_limits = (data['coord_dec'].min(), data['coord_dec'].max())
        self.ra = data['coord_ra']
        self.dec = data['coord_dec']
        self.basis1 = data['coord_ra']
        self.basis2 = data['coord_dec']

        ## Rubin bands
        self.g = Band(data['g_psfFlux'], data['g_psfFluxErr'], 'g')
        self.r = Band(data['r_psfFlux'], data['r_psfFluxErr'], 'r')
        self.i = Band(data['i_psfFlux'], data['i_psfFluxErr'], 'i')
        self.z = Band(data['z_psfFlux'], data['z_psfFluxErr'], 'z')
        #then because I have so many functions already defined, some retroactive definitions:
        self.g_mag = self.g.mag
        self.g_magerr = self.g.magerr
        self.r_mag = self.r.mag
        self.r_magerr = self.r.magerr
        self.i_mag = self.i.mag
        self.i_magerr = self.i.magerr
        self.z_mag = self.z.mag
        self.z_magerr = self.z.magerr

    ## morphology
    def band_psfmincmodel(self, band):
        psf_flux = utils.flux2mag(self.data[f'{band}_psfFlux'])
        cmodel_flux = utils.flux2mag(self.data[f'{band}_cModelFlux'])
        return psf_flux - cmodel_flux
    def band_psfdivcmodel(self, band):
        psf_flux = utils.flux2mag(self.data[f'{band}_psfFlux'])
        cmodel_flux = utils.flux2mag(self.data[f'{band}_cModelFlux'])
        return psf_flux / cmodel_flux

    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return LSSTData(new_data, self.survey)

class EuclidData(Data):
    def __init__(self, data, euclid_survey='',**kwargs):
        super().__init__(data, **kwargs)
        self.release = euclid_survey
        self.survey = euclid_survey #I realized release might be more confusing than survey? 
                                    #will try to fix where I use .release attribute

        ## coordinates
        self.ra = data['RIGHT_ASCENSION']
        self.dec = data['DECLINATION']
        self.basis1 = data['RIGHT_ASCENSION']
        self.basis2 = data['DECLINATION']

        ## Euclid bands (flux given in mu_Jy)
        num = 2 #as suggested in Zerjal et al
        #convert fluxes to nJy, that's what the flux -> mag functions assume
        self.VIS = Band(data[f'FLUX_VIS_2FWHM_APER']*(10**3),
                        data[f'FLUXERR_VIS_2FWHM_APER']*(10**3),
                        'VIS')
        self.H = Band(data[f'FLUX_H_2FWHM_APER']*(10**3),
                      data[f'FLUXERR_H_2FWHM_APER']*(10**3),
                      'H')
        self.Y = Band(data[f'FLUX_Y_2FWHM_APER']*(10**3),
                      data[f'FLUXERR_Y_2FWHM_APER']*(10**3),
                      'Y')
        self.J = Band(data[f'FLUX_J_2FWHM_APER']*(10**3),
                      data[f'FLUXERR_J_2FWHM_APER']*(10**3),
                      'J')
        #then because I have so many functions already defined, some retroactive definitions:
        self.VIS_mag = self.VIS.mag
        self.VIS_magerr = self.VIS.magerr
        self.H_mag = self.H.mag
        self.H_magerr = self.H.magerr
        self.Y_mag = self.Y.mag
        self.Y_magerr = self.Y.magerr
        self.J_mag = self.J.mag
        self.J_magerr = self.J.magerr

        ## morphology
        self.pointlikeprob = data['POINT_LIKE_PROB']
        self.ellipticity = data['ELLIPTICITY']
        self.mumax_minus_mag = self.data['MUMAX_MINUS_MAG']

    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return EuclidData(new_data, self.release)

class DESData(Data):
    def __init__(self, data, des_survey='',**kwargs):
        super(DESData, self).__init__(data,**kwargs)
        self.release = des_survey
        self.survey = des_survey
        
        ## coordinates
        self.ra_limits = (data['alphawin_j2000'].min(), data['alphawin_j2000'].max())
        self.dec_limits = (data['deltawin_j2000'].min(), data['deltawin_j2000'].max())
        self.ra = data['alphawin_j2000']
        self.dec = data['deltawin_j2000']
        self.basis1 = data['alphawin_j2000']
        self.basis2 = data['deltawin_j2000']

        ## DES bands
        self.g = Band(data['psf_mag_aper_8_g_corrected'], data['psf_mag_err_aper_8_g'], 'g', input_type='mag')
        self.r = Band(data['psf_mag_aper_8_r_corrected'], data['psf_mag_err_aper_8_r'], 'r', input_type='mag')
        self.i = Band(data['psf_mag_aper_8_i_corrected'], data['psf_mag_err_aper_8_i'], 'i', input_type='mag')
        self.z = Band(data['psf_mag_aper_8_z_corrected'], data['psf_mag_err_aper_8_z'], 'z', input_type='mag')
        self.y = Band(data['psf_mag_aper_8_y_corrected'], data['psf_mag_err_aper_8_y'], 'y', input_type='mag')
        #then because I have so many functions already defined, some retroactive definitions:
        self.g_mag = self.g.mag
        self.g_magerr = self.g.magerr
        self.r_mag = self.r.mag
        self.r_magerr = self.r.magerr
        self.i_mag = self.i.mag
        self.i_magerr = self.i.magerr
        self.z_mag = self.z.mag
        self.z_magerr = self.z.magerr
        
    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return DESData(new_data, self.survey)


class LSSTnEuclidData(LSSTData, EuclidData):
    def __init__(self, merged_data, lsst_survey='', euclid_survey='', coord_choice='LSST', **kwargs):
        LSSTData.__init__(self, data=merged_data, lsst_survey=lsst_survey, euclid_survey=euclid_survey, **kwargs)
        EuclidData.__init__(self, data=merged_data, euclid_survey=euclid_survey)

        if coord_choice=='LSST':
            self.ra = merged_data['coord_ra']
            self.dec = merged_data['coord_dec']
        else:
            self.ra = merged_data['RIGHT_ASCENSION']
            self.dec = merged_data['DECLINATION']
        self.coord_choice = coord_choice
        self.release = f'{lsst_survey}_{euclid_survey}'
        self.survey = f'{lsst_survey}_{euclid_survey}'
        self.lsst_release = lsst_survey
        self.euclid_release = euclid_survey
        self.lsst_survey = lsst_survey
        self.euclid_survey = euclid_survey

    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return LSSTnEuclidData(new_data, self.lsst_survey, self.euclid_survey, coord_choice=self.coord_choice)

class DESnEuclidData(DESData, EuclidData):
    def __init__(self, merged_data, des_survey='', euclid_survey='', coord_choice='DES',**kwargs):
        DESData.__init__(self, data=merged_data, des_survey=des_survey, euclid_survey=euclid_survey,**kwargs)
        EuclidData.__init__(self, merged_data, euclid_survey)

        if coord_choice=='DES':
            self.ra = merged_data['alphawin_j2000']
            self.dec = merged_data['deltawin_j2000']
        else:
            self.ra = merged_data['RIGHT_ASCENSION']
            self.dec = merged_data['DECLINATION']
        self.coord_choice = coord_choice
        self.release = f'{des_survey}_{euclid_survey}'
        self.survey = f'{des_survey}_{euclid_survey}'
        self.des_release = des_survey
        self.euclid_release = euclid_survey
        self.des_survey = des_survey
        self.euclid_survey = euclid_survey

    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return DESnEuclidData(new_data, self.des_survey, self.euclid_survey, coord_choice=self.coord_choice)

class Peak(): #TO BUILD
    def __init__(self, results_T):
        #results_T = ra_peak, dec_peak, r_peak, sig_peak, distance_modulus, n_obs_peak, n_obs_half_peak, n_model_peak
        self.ra = results_T[0]
        self.dec = results_T[1]
        self.r = results_T[2]
        self.sig = results_T[3]
        self.distance_modulus = results_T[4]
        self.distance = projector.distanceModulusToDistance(results_T[4])
        self.n_obs = results_T[5]
        self.n_obs_half = results_T[6]
        self.n_model = results_T[7]

        
'''
class LSSTnEuclidData(LSSTData):
    def __init__(self, data, lsst_survey, euclid_survey, field):
        super(LSSTnEuclidData, self).__init__(data, lsst_survey, field)
        self.euclid_survey = euclid_survey
        
        ## coordinates
        self.euclid_ra = data['right_ascension']
        self.euclid_dec = data['declination']
    
        ## Euclid bands (flux given in mu_Jy)
        num = 2 #as suggested in Zerjal et al
        #convert fluxes to nJy, that's what the flux -> mag functions assume
        self.VIS = Band(data[f'FLUX_VIS_{num}FWHM_APER'.lower()]*(10**3), 
                        data[f'FLUXERR_VIS_{num}FWHM_APER'.lower()]*(10**3),
                        'VIS')
        self.H = Band(data[f'FLUX_H_{num}FWHM_APER'.lower()]*(10**3), 
                      data[f'FLUXERR_H_{num}FWHM_APER'.lower()]*(10**3),
                      'H')
        self.Y = Band(data[f'FLUX_Y_{num}FWHM_APER'.lower()]*(10**3), 
                      data[f'FLUXERR_Y_{num}FWHM_APER'.lower()]*(10**3),
                      'Y')
        self.J = Band(data[f'FLUX_J_{num}FWHM_APER'.lower()]*(10**3), 
                      data[f'FLUXERR_J_{num}FWHM_APER'.lower()]*(10**3),
                      'J')
        #then because I have so many functions already defined, some retroactive definitions:
        self.VIS_mag = self.VIS.mag
        self.VIS_magerr = self.VIS.magerr
        self.H_mag = self.H.mag
        self.H_magerr = self.H.magerr
        self.Y_mag = self.Y.mag
        self.Y_magerr = self.Y.magerr
        self.J_mag = self.J.mag
        self.J_magerr = self.J.magerr
        
        ## morphology
        self.pointlikeprob = data['point_like_prob']
        self.ellipticity = data['ellipticity']
        self.mumax_minus_mag = self.data['mumax_minus_mag']
        
    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return LSSTnEuclidData(new_data, self.lsst_survey, self.euclid_survey, self.tract)
'''


    

