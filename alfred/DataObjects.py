from alfred import utils

class Data():
    def __init__(self, data):
        self.data = data

    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return Data(new_data)

class Band():
    def __init__(self, flux, fluxerr, name):
        self.flux = flux
        self.fluxerr = fluxerr
        self.mag = utils.flux2mag(flux)
        self.magerr = utils.fluxerr2magerr(flux, fluxerr)
        self.str = name

class LSSTData(Data):
    def __init__(self, data, lsst_release):
        super(LSSTData, self).__init__(data)
        self.release = lsst_release
        #self.tract = tract
        #self.field = utils.get_field(tract)
        
        ## coordinates
        self.ra_limits = (data['coord_ra'].min(), data['coord_ra'].max())
        self.dec_limits = (data['coord_dec'].min(), data['coord_dec'].max())
        self.ra = data['coord_ra']
        self.dec = data['coord_dec']

        ## Rubin bands
        self.g = Band(data['g_psfFlux'], data['g_psfFlux'], 'g')
        self.r = Band(data['r_psfFlux'], data['r_psfFlux'], 'r')
        self.i = Band(data['i_psfFlux'], data['i_psfFlux'], 'i')
        self.z = Band(data['z_psfFlux'], data['z_psfFlux'], 'z')
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
    def band_psfmincmodel(band, self):
        psf_flux = utils.flux2mag(self.data[f'{band}_psfFlux'])
        cmodel_flux = utils.flux2mag(self.data[f'{band}_cModelFlux'])
        return psf_flux - cmodel_flux
    def band_psfdivcmodel(band, self):
        psf_flux = utils.flux2mag(self.data[f'{band}_psfFlux'])
        cmodel_flux = utils.flux2mag(self.data[f'{band}_cModelFlux'])
        return psf_flux / cmodel_flux
        
    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return LSSTData(new_data, self.release)

class EuclidData(Data):
    def __init__(self, data, euclid_survey):
        super(EuclidData, self).__init__(data)
        self.release = euclid_survey
        
        ## coordinates
        self.ra = data['RIGHT_ASCENSION']
        self.dec = data['DECLINATION']
    
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
    def __init__(self, data, des_survey):
        super(DESData, self).__init__(data)
        self.release = des_survey
        #self.ra = 
        #self.dec = 
        #self.... insert all the bands and morphology here

        
class LSSTnEuclidData(LSSTData, EuclidData):
    def __init__(self, merged_data, lsst_release, euclid_release, coord_choice='LSST'):
        LSSTData.__init__(self, merged_data, lsst_release)
        EuclidData.__init__(self, merged_data, euclid_release)

        if coord_choice=='LSST':
            self.ra = LSSTData.ra
            self.dec = LSSTData.dec
        else:
            self.ra = EuclidData.ra
            self.dec = EuclidData.dec
        self.coord_choice = coord_choice
        self.release = f'{LSSTData.release}_{EuclidData.release}'

    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return LSSTnEuclidData(new_data, self.release, self.release, coord_choice=self.coord_choice)

class DESnEuclidData(DESData, EuclidData):
    def __init__(self, merged_data, des_release, euclid_release, coord_choice='DES'):
        DESData.__init__(self, merged_data, des_release)
        EuclidData.__init__(self, merged_data, euclid_release)

        if coord_choice=='DES':
            self.ra = DESData.ra
            self.dec = DESData.dec
        else:
            self.ra = EuclidData.ra
            self.dec = EuclidData.dec
        self.coord_choice = coord_choice
        self.release = f'{DESData.release}_{EuclidData.release}'

    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return DESnEuclidData(new_data, coord_choice=self.coord_choice)
        
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

class Peaks(): #TO BUILD
    def __init__(self, x, y, angsep):
        self.x = x
        self.y = y
        self.angsep = angsep
        self.ra = 0
        self.dec = 0
        self.r = 0
        self.sig = 0
        self.dist = 0
        self.n_obs = 0
        self.n_obs_half = 0
        self.n_model = 0

    def compute_local_char_density():
        # takes in nside, data, characteristic density, ra, dec, mag_max, and a fracdet map (x, y needed but already attributes)
        # finds and returns local characteristic density
        return None
        
    def fit_aperature():
        #takes in a projection, distance, local characteristic density, (x, y, angsep needed but already attributes)
        #finds and returns ra_peaks, dec_peaks, r_peaks, sig_peaks, distance_moduli, n_obs_peaks, n_obs_half_peaks, n_model_peaks (makes these attributes)
        #but there also may be multiples so i'd have to think how to handle those....
        return None

#eventually I want to have an array/list of Peak objects

    

