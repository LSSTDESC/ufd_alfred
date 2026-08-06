## NEED TO TEST THE ATTRIBUTES ARE ASSIGNED CORRECTLY    
    ## write new quality cuts functions - maybe not within the classes but just using the apply_mask
from alfred.utils import *

class Data():
    def __init__(self, data):
        self.data = data
        
    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return Data(new_data)


class LSSTData(Data):
    def __init__(self, data, lsst_release, field):
        super(LSSTData, self).__init__(data)
        self.lsst_survey = lsst_release
        self.field = field
        
        ## coordinates
        self.ra_limits = (data['coord_ra'].min(), data['coord_ra'].max())
        self.dec_limits = (data['coord_dec'].min(), data['coord_dec'].max())
        self.rubin_ra = data['coord_ra']
        self.rubin_dec = data['coord_dec']

        ## Rubin bands
        self.g_mag = flux2mag(data['g_psfFlux'])
        self.r_mag = flux2mag(data['r_psfFlux'])
        self.i_mag = flux2mag(data['i_psfFlux'])
        self.z_mag = flux2mag(data['z_psfFlux'])
        self.g_magerr = fluxerr2magerr(self.g_mag, flux2mag(data['g_psfFluxErr']))
        self.r_magerr = fluxerr2magerr(self.r_mag, flux2mag(data['r_psfFluxErr']))
        self.i_magerr = fluxerr2magerr(self.i_mag, flux2mag(data['i_psfFluxErr']))
        self.z_magerr = fluxerr2magerr(self.z_mag, flux2mag(data['z_psfFluxErr']))
        
    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return LSSTData(new_data, self.lsst_release, self.field)


class LSSTnEuclidData(LSSTData):
    def __init__(self, data, lsst_survey, euclid_survey, field):
        super(LSSTnEuclidData, self).__init__(data, lsst_survey, field)
        self.euclid_survey = euclid_survey
        
        ## coordinates
        self.euclid_ra = data['right_ascension']
        self.euclid_dec = data['declination']
    
        ## Euclid bands (flux given in mu_Jy)
        num = 2 #as suggested in Zerjal et al
        self.H_mag = flux2mag(data[f'FLUX_H_{num}FWHM_APER'.lower()]*(10**3)) #convert to nJy
        self.Y_mag = flux2mag(data[f'FLUX_Y_{num}FWHM_APER'.lower()]*(10**3)) #convert to nJy
        self.J_mag = flux2mag(data[f'FLUX_J_{num}FWHM_APER'.lower()]*(10**3)) #convert to nJy
        self.H_magerr = fluxerr2magerr(self.H_mag, flux2mag(data[f'FLUXERR_H_{num}FWHM_APER'.lower()]*(10**3)))
        self.Y_magerr = fluxerr2magerr(self.Y_mag, flux2mag(data[f'FLUXERR_Y_{num}FWHM_APER'.lower()]*(10**3)))
        self.J_magerr = fluxerr2magerr(self.J_mag, flux2mag(data[f'FLUXERR_J_{num}FWHM_APER'.lower()]*(10**3)))

        ## morphology
        self.pointlikeprob = data['point_like_prob']
        self.ellipticity = data['ellipticity']
        self.mumax_minus_mag = self.data['mumax_minus_mag']
        
    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another Data object
        new_data = self.data[mask]
        return LSSTnEuclidData(new_data, self.lsst_release, self. euclid_release, self.field)
