class DataFrame():
    def __init__(self, survey, data, field, qual_cuts_dict={}):
        self.survey = survey
        self.data = data
        self.field = field
        self.quality_cuts = qual_cuts_dict
        
        ## coordinates
        self.ra_limits = (data['coord_ra'].min(), data['coord_ra'].max())
        self.dec_limits = (data['coord_dec'].min(), data['coord_dec'].max())
        self.rubin_ra = data['coord_ra']
        self.rubin_dec = data['coord_dec']
        self.euclid_ra = data['RIGHT_ASCENSION']
        self.euclid_dec = data['DECLINATION']

        ## Euclid bands (flux given in mu_Jy)
        num = 2 #as suggested in Zerjal et al
        self.H_mag = flux2mag(data[f'FLUX_H_{num}FWHM_APER']*(10**3)) #convert to nJy
        self.Y_mag = flux2mag(data[f'FLUX_Y_{num}FWHM_APER']*(10**3)) #convert to nJy
        self.J_mag = flux2mag(data[f'FLUX_J_{num}FWHM_APER']*(10**3)) #convert to nJy
        self.H_magerr = fluxerr2magerr(self.H_mag, flux2mag(data[f'FLUXERR_H_{num}FWHM_APER']*(10**3)))
        self.Y_magerr = fluxerr2magerr(self.Y_mag, flux2mag(data[f'FLUXERR_Y_{num}FWHM_APER']*(10**3)))
        self.J_magerr = fluxerr2magerr(self.J_mag, flux2mag(data[f'FLUXERR_J_{num}FWHM_APER']*(10**3)))

        ## Rubin bands
        self.g_mag = flux2mag(data['g_psfFlux'])
        self.r_mag = flux2mag(data['r_psfFlux'])
        self.i_mag = flux2mag(data['i_psfFlux'])
        self.z_mag = flux2mag(data['z_psfFlux'])
        self.g_magerr = fluxerr2magerr(self.g_mag, flux2mag(data['g_psfFluxErr']))
        self.r_magerr = fluxerr2magerr(self.r_mag, flux2mag(data['r_psfFluxErr']))
        self.i_magerr = fluxerr2magerr(self.i_mag, flux2mag(data['i_psfFluxErr']))
        self.z_magerr = fluxerr2magerr(self.z_mag, flux2mag(data['z_psfFluxErr']))

        self.pointlikeprob = data['POINT_LIKE_PROB']
        self.ellipticity = data['ELLIPTICITY']
        self.mumax_minus_mag = self.data['MUMAX_MINUS_MAG']

    def Zerjal_cut(self):
        stars_mask = (self.ellipticity < 0.2) 
        stars_mask &= (self.mumax_minus_mag< -2.7) & (self.mumax_minus_mag> -3.25)
        return stars_mask, self.apply_mask(stars_mask)

    def POINT_LIKE_PROB_cut(self, threshold):
        stars_mask = (self.pointlikeprob > threshold)
        return stars_mask, self.apply_mask(stars_mask)

    def POINT_LIKE_FLAG_cut(self):
        stars_mask = (self.data['POINT_LIKE_FLAG'] == 1)
        return stars_mask, self.apply_mask(stars_mask)

    def apply_mask(self, mask):
        ## takes in a mask, applies it to the df, then returns another DataFrame object
        new_data = self.data[mask]
        return DataFrame(self.survey, new_data, self.field, self.quality_cuts)
        
    def clean_data(self, lsst_bands, euclid_bands, lsst_sn_bands, euclid_sn_bands, snr=5, flags = True, det_qual_flags = [0, 2]):
        '''
        inputs
        --------
        lsst_bands : str, bands (from ugrizy) in which to enforce quality flags
        euclid_bands : list of str, bands (from VIS, Y, J, H) in which to enforce quality flags
        lsst_sn_bands : str, bands (from ugrizy) in which to enforce source to noise cuts
        euclid_sn_bands : list of str, bands (from VIS, Y, J, H) in which to enforce source to noise cuts
        snr : the source to noise ratio you want to enforce

        returns
        --------
        new DataFrame object but with the cuts applied and an updated dictionary to keep track
        '''
        df = self.data
        ## 1. flags
        if flags == True:
            mask = (df['detect_isIsolated'] == True) #lsst
            mask &= (df['SPURIOUS_FLAG'] == 0) #euclid
            det_qual_mask = (df['DET_QUALITY_FLAG'] == 0)
            for i in det_qual_flags:
                if i == 0:
                    pass
                else:
                    det_qual_mask |= (df['DET_QUALITY_FLAG'] == i)
            mask &= det_qual_mask
        ## 2. per band S/N
        for band in lsst_sn_bands:
            mask &= (df[f'{band}_psfFlux']/df[f'{band}_psfFluxErr'] > snr)
        for band in euclid_sn_bands:
                #if band == 'VIS':
            #    mask &= (merged_df[f'FLUX_VIS_PSF']/merged_df[f'FLUXERR_VIS_PSF'] > snr)
            #else:
            ## I actually don't know which one to use for VIS
            num = 2
            mask &= (df[f'FLUX_{band}_{num}FWHM_APER']/df[f'FLUXERR_{band}_{num}FWHM_APER'] > snr)
        ## 3. per band flag
        for band in lsst_bands:
            mask &= (df[f'{band}_psfFlux_flag'] == 0)
            
        ## see https://euclid.esac.esa.int/dr/q1/dpdd/merdpd/merphotometrycookbook.html abt which flags to enforce
        #for band in euclid_bands:
        #    mask &= (merged_df[f'FLAG_{band}']==0)
        
        ## 4. enforcing that euclid sources having FWHM above 1.5" is spurious (Zerjal et al. suggestion)
        f = 1.5
        mask &= (df['FWHM'] <= f)

        ## updates the quality_cuts attribute along the way, so I can keep track of these things
        cut_dict = {'detect_isIsolated' : True, 'SPURIOUS_FLAG' : 0, 'DET_QUALITY_FLAG' : det_qual_flags, 
                    'LSST S/N' : f'{lsst_sn_bands} > {snr} ', 'Euclid S/N' : f'{euclid_sn_bands} > {snr} ', 
                    'LSST per band flag' : lsst_bands, 'FWHM' : f}
        if len(self.quality_cuts.keys()) != 0:
            for key in self.quality_cuts.keys():
                cut_dict[key] = cut_dict[key] + self.quality_cuts[key]

        return DataFrame(self.survey, df[mask], self.field, cut_dict)
