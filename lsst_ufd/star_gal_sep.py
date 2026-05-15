def stellar_catalog(survey, func, bands_list):
    '''
    accepts which classifier to use and a Data object
    returns selected stars'''
    if survey == 'dp1':
        #i'm sure there's a better way to do this
        if func == 'i_fluxratioerr':
            i_psfFlux, i_cModelFlux, i_cModelFluxErr, i_psfFluxErr = bands_list

            i_flux_ratio = i_psfFlux / i_cModelFlux
            i_flux_ratio_err = np.sqrt((i_psfFluxErr / i_cModelFlux)**2
                         + ((i_psfFlux / i_cModelFlux**2)*i_cModelFluxErr)**2)
            return (1 - i_flux_ratio) + c*i_flux_ratio_err
