import numpy as np

def flux_ratio(data, band, c):
    psfFlux = data[f'{band}_psfFlux']
    psfFluxErr = data[f'{band}_psfFluxErr']
    cModelFlux = data[f'{band}_cModelFlux']
    cModelFluxErr = data[f'{band}_cModelFluxErr']

    flux_ratio = psfFlux / cModelFlux
    flux_ratio_err = np.sqrt((psfFluxErr / cModelFlux)**2
                         + ((psfFlux / cModelFlux**2)*cModelFluxErr)**2)
    return (1 - flux_ratio) + c*flux_ratio_err



def stellar_catalog(data, survey, func, band, lt_threshold, c=5/2):
    '''
    accepts which classifier to use, which threshold you want the stars to be less than
    returns selected stars'''
    if survey == 'dp1':
        #i'm sure there's a better way to do this
        if func == 'fluxratioerr':
            data['starClassifier'] = flux_ratio(data, band, c)
    if survey == 'dp2':
        if func == 'model_extendedness':
            data['starClassifier'] = data[f'{band}_model_extendedness']

    stars = data[data['starClassifier'] < lt_threshold]
    return stars
