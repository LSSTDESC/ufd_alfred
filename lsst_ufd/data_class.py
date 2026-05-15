import star_gal_sep

class Data:

    def __init__(self, survey, data):
        self.survey = survey
        self.ra = data['coord_ra']
        self.dec = data['coord_dec']
        self.data = data


    def star_filter(self, func, threshold, bands_list):
        '''
        Accepts a list of bands to be used for the classifier,
            that will change from survey to survey and from classifier to classifier
        Returns the masked data for only the stellar catalog
        '''
        self.data['star_clsfr'] = stellar_catalog(self.survey, func, bands_list)
        stars = self.data[self.data['star_clsfr'] < threshold]
        return stars
