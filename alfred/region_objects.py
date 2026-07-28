class Tract():
    def __init__(self, tract, SkyMap):
        '''
        tract = int
        SkyMap = skyMap object, generated from butler
        '''
        self.tract = tract
        self.center = SkyCoord(SkyMap.generateTract(tract).getCtrCoord().getRa().asDegrees()*u.deg, 
                               SkyMap.generateTract(tract).getCtrCoord().getDec().asDegrees()*u.deg, 
                               frame='icrs')
        self.center_SpherePoint = SkyMap.generateTract(tract).getCtrCoord()
        ras = [SkyMap.getRaDecRange(tract)[0].asDegrees()*u.deg, SkyMap.getRaDecRange(tract)[1].asDegrees()*u.deg]
        self.ra_range = ras
        decs = [SkyMap.getRaDecRange(tract)[2].asDegrees()*u.deg, SkyMap.getRaDecRange(tract)[3].asDegrees()*u.deg]
        self.dec_range = decs
        corners = [SkyCoord(ra,dec,frame='icrs') for ra in ras for dec in decs]
        self.corners = corners
        self.corners_Angle = SkyMap.getRaDecRange(tract)

    def rubin_query(self):
        full_tract = butler.get('object', 
                                dataId={'skymap': skymap, 'tract': tract}, 
                                collections=dp2_collection, parameters={"columns":INCOLS})
        return full_tract

## CURRENTLY TRYING AND FAILING TO GET THE EUCLID QUERY WORKING
"""
    def euclid_query(self):
        '''
        returns results Table of a 1.7 deg circle around the center of tract coord 
        (which hopefully is big enough to capture the whole tract)
        '''
        
        # I don't know how big to make the Euclid query. I suppose big and then the merged catalog scales it down?
        # did the math and I think the dist from center to corner of tract is ~1.7 deg
        # assuming a patch is 13.7"
        radius = 1.7
        
        file_path = Path(my_path + f'/euclid_data/q1/{tract}_euclid_circle.parquet')
        
        if file_path.is_file():
            print('File exists.')
            results_table = Table.read(str(file_path))
        else:
            print('File not found. Querying now')
            query = '''
                    SELECT right_ascension, declination, point_like_prob, point_like_flag,
                    ellipticity, mumax_minus_mag, flux_vis_psf, fluxerr_vis_psf, spurious_flag,
                    det_quality_flag, fwhm, segmentation_map_id
                    '''
            num = 2
            for band in ['VIS', 'Y', 'J', 'H']:
                query += f", FLAG_{band}, FLUX_{band}_{num}FWHM_APER, FLUXERR_{band}_{num}FWHM_APER".lower()
            query += f" FROM mer_catalogue WHERE DISTANCE({ra}, {dec}, right_ascension, declination) < {radius}"
            #results_table = Euclid.cone_search(coordinate=coord, radius=1 * u.degree).get_results()#, columns = ['tileId'])
            
            results_table = Euclid.launch_job_async(query, verbose=False).get_results()
            results_table.write(my_path + f'/euclid_data/q1/{tract}_euclid_circle.parquet', 
                                   format='parquet', overwrite = True)
            print('File saved.')
            
        tile_ids = np.unique(results_table['segmentation_map_id'] // 10**6)  #from final catalog documentation
        self.euclid_tiles = tile_ids

        return results_table

    def merged_catalog(self):

        return table

    def euclid_maps(self, map_type):

        return euclid_healpix_map

    def rubin_maps(self, map_type):
        #look at this https://pipelines.lsst.io/py-api/lsst.analysis.tools.actions.plot.PerTractPropertyMapPlot.html

        return rubin_healsparse_map
"""
