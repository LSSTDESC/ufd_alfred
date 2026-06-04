def isochrone_search(star_data, distance, age=12.0, Z=0.0002, plots_dir):
  #distance in kpc
  distance_modulus = coordinate_tools.distanceToDistanceModulus(distance)
  
  iso = isochrone.Isochrone(
          age=age,
          metallicity=Z,
          distance_modulus=distance_modulus,
          survey= 'lsst',
          band_1= 'g',
          band_2= 'r')
  star_data['g mag'] = flux2mag(star_data['g_psfFlux'])
  star_data['r mag'] = flux2mag(star_data['r_psfFlux'])
  star_data['g mag err'] = -2.5/np.log(10)*(star_data['g_psfFluxErr']/star_data['g_psfFlux'])
  star_data['r mag err'] = -2.5/np.log(10)*(star_data['r_psfFluxErr']/star_data['r_psfFlux'])
  star_data['g mag err'][~np.isfinite(star_data['g mag err'])] = np.nan
  star_data['r mag err'][~np.isfinite(star_data['r mag err'])] = np.nan
  
  cut = cut_isochrone_path(star_data['g mag'], star_data['r mag'], star_data['g mag err'], star_data['r mag err'], iso)
  star_data = star_data[cut]
  
  fig, ax = plt.subplots(1,1, figsize=(6,6))
  index = np.min(np.where(iso.stage == iso.hb_stage)[0]) + 1
  ax.set(xlabel = 'g-r', ylabel = 'g', xlim = (-1,4), ylim = (28,18))
  ax.plot(iso.mag_1[0:index] - iso.mag_2[0:index], iso.mag_1[0:index] + distance_modulus)
  ax.plot(iso.mag_1[index:] - iso.mag_2[index:], iso.mag_1[index:] + distance_modulus)
  ax.scatter(data['g mag'] - data['r mag'], 
             data['g mag'])
  
  plt.savefig(plots_dir + '/iso_test.png')
