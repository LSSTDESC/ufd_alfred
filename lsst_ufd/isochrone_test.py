import sys
sys.path.append("../")
import simple_adl.simple_adl.isochrone as isochrone
import simple_adl.simple_adl.coordinate_tools as coordinate_tools
import matplotlib.pyplot as plt
import numpy as np

distance = 300 # kpc
#distance_modulus = ugali.utils.projector.distanceToDistanceModulus(distance)
distance_modulus = coordinate_tools.distanceToDistanceModulus(distance)

survey = 'LSST'
bands = ['g', 'r']

iso = isochrone.Isochrone(
        age=12.0,
        metallicity=0.0002,
        distance_modulus=distance_modulus,
        survey= survey.lower(),
        band_1= bands[0],
        band_2= bands[1])

index = np.min(np.where(iso.stage == iso.hb_stage)[0]) + 1
plt.plot(iso.mag_1[0:index] - iso.mag_2[0:index], iso.mag_1[0:index] + distance_modulus)
plt.plot(iso.mag_1[index:] - iso.mag_2[index:], iso.mag_1[index:] + distance_modulus)
plt.savefig('/global/u2/k/kexcell/ultrafaints/plots/iso_test.png')
