"""Runs commands to produce convolved predicted counts map in current directory
"""

from gammapy.spectral_cube import compute_npred_cube, GammaSpectralCube, convolve_npred_cube
from astropy.units import Quantity
from gammapy.datasets.load import get_fermi_diffuse_background_model
from gammapy.datasets import FermiVelaRegion
from gammapy.stats import significance
from gammapy.image.utils import disk_correlate
from astropy.io import fits
from astropy.wcs import WCS
from astropy.units import Quantity
import numpy as np
import matplotlib.pyplot as plt

background_file = get_fermi_diffuse_background_model()
exposure_file = FermiVelaRegion.filenames()['exposure_cube']
counts_file = FermiVelaRegion.filenames()['counts']

background_model = GammaSpectralCube.read(background_file)
exposure_cube = GammaSpectralCube.read(exposure_file)
repro_bg_cube = background_model.reproject_to(exposure_cube)
energies = Quantity([10000, 10000000], 'MeV')
npred_cube = compute_npred_cube(repro_bg_cube, exposure_cube, energies)
convolved_npred_cube = convolve_npred_cube(npred_cube, 3, 0.1)

counts_data = fits.open(counts_file)[0].data
counts_wcs = WCS(fits.open(counts_file)[0].header)
counts_cube = GammaSpectralCube(data = Quantity(counts_data, ''), wcs = counts_wcs, energy = energies)
counts_cube = counts_cube.reproject_to(npred_cube)

counts = np.nan_to_num(counts_cube.data[0])
model = np.nan_to_num(convolved_npred_cube.data[0])

correlation_radius = 3
correlated_counts = disk_correlate(counts, correlation_radius)
correlated_model = disk_correlate(model, correlation_radius)

significance = significance(correlated_counts, correlated_model, method='lima')

fig = plt.imshow(significance)
fig.colorbar
fig.axes.set_visible('False')

# Then plots results

