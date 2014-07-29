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

background_file = get_fermi_diffuse_background_model(filename='gll_iem_v05_rev1.fit')
#background_file = FermiVelaRegion.filenames()['background_cube']
exposure_file = FermiVelaRegion.filenames()['exposure_cube']
counts_file = FermiVelaRegion.filenames()['counts']
# BACKGROUND MODEL .XML FILE SHOULD JUST HAVE THE FERMI DIFFUSE MODEL...
background_model = GammaSpectralCube.read(background_file)
exposure_cube = GammaSpectralCube.read(exposure_file)
repro_bg_cube = background_model.reproject_to(exposure_cube)
energies = Quantity([10000, 500000], 'MeV')
npred_cube = compute_npred_cube(repro_bg_cube, exposure_cube, energies)
convolved_npred_cube = convolve_npred_cube(npred_cube, 3, 1)

counts_data = fits.open(counts_file)[0].data
counts_wcs = WCS(fits.open(counts_file)[0].header)
counts_cube = GammaSpectralCube(data = Quantity(counts_data, ''), wcs = counts_wcs, energy = energies)
counts_cube = counts_cube.reproject_to(npred_cube)

counts = np.nan_to_num(counts_cube.data[0])
model = np.nan_to_num(convolved_npred_cube.data[0])

correlation_radius = 3
gtmodel = fits.open(FermiVelaRegion.filenames()['background_cube'])[0].data.astype(float)
correlated_gtmodel = disk_correlate(gtmodel, correlation_radius)
correlated_counts = disk_correlate(counts, correlation_radius)
correlated_model = disk_correlate(model, correlation_radius)

# Fermi significance
fermi_significance = np.nan_to_num(significance(correlated_counts, correlated_gtmodel, method='lima'))
# Gammapy significance
significance = np.nan_to_num(significance(correlated_counts, correlated_model, method='lima'))
# Ratio values
ratio = np.nan_to_num(correlated_model/correlated_gtmodel)

# Plotting
vmin, vmax = 0, 10
fig, axes = plt.subplots(nrows=1, ncols=3)
im = axes.flat[0].imshow(correlated_model,
                         interpolation='nearest',
                         origin="lower", vmin=vmin, vmax=vmax,
                         cmap=plt.get_cmap())
axes.flat[0].set_title('Gammapy Significance', fontsize=12)
im = axes.flat[1].imshow(correlated_gtmodel,
                         interpolation='nearest',
                         origin="lower", vmin=vmin, vmax=vmax,
                         cmap=plt.get_cmap())
axes.flat[1].set_title('Fermi Tools Significance', fontsize=12)
im = axes.flat[2].imshow(ratio,
                         interpolation='nearest',
                         origin="lower", vmin=vmin, vmax=vmax,
                         cmap=plt.get_cmap())
axes.flat[2].set_title('Ratio: \n Gammapy/Fermi Tools', fontsize=12)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.3, 0.025, 0.4])
fig.colorbar(im, cax=cbar_ax)
a = fig.get_axes()[0]
b = fig.get_axes()[1]
c = fig.get_axes()[2]
a.set_axis_off()
b.set_axis_off()
c.set_axis_off()
