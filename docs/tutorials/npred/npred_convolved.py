"""Runs commands to produce convolved predicted counts map in current directory
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.units import Quantity
from astropy.io import fits
from astropy.wcs import WCS
from gammapy.spectral_cube import (compute_npred_cube, GammaSpectralCube,
                                   convolve_npred_cube)
from gammapy.datasets import FermiVelaRegion
from gammapy.irf import EnergyDependentTablePSF

# Reads in data
background_file = FermiVelaRegion.filenames()['diffuse_model']
exposure_file = FermiVelaRegion.filenames()['exposure_cube']
counts_file = FermiVelaRegion.filenames()['counts_cube']
background_model = GammaSpectralCube.read(background_file)
exposure_cube = GammaSpectralCube.read(exposure_file)

# Reproject background cube
repro_bg_cube = background_model.reproject_to(exposure_cube)

# Define energy band required for output
energies = Quantity([10, 500], 'GeV')

# Compute the predicted counts cube
npred_cube = compute_npred_cube(repro_bg_cube, exposure_cube, energies)

# Convolve with Energy-dependent Fermi LAT PSF
psf = EnergyDependentTablePSF.read(FermiVelaRegion.filenames()['psf'])
convolved_npred_cube = convolve_npred_cube(npred_cube, psf, 3, 1)

# Counts data
counts_data = fits.open(counts_file)[0].data
counts_wcs = WCS(fits.open(counts_file)[0].header)
counts_cube = GammaSpectralCube(data=Quantity(counts_data, ''), wcs=counts_wcs,
                                energy=energies)
counts_cube = counts_cube.reproject_to(npred_cube)

counts = counts_cube.data[0]
model = convolved_npred_cube.data[0]

# Load Fermi tools gtmodel result
gtmodel = fits.open(FermiVelaRegion.filenames()['background_image'])[0].data.astype(float)

# Ratio for the two background images
ratio = np.nan_to_num(model / gtmodel)

# Plotting
vmin, vmax = 0, 1

fig, axes = plt.subplots(nrows=1, ncols=3)

results = [model, gtmodel, ratio]
titles = ['Gammapy Background', 'Fermi Tools Background', 'Ratio: \n Gammapy/Fermi Tools']

for i in np.arange(3):
    im = axes.flat[i].imshow(results[i],
                             interpolation='nearest',
                             origin="lower", vmin=vmin, vmax=vmax,
                             cmap=plt.get_cmap())

    axes.flat[i].set_title(titles[i], fontsize=12)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.3, 0.025, 0.4])
fig.colorbar(im, cax=cbar_ax)
a = fig.get_axes()[0]
b = fig.get_axes()[1]
c = fig.get_axes()[2]
a.set_axis_off()
b.set_axis_off()
c.set_axis_off()
