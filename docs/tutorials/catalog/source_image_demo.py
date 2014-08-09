"""Produces an image from 1FHL catalog point sources.
"""
import matplotlib.pyplot as plt
from gammapy.datasets import FermiGalacticCenter
from gammapy.image import make_empty_image
from gammapy.image.catalog import catalog_image
from gammapy.irf import EnergyDependentTablePSF

# Define image size
nxpix = 1000
nypix = 1000
binsz = 0.1

reference = make_empty_image(nxpix, nypix, binsz)
psf = EnergyDependentTablePSF.read(psf_file)

# Create image
image = catalog_image(reference, psf, catalog='1FHL', source_type = 'point',
                  total_flux='True')
image_hdus = image.to_fits()

# Plot
image = image_hdus[0].data
plt.imshow(image)