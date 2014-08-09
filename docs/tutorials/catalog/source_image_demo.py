"""Produces an image from 1FHL catalog point sources.
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from gammapy.datasets import FermiGalacticCenter
from gammapy.image import make_empty_image
from gammapy.image.catalog import catalog_image
from gammapy.irf import EnergyDependentTablePSF

# Define image size
nxpix = 100
nypix = 100
binsz = 1

reference = make_empty_image(nxpix, nypix, binsz)
psf_file = FermiGalacticCenter.filenames()['psf']
psf = EnergyDependentTablePSF.read(psf_file)

# Create image
image = catalog_image(reference, psf, catalog='1FHL', source_type='point',
                  total_flux='True')

# Plot
image_hdus = image.to_fits()
image = image_hdus[0].data
fig = plt.imshow(image, cmap = cm.Greys_r)
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
