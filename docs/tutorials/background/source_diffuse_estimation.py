"""Estimate a diffuse emission model from Fermi LAT data.
"""
import numpy as np
from astropy.io import fits
from gammapy.datasets import FermiGalacticCenter
from gammapy.background.kernel import GammaImages, IterativeBackgroundEstimator
from gammapy.irf import EnergyDependentTablePSF
from gammapy.image import make_empty_image, catalog_image, binary_disk
from gammapy.image.utils import cube_to_image, solid_angle
from aplpy import FITSFigure

# *** PREPARATION ***

# Parameters
TOTAL_COUNTS = 1e6
SOURCE_FRACTION = 0.2

CORRELATION_RADIUS = 3
SIGNIFICANCE_THRESHOLD = 5
MASK_DILATION_RADIUS = 3 # pix
NUMBER_OF_ITERATIONS = 4

# Derived parameters
DIFFUSE_FRACTION = 1. - SOURCE_FRACTION

psf_file = FermiGalacticCenter.filenames()['psf']
psf = EnergyDependentTablePSF.read(psf_file)

# Load/create example model images
filename = FermiGalacticCenter.filenames()['diffuse_model']
# May only put in a 1D image
diffuse_image_hdu = cube_to_image(fits.open(filename)[0], 0)
solid_angle_image = solid_angle(diffuse_image_hdu)
diffuse_image_true = diffuse_image_hdu.data * solid_angle_image.value
reference = make_empty_image(nxpix=61, nypix=21, binsz=0.5)
sources = catalog_image(reference, psf, catalog='1FHL',
                        source_type='point', total_flux='True')
source_image_true = sources.data  
total_image_true = source_image_true + diffuse_image_true

# This is a flux image. Need to create a counts image.
flux_data = total_image_true
exposure_filename = FermiGalacticCenter.filenames()['exposure_cube']
# Assume uniform exposure as an approximation (saves reprojecting the cube)
exposure_value = fits.open(exposure_filename)[0].data.mean()
counts_data = flux_data * exposure_value

# *** LOADING INPUT ***

# Counts must be provided as an ImageHDU
counts = fits.ImageHDU(data=counts_data, header=reference.header)
# Start with flat background estimate
# Background must be provided as an ImageHDU
background_data=np.ones_like(counts_data, dtype=float)
background = fits.ImageHDU(data=background_data, header=reference.header)
images = GammaImages(counts=counts, background=background)

source_kernel = np.ones((3, 3))
#source_kernel = binary_disk(CORRELATION_RADIUS).astype(float)
source_kernel /= source_kernel.sum()

background_kernel = np.ones((5, 5))

# *** ITERATOR ***

ibe = IterativeBackgroundEstimator(images=images,
                                   source_kernel=source_kernel,
                                   background_kernel=background_kernel,
                                   significance_threshold=SIGNIFICANCE_THRESHOLD,
                                   mask_dilation_radius=MASK_DILATION_RADIUS
                                   )

mask, background = ibe.run()

fig = FITSFigure(background)
fig.show_colorscale(stretch='linear', interpolation='none')
fig.add_colorbar()