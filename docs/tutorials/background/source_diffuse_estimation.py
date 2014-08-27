"""Estimate a diffuse emission model from Fermi LAT data.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from gammapy.datasets import FermiGalacticCenter
from gammapy.background import IterativeKernelBackgroundEstimator, GammaImages
from gammapy.irf import EnergyDependentTablePSF
from gammapy.image import make_empty_image, catalog_image, binary_disk
from gammapy.image.utils import cube_to_image, solid_angle

# *** PREPARATION ***

# Parameters

CORRELATION_RADIUS = 0.3
SIGNIFICANCE_THRESHOLD = 4
MASK_DILATION_RADIUS = 0.3 # pix

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

source_kernel = binary_disk(CORRELATION_RADIUS).astype(float)
source_kernel /= source_kernel.sum()

background_kernel = np.ones((5, 5))

# *** ITERATOR ***

ibe = IterativeKernelBackgroundEstimator(images=images,
                                         source_kernel=source_kernel,
                                         background_kernel=background_kernel,
                                         significance_threshold=SIGNIFICANCE_THRESHOLD,
                                         mask_dilation_radius=MASK_DILATION_RADIUS,
                                         save_intermediate_results=True
                                         )

n_iterations = 5

# *** RUN & PLOT ***

for iteration in range(n_iterations):
    ibe.run_iteration()
    mask_hdu = ibe.mask_image_hdu
    mask = mask_hdu.data

    plt.subplot(n_iterations, 2, 2 * iteration + 1)
    background_hdu = ibe.background_image_hdu
    data = background_hdu.data
    plt.imshow(data)
    plt.contour(mask, levels=[0], linewidths=2, colors='white')
    plt.axis('off')
    
    plt.subplot(n_iterations, 2, 2 * iteration + 2)
    significance_hdu = ibe.significance_image_hdu
    data = significance_hdu.data
    plt.imshow(data, vmin=-3, vmax=5)
    plt.contour(mask, levels=[0], linewidths=2, colors='white')
    plt.axis('off')

plt.tight_layout()
