from gammapy.image.catalog import catalog_image
from gammapy.image import make_empty_image
from gammapy.irf import EnergyDependentTablePSF
from gammapy.datasets import FermiGalacticCenter
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


psf_file = FermiGalacticCenter.filenames()['psf']

nxpix = 1000
nypix = 1000
binsz = 0.1

reference = make_empty_image(nxpix, nypix, binsz)
psf = EnergyDependentTablePSF.read(psf_file)

image = catalog_image(reference, psf, catalog='1FHL', source_type = 'point',
                  total_flux='True')
image_hdus = image.to_fits()

image_hdu = image_hdus[0]