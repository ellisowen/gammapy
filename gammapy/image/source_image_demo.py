from gammapy.image.catalog import catalog_image
from gammapy.image import make_empty_image
from gammapy.irf import EnergyDependentTablePSF
from gammapy.datasets import FermiGalacticCenter
from astropy.io import fits

psf_file = FermiGalacticCenter.filenames()['psf']

nxpix = 1000
nypix = 1000
binsz = 0.1

reference = make_empty_image(nxpix, nypix, binsz)
psf = EnergyDependentTablePSF.read(psf_file)

image = catalog_image(reference, psf, catalog='1FHL', source_type = 'point',
                  total_flux='True')

a = image.wcs
image_hdu = fits.ImageHDU(data = image.data, header = a.to_header())

import IPython; IPython.embed()