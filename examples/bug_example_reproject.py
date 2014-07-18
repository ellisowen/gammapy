from FITS_tools.cube_regrid import regrid_cube_hdu
from gammapy.datasets import FermiGalacticCenter
from astropy.io import fits
from gammapy.spectral_cube.core import compute_npred_cube

diffuse = fits.open(FermiGalacticCenter.filenames()['diffuse_model'])[0]
exposure = fits.open(FermiGalacticCenter.filenames()['exposure_cube'])[0]
import IPython; IPython.embed()
compute_npred_cube(diffuse, exposure, desired_energy=None,
                       convolve='Yes', max_convolution_offset=5)
#regrid_cube_hdu(exposure, diffuse.header, smooth=False)