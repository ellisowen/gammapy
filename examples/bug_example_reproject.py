from FITS_tools.cube_regrid import regrid_cube_hdu
from gammapy.datasets import FermiGalacticCenter
from astropy.io import fits

diffuse = fits.open(FermiGalacticCenter.filenames()['diffuse_model'])[0]
exposure = fits.open(FermiGalacticCenter.filenames()['exposure_cube'])[0]

regrid_cube_hdu(exposure, diffuse.header, smooth=False)