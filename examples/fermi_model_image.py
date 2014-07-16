"""Script to generate fermi model image from counts cube, exposure cube and psf
"""

from gammapy.spectral_cube import GammaSpectralCube
from astropy.units import Quantity
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
dir = "/home/ellis/vela_analysis/"
filename = dir + "gll_iem_v05.fits"
energy = Quantity([1, 10000000], 'MeV')
wcs = WCS(filename)
#read both in
#resample to same energy
#multiply exposure and flux cubes to get npred (predicted counts)
#USE small region - vela from tuesday(!)
from gammapy.image.utils import coordinates
from gammapy.image.utils import cube_to_image
image = fits.ImageHDU(data = fits.open(filename)[0].data, header = fits.open(filename)[0].header)
image = cube_to_image(image)
lat, lon = coordinates(image, world=True, radians=True)
lon = Quantity(lon, 'rad')
lat = Quantity(lat, 'rad')
cube = GammaSpectralCube(filename, wcs, energy)
new_energy = Quantity(10, 'MeV')
cube.flux(lon, lat, new_energy)

