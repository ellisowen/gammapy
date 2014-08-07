"""Runs commands to produce convolved predicted flux map of sources based on input catalog
"""
from gammapy.image.catalog import catalog_image
from gammapy.irf import EnergyDependentTablePSF
from gammapy.image.utils import make_empty_image
import numpy as np
from gammapy.datasets import FermiGalacticCenter
from astropy.units import Quantity
import matplotlib.pyplot as plt

resolution=0.1
center=[0, 0]
lat_range=[0, 180]
lon_range=[0, 360]

lat_pix = (lat_range[1]-lat_range[0])/resolution
lon_pix = (lon_range[1]-lon_range[0])/resolution
reference = make_empty_image(lon_pix, lat_pix, resolution, xref=center[0], yref=center[1], fill=0,
                             proj='CAR', coordsys='GAL', xrefpix=np.floor(0.5*lon_pix)+0.5,
                             yrefpix=np.floor(0.5*lat_pix)+0.5, dtype='float32')

psf = EnergyDependentTablePSF.read(FermiGalacticCenter.filenames()['psf'])

image = catalog_image(reference, psf, catalog='1FHL', source_type = 'ExtendedSource',
                  total_flux='False', filename='1fhl_fermi_psf.fits',
                  energy = Quantity([10, 500], 'GeV'))
import IPython; IPython.embed()