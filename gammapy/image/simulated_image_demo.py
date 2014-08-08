from gammapy.image.catalog import catalog_image
from gammapy.image import make_empty_image
from gammapy.irf import EnergyDependentTablePSF
from gammapy.datasets import FermiGalacticCenter
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
# Spatial distribution using Lorimer (2006) model
from gammapy.utils.random import sample_powerlaw
from gammapy.astro import population
from gammapy.astro.population import add_observed_parameters
psf_file = FermiGalacticCenter.filenames()['psf']

nxpix = 1000
nypix = 1000
binsz = 0.1

reference = make_empty_image(nxpix, nypix, binsz)
psf = EnergyDependentTablePSF.read(psf_file)

rho_sun = 3 # source density at the sun (sources kpc^-1)
n_sources = int(5e2) # TODO: I guess this should be computed from rho_sun?
rad_dis = 'L06'
vel_dis = 'F06B'
spiralarms = True
table = population.make_cat_gal(nsources=n_sources, rad_dis=rad_dis, vel_dis=vel_dis,
                                max_age=1e6, spiralarms=spiralarms)
# TODO: Up to here...
luminosity_min  = 4e34 # Minimum source luminosity (ph s^-1)
luminosity_max = 4e37 # Maximum source luminosity (ph s^-1)
luminosity_index = 1.5 # Luminosity function differential power-law index

luminosity = sample_powerlaw(luminosity_min, luminosity_max, luminosity_index, n_sources)
table['luminosity'] = luminosity

table = add_observed_parameters(table)
import IPython; IPython.embed()
# pass table into function

image = catalog_image(reference, psf, catalog='simulation', source_type = 'point',
                  total_flux=True, sim_hdu_list=table)