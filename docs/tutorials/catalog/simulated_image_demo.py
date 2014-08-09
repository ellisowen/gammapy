"""Simulates a galaxy of point sources and produces an image.
"""
import matplotlib.pyplot as plt
from gammapy.astro import population
from gammapy.astro.population import add_observed_parameters
from gammapy.datasets import FermiGalacticCenter
from gammapy.image import make_empty_image
from gammapy.image.catalog import catalog_image
from gammapy.irf import EnergyDependentTablePSF
from gammapy.utils.random import sample_powerlaw

# Define image size
nxpix = 1000
nypix = 1000
binsz = 0.1

reference = make_empty_image(nxpix, nypix, binsz)

psf_file = FermiGalacticCenter.filenames()['psf']
psf = EnergyDependentTablePSF.read(psf_file)

# Simulation Parameters

# source density at the sun (sources kpc^-1)
rho_sun = 3 
# number of sources
n_sources = int(5e2)
# Spatial distribution using Lorimer (2006) model
rad_dis = 'L06'
# Velocity dispersion
vel_dis = 'F06B'
# Includes spiral arms
spiralarms = True
# Creates table
table = population.make_cat_gal(nsources=n_sources, rad_dis=rad_dis, vel_dis=vel_dis,
                                max_age=1e6, spiralarms=spiralarms)

# Minimum source luminosity (ph s^-1)
luminosity_min  = 4e34
# Maximum source luminosity (ph s^-1)
luminosity_max = 4e37
# Luminosity function differential power-law index
luminosity_index = 1.5

# Assigns luminosities to sources
luminosity = sample_powerlaw(luminosity_min, luminosity_max, luminosity_index, n_sources)
table['luminosity'] = luminosity

# Adds parameters to table: distance, glon, glat, flux, angular_extension
table = add_observed_parameters(table)

# Create image
image = catalog_image(reference, psf, catalog='simulation', source_type = 'point',
                  total_flux=True, sim_table=table)

# Plot
image_hdus = image.to_fits()
image = image_hdus[0].data
plt.imshow(image)