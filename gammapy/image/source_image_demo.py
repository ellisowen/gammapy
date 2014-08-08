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



import numpy as np
import matplotlib.pyplot as plt



# Spatial distribution using Lorimer (2006) model
from gammapy.astro import population

rho_sun = 3 # source density at the sun (sources kpc^-1)
n_sources = int(5e2) # TODO: I guess this should be computed from rho_sun?
rad_dis = 'L06'
vel_dis = 'F06B'
spiralarms = True
table = population.make_cat_gal(nsources=n_sources, rad_dis=rad_dis, vel_dis=vel_dis,
                                max_age=1e6, spiralarms=spiralarms)



# This is done in `make_cat_gal`:
#from gammapy.astro.population import L06
#from gammapy.utils.distributions import draw, GeneralRandom
#galactic_radius = np.linspace(0, 50, 1000)
#density = L06(galactic_radius)
#galactic_angle = np.random.uniform(0, 360, n_sources)
#plt.plot(galactic_radius, density)



luminosity_min  = 4e34 # Minimum source luminosity (ph s^-1)
luminosity_max = 4e37 # Maximum source luminosity (ph s^-1)
luminosity_index = 1.5 # Luminosity function differential power-law index

from gammapy.utils.random import sample_powerlaw

luminosity = sample_powerlaw(luminosity_min, luminosity_max, luminosity_index, n_sources)
table['luminosity'] = luminosity



from gammapy.astro.population import add_observed_parameters
table = add_observed_parameters(table)

table.write('simulated_galaxy_3.fits', overwrite=True)

plt.hist(np.log10(table['flux']), bins=25, log=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlabel(r'Flux, photons cm$^{-2}$s$^{-1}$')
plt.ylabel(r'Source Number Counts')
plt.savefig('galaxy_population_plot3.pdf')