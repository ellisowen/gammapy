"""Runs commands to produce convolved predicted counts map in current directory"""
from gammapy.spectral_cube import compute_npred_cube, GammaSpectralCube, convolve_npred_cube
from astropy.units import Quantity
background_model = GammaSpectralCube.read('/home/eowen/analyses/galactic_high_energy/gll_iem_v05.fits')
exposure_cube = GammaSpectralCube.read('exposure_cube_lowres.fits')
print('Reprojecting Background Cube')
repro_bg_cube = background_model.reproject_to(exposure_cube)

energies = Quantity([10, 30, 100, 500], 'GeV')
print('Generating Predicted Counts Cube')
npred_cube = compute_npred_cube(repro_bg_cube, exposure_cube, energies)
npred_cube.write_to_fits('npred_cube.fits')
print('Convolving npred Cube')
convolved_npred_cube = convolve_npred_cube(npred_cube, 3, 0.1)
convolved_npred_cube.write_to_fits('convolved_npred_cube.fits')
print('Done')
