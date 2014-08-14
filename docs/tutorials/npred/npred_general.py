"""Runs commands to produce convolved predicted counts map in current directory
"""

import numpy as np
from astropy.coordinates import Angle
from astropy.units import Quantity
from astropy.io import fits
from astropy.wcs import WCS
from gammapy.spectral_cube import (compute_npred_cube, GammaSpectralCube,
                                   convolve_npred_cube)
from gammapy.datasets import FermiVelaRegion
from gammapy.irf import EnergyDependentTablePSF

__all__ = ['prepare_images']

def prepare_images():

    # Reads in data
    background_file = "/home/ellis/software/python/gammapy-extra/datasets/vela_region/gll_iem_v05_rev1_cutout.fits"#FermiVelaRegion.filenames()['diffuse_model']
    exposure_file = "/home/ellis/software/python/gammapy-extra/datasets/vela_region/exposure_vela.fits"#FermiVelaRegion.filenames()['exposure_cube']
    counts_file = "/home/ellis/software/python/gammapy-extra/datasets/vela_region/counts_vela.fits"#FermiVelaRegion.filenames()['counts_cube']
    background_model = GammaSpectralCube.read(background_file)
    exposure_cube = GammaSpectralCube.read(exposure_file)
    # Add correct units
    exposure_cube.data = Quantity(exposure_cube.data, '1/(cm2 deg2 s GeV)')
    # Reproject background cube
    repro_bg_cube = background_model.reproject_to(exposure_cube)
    # Define energy band required for output
    energies = Quantity([10, 500], 'GeV')

    # Compute the predicted counts cube
    npred_cube = compute_npred_cube(repro_bg_cube, exposure_cube, energies)

    # Convolve with Energy-dependent Fermi LAT PSF
    psf = EnergyDependentTablePSF.read("/home/ellis/software/python/gammapy-extra/datasets/vela_region/psf_vela.fits")
    convolved_npred_cube = convolve_npred_cube(npred_cube, psf,
                                               Angle(3, 'deg'), Angle(1, 'deg'))

    # Counts data
    counts_data = fits.open(counts_file)[0].data
    counts_wcs = WCS(fits.open(counts_file)[0].header)
    counts_cube = GammaSpectralCube(data=Quantity(counts_data, ''), wcs=counts_wcs,
                                    energy=energies)
    counts_cube = counts_cube.reproject_to(npred_cube)

    counts = counts_cube.data[0]
    model = convolved_npred_cube.data[0]

    # Load Fermi tools gtmodel result
    gtmodel = fits.open("/home/ellis/software/python/gammapy-extra/datasets/vela_region/gtmodel_nosrc.fits")[0].data.astype(float)

    # Ratio for the two background images
    ratio = np.nan_to_num(model / gtmodel)

    return model, gtmodel, ratio, counts
