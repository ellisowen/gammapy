"""Script to generate fermi model image from counts cube, exposure cube and psf
"""

from gammapy.spectral_cube import GammaSpectralCube
from astropy.units import Quantity
from gammapy.image import coordinates
from astropy.coordinates import Angle
from gammapy.datasets import FermiGalacticCenter
from gammapy.irf import EnergyDependentTablePSF
from scipy.ndimage import convolve


def correlate_fermi_psf(image, max_offset, resolution=0.1, energy = 'None', energy_band=[10, 500]):
    # Parameters
    filename = FermiGalacticCenter.filenames()['psf']
    pixel_size = Angle(resolution, 'deg')
    offset_max = Angle(max_offset, 'deg')
    if energy == 'None':
        energy_band = Quantity(energy_band, 'GeV')
        fermi_psf = EnergyDependentTablePSF.read(filename)
        psf = fermi_psf.table_psf_in_energy_band(energy_band=energy_band, spectral_index=2.5)
    else:
        energy = Quantity(energy, 'GeV')
        fermi_psf = EnergyDependentTablePSF.read(filename)
        psf = fermi_psf.table_psf_at_energy(energy=energy)
    psf.normalize()
    kernel = psf.kernel(pixel_size=pixel_size, offset_max=offset_max)
    kernel_image = kernel.value/kernel.value.sum()
    
    # TODO: Write unit test (this will be useful):
    
    #kernel_image_integral = kernel_image.sum() * pixel_size.to('radian').value ** 2
    #print('Kernel image integral: {0}'.format(kernel_image_integral))
    #print('shape: {0}'.format(kernel_image.shape))
    return convolve(image, kernel_image, mode='constant')

def interp_cube(filename, energy_band, new_energy):
    cube = GammaSpectralCube.read(filename)
    int_flux_image = cube.integral_flux_image(energy_band)
    lat, lon = coordinates(int_flux_image, world=True, radians=True)
    

    lat, lon = coordinates(int_flux_image, world=True, radians=True)

    lat = Quantity(lat, 'rad')
    lon = Quantity(lon, 'rad')

    array = cube.flux(lat, lon, new_energy)

    return array.reshape(lat.shape)

def exposure_map_from_cube(filename, energy):
    from astropy.io import fits
    import numpy as np
    cube = fits.open(filename)
    
    max_energy = np.array(Quantity([cube[1].data['Energy']], 'MeV')).max()
    import IPython; IPython.embed()
    if energy >= max_energy:
        exposure = cube[0].data[30]
    
    else:
        raise NotImplementedError
    return exposure


if __name__=='__main__':

    # Integral flux cube script

    filename = "/home/eowen/software/python/gammapy/gammapy/datasets/data/fermi/gll_iem_v02_cutout.fits"
    cube = GammaSpectralCube.read(filename)
    energy_band = Quantity([10, 100], 'MeV')
    int_flux_image = cube.integral_flux_image(energy_band)
    print(int_flux_image.data)

    # Predicted counts script
    
    #filename = "/home/eowen/software/python/gammapy/gammapy/datasets/data/fermi/gll_iem_v02_cutout.fits"
    #energy_band = Quantity([10, 100], 'MeV')
    #new_energy = Quantity(1, 'TeV')
    #test_image = interp_cube(filename, energy_band, new_energy)
    
    source_filename = "/home/eowen/software/python/gammapy/gammapy/datasets/data/fermi/gll_iem_v02_cutout.fits"
    background_filename = "/home/eowen/software/python/gammapy/gammapy/datasets/data/fermi/gll_iem_v02_cutout.fits"
    source_energy_band = Quantity([10, 100], 'MeV')
    background_energy_band = Quantity([10, 100], 'MeV')
    new_energy = Quantity(1, 'TeV')
    
    point_source_flux = interp_cube(source_filename, source_energy_band, new_energy)
    background_model_flux = interp_cube(background_filename, background_energy_band, new_energy)

    total_flux = point_source_flux + background_model_flux
    flux_psf_conv = correlate_fermi_psf(total_flux, max_offset=3, resolution=0.1, energy = new_energy, energy_band=energy_band)
    
    exposure_filename = "/home/eowen/software/python/gammapy/gammapy/datasets/data/fermi/fermi_exposure_cutout.fits"
    exposure_map = exposure_map_from_cube(exposure_filename, new_energy)
    import IPython; IPython.embed()
    counts_map = exposure_map * flux_psf_conv