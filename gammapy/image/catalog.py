# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" make an image from a source catalog, or simulated catalog, e.g 1FHL 2FGL etc
"""
__all__ = ['catalog_image', 'catalog_table']

from astropy.io import fits
from astropy.modeling.models import Gaussian2D, Disk2D
from gammapy.image.utils import make_header, disk_correlate
from gammapy.datasets.load import fetch_fermi_extended_sources, fetch_fermi_catalog
from gammapy.datasets import FermiGalacticCenter
from gammapy.spectral_cube import GammaSpectralCube
from gammapy.image import coordinates
from astropy.wcs import WCS
from astropy.units import Quantity
from astropy.table import Table
from gammapy.irf import EnergyDependentTablePSF
import numpy as np
from astropy.coordinates import Angle


def _extended_image(catalog, energy, reference_cube):
    # Note that the first extended source fits file is unreadable...
    hdu_list = fetch_fermi_extended_sources(catalog)[1:]
    for source in hdu_list:
        source_wcs = WCS(source.header)
        source_spec_cube = GammaSpectralCube(data = Quantity(np.array([source.data]), ''),
                                                 wcs=source_wcs, energy=energy)
        new_source_cube = source_spec_cube.reproject_to(reference_cube)
        # TODO: Fix this hack
        reference_cube.data = reference_cube.data + np.nan_to_num(new_source_cube.data * 1e-30)
    return reference_cube.data[0]


def _source_image(catalog, energy, reference_cube):
    new_image = np.zeros_like(reference_cube.data, dtype=np.float64)
    source_table = catalog_table(catalog, ebands='No')
    for source in np.arange(len(source_table['Flux'])):
        lon = source_table['GLON'][source]
        lat = source_table['GLAT'][source]
        flux = source_table['Flux'][source]
        wcs = reference_cube.wcs
        x, y = wcs.wcs_world2pix(lon, lat, 0)
        xi, yi = x.astype(int), y.astype(int)
        new_image[0][yi, xi] = new_image[0][yi, xi] + flux
    return new_image[0]


def catalog_image(reference, psf, catalog='1FHL', source_type = 'All',
                  total_flux='False', filename='1fhl_fermi_psf.fits',
                  energy = Quantity([10, 500], 'GeV')):  
    """TODO
    """
    lons, lats = coordinates(reference)
    wcs = WCS(reference.header)
    reference_cube = GammaSpectralCube(data = Quantity(np.array([reference.data]), ''),
                                          wcs = wcs, energy = energy)
    if source_type == 'ExtendedSource':
        new_image = _extended_image(catalog, energy, reference_cube)
    elif source_type == 'PointSource':
        new_image = _source_image(catalog, energy, reference_cube)
    elif source_type == 'All':
        new_image = _extended_image(catalog, energy, reference_cube) + _source_image(catalog, energy, reference_cube)
    else:
        raise ValueError
    total_point_image = GammaSpectralCube(data = new_image, wcs = wcs, energy = energy)
    from scipy.ndimage import convolve
    psf_object = psf
    convolved_cube = new_image.copy()
    psf = psf_object.table_psf_in_energy_band(Quantity([energy[0].value,
                                                            energy[1].value],
                                                           energy.unit))
    resolution = abs(reference.header['CDELT1'])
    kernel_array = psf.kernel(Angle(resolution, 'deg'), Angle(5, 'deg'))
    kernel_image = kernel_array / kernel_array.sum()
    convolved_cube = convolve(new_image, kernel_image,
                              mode='constant')
    out_cube = GammaSpectralCube(data=convolved_cube, wcs=total_point_image.wcs,
                                     energy=energy)
    
    if total_flux == 'True':
        factor = source_table['Flux'].sum()
        out_cube.data = ((out_cube.data / out_cube.data.sum()) * factor)
    else:
        out_cube.data = out_cube.data
    return out_cube
 

def catalog_table(catalog, ebands='No'):
    """ TODO
    """
    from scipy.stats import gmean
    data = []
    cat_table = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')
    for source in np.arange(len(cat_table)):
        glon = cat_table['GLON'][source]
        glat = cat_table['GLAT'][source]
        source_name = cat_table['Source_Name'][source]
        # Different from here between each catalog because of different catalog header names
        if catalog == '1FHL':
            energy = Quantity([10, 30, 100, 500], 'GeV')
            if ebands == 'No':
                flux_bol = cat_table['Flux'][source]
                row = dict(Source_Type='PointSource', Source_Name=source_name,
                           GLON=glon, GLAT=glat, Flux=flux_bol)
            else:
                Flux_10_30 = cat_table['Flux10_30GeV'][source]
                Flux_30_100 = cat_table['Flux30_100GeV'][source]
                Flux_100_500 = cat_table['Flux100_500GeV'][source]
                row = dict(Source_Type='PointSource', Source_Name=source_name,
                           GLON=glon, GLAT=glat, Flux10_30=Flux10_30,
                           Flux30_100=Flux30_100, Flux100_500=Flux100_500)
        elif catalog == '2FGL':
            energy = Quantity([30, 100, 300, 1000, 3000, 10000, 100000], 'GeV') 
            if ebands == 'No':
                flux_bol = cat_table['Flux_Density'][source]
                row = dict(Source_Type='PointSource', Source_Name=source_name,
                           GLON=glon, GLAT=glat, Flux=flux_bol)
            else:
                Flux_30_100 = cat_table['Flux30_100'][source]
                Flux_100_300  = cat_table['Flux100_300'][source]
                Flux_300_1000  = cat_table['Flux300_1000'][source]
                Flux_1000_3000  = cat_table['Flux1000_3000'][source]
                Flux_3000_10000  = cat_table['Flux3000_10000'][source]
                Flux_10000_100000  = cat_table['Flux10000_100000'][source]
                row = dict(Source_Type='PointSource', Source_Name=source_name,
                           GLON=glon, GLAT=glat, Flux_30_100=Flux_30_100,
                           Flux_100_300=Flux_100_300, Flux_300_1000=Flux_300_1000,
                           Flux_1000_3000=Flux_1000_3000, Flux_3000_10000=Flux_3000_10000,
                           Flux_10000_100000=Flux_10000_100000)
        data.append(row)
    table = Table(data)
    table.meta['Energy Bins'] = energy
    return table
