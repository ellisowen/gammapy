# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" Make an image from a source catalog, or simulated catalog, e.g 1FHL 2FGL etc
"""
import numpy as np
from astropy.coordinates import Angle
from astropy.wcs import WCS
from astropy.units import Quantity
from astropy.table import Table
from ..datasets.load import fetch_fermi_extended_sources, fetch_fermi_catalog
from . import coordinates
from ..spectral_cube import GammaSpectralCube

__all__ = ['catalog_image', 'catalog_table']


def _extended_image(catalog, reference_cube):
    """Reprojects and adds extended source images to a larger survey image.
    """
    # Note that the first extended source fits file is unreadable...
    hdu_list = fetch_fermi_extended_sources(catalog)[1:]
    for source in hdu_list:
        source_wcs = WCS(source.header)
        source_spec_cube = GammaSpectralCube(data=Quantity(np.array([source.data]), ''),
                                                 wcs=source_wcs, energy=energy)
        new_source_cube = source_spec_cube.reproject_to(reference_cube)
        # TODO: Fix this hack
        reference_cube.data = reference_cube.data + np.nan_to_num(new_source_cube.data * 1e-30)
    return reference_cube.data[0]


def _source_image(catalog, reference_cube, sim_table=None, total_flux=True):
    """Adds point sources to a larger survey image.
    """
    new_image = np.zeros_like(reference_cube.data, dtype=np.float64)
    if sim_table == None:
        source_table = catalog_table(catalog, ebands=False)
    else:
        source_table = sim_table
    energies = source_table.meta['Energy Bins']
    wcs_reference = reference_cube.wcs
    footprint = wcs_reference.calc_footprint()
    glon_max, glon_min = footprint[0][0], footprint[2][0] - 360
    glat_min, glat_max = footprint[0][1], footprint[1][1]
    for source in np.arange(len(source_table['flux'])):
        lon = source_table['GLON'][source]
        if lon >= 180:
            lon = lon - 360
        if (glon_min < lon) & (lon < glon_max):            
            lat = source_table['GLAT'][source]
            if (glat_min < lat) & (lat < glat_max):
                flux = source_table['flux'][source]
                wcs = reference_cube.wcs
                x, y = wcs.wcs_world2pix(lon, lat, 0)
                xi, yi = x.astype(int), y.astype(int)
                new_image[yi, xi] = new_image[yi, xi] + flux
    if total_flux == True:
        factor = source_table['flux'].sum() / new_image.sum()
        return new_image * factor, energies
    else:
        return new_image, energies


def catalog_image(reference, psf, catalog='1FHL', source_type='point',
                  total_flux=False, sim_table=None):  
    """Creates an image from a simulated catalog, or from 1FHL or 2FGL sources.

    Parameters
    ----------
    reference : `~fits.ImageHDU`
        Reference Image HDU. The output takes the shape and resolution of this.
    psf : `~gammapy.irf.EnergyDependentTablePSF`
        Energy dependent Table PSF object for image convolution.
    catalog : str
        Flag which source catalog is to be used to create the image.
        Options: {'1FHL', '2FGL', 'simulation'}
        If 'simulation' is used, sim_table must also be provided.
    source_type : str
        Specify whether point or extended sources should be included, or both.
        Options: {'point', 'extended', 'all'}
        TODO: Currently only 'point' is implemented.
    total_flux : bool
        Specify whether to conserve total flux.
    sim_table : `~astropy.table.Table`
        Table of simulated point sources. Only required if catalog = 'simulation'

    Returns
    -------
    out_cube : `~gammapy.spectral_cube.GammaSpectralCube`
        2D Spectral cube containing the image.

    Notes
    -----
    This is currently only implemented for a single energy band.
    """
    from scipy.ndimage import convolve
    lons, lats = coordinates(reference)
    wcs = WCS(reference.header)
    # Uses dummy energy for now to construct spectral cube
    # TODO : Fix this hack
    reference_cube = GammaSpectralCube(data=Quantity(np.array(reference.data), ''),
                                          wcs=wcs, energy=Quantity([0, 1], 'GeV'))
    if source_type == 'extended':
        raise NotImplementedError
        # TODO: Currently fluxes are not correct for extended sources.
        new_image = _extended_image(catalog, reference_cube)
    elif source_type == 'point':
        new_image, energy = _source_image(catalog, reference_cube, sim_table, total_flux)
    elif source_type == 'all':
        raise NotImplementedError
        # TODO: Currently Extended Sources do not work
        extended = _extended_image(catalog, reference_cube)
        point_source = _source_image(catalog, reference_cube, total_flux=True)[0]
        new_image = extended + point_source
    else:
        raise ValueError
    total_point_image = GammaSpectralCube(data=new_image, wcs=wcs, energy=energy)
    convolved_cube = new_image.copy()
    psf = psf.table_psf_in_energy_band(Quantity([np.min(energy).value,
                                        np.max(energy).value], energy.unit))
    resolution = abs(reference.header['CDELT1'])
    kernel_array = psf.kernel(pixel_size=Angle(resolution, 'deg'),
                              offset_max=Angle(5, 'deg'), normalize=True)
    convolved_cube = convolve(new_image, kernel_array, mode='constant')
    out_cube = GammaSpectralCube(data=convolved_cube, wcs=total_point_image.wcs,
                                     energy=energy)
    return out_cube
 

def catalog_table(catalog, ebands=False):
    """ Creates catalog table from published source catalog.

    Parameters
    ----------
    catalog : str
        Catalog to load.
        Options: {'1FHL', '2FGL'}
    ebands : bool
        Whether to return catalog in energy bands.

    Returns
    -------
    table : `~astropy.table.Table`
        Point source catalog table.
    """
    data = []
    cat_table = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')
    for source in np.arange(len(cat_table)):
        glon = cat_table['GLON'][source]
        glat = cat_table['GLAT'][source]
        # Different from here between each catalog because of different catalog header names
        if catalog == ('1FHL' or 'simulation'):
            energy = Quantity([10, 30, 100, 500], 'GeV')
            if ebands == False:
                flux_bol = cat_table['Flux'][source]
                row = dict(Source_Type='PointSource',
                           GLON=glon, GLAT=glat, flux=flux_bol)
            else:
                Flux_10_30 = cat_table['Flux10_30GeV'][source]
                Flux_30_100 = cat_table['Flux30_100GeV'][source]
                Flux_100_500 = cat_table['Flux100_500GeV'][source]
                row = dict(Source_Type='PointSource',
                           GLON=glon, GLAT=glat, Flux10_30=Flux10_30,
                           Flux30_100=Flux30_100, Flux100_500=Flux100_500)
        elif catalog == '2FGL':
            energy = Quantity([30, 100, 300, 1000, 3000, 10000, 100000], 'GeV') 
            if ebands == False:
                flux_bol = cat_table['Flux_Density'][source]
                row = dict(Source_Type='PointSource',
                           GLON=glon, GLAT=glat, flux=flux_bol)
            else:
                Flux_30_100 = cat_table['Flux30_100'][source]
                Flux_100_300 = cat_table['Flux100_300'][source]
                Flux_300_1000 = cat_table['Flux300_1000'][source]
                Flux_1000_3000 = cat_table['Flux1000_3000'][source]
                Flux_3000_10000 = cat_table['Flux3000_10000'][source]
                Flux_10000_100000 = cat_table['Flux10000_100000'][source]
                row = dict(Source_Type='PointSource', Source_Name=source_name,
                           GLON=glon, GLAT=glat, Flux_30_100=Flux_30_100,
                           Flux_100_300=Flux_100_300, Flux_300_1000=Flux_300_1000,
                           Flux_1000_3000=Flux_1000_3000, Flux_3000_10000=Flux_3000_10000,
                           Flux_10000_100000=Flux_10000_100000)
        data.append(row)
    table = Table(data)
    table.meta['Energy Bins'] = energy
    return table
