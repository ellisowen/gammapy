# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Gamma-ray spectral cube.

TODO: split `GammaSpectralCube` into a base class ``SpectralCube`` and a few sub-classes:

* ``GammaSpectralCube`` to represent functions evaluated at grid points (diffuse model format ... what is there now).
* ``ExposureCube`` should also be supported (same semantics, but different units / methods as ``GammaSpectralCube`` (``gtexpcube`` format)
* ``SpectralCubeHistogram`` to represent model or actual counts in energy bands (``gtbin`` format)
"""
from __future__ import print_function, division
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.units import Quantity
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import Angle
from ..spectrum import (LogEnergyAxis,
                        energy_bounds_equal_log_spacing,
                        energy_bin_centers_log_spacing
                        )
from ..spectrum import powerlaw
from ..image.utils import coordinates
from ..irf import EnergyDependentTablePSF
from ..image import cube_to_image

__all__ = ['GammaSpectralCube', 'compute_npred_cube', 'convolve_npred_cube']


class GammaSpectralCube(object):
    """Spectral cube for gamma-ray astronomy.

    Can be used e.g. for Fermi or GALPROP diffuse model cubes.

    Note that there is a very nice ``SpectralCube`` implementation here:
    http://spectral-cube.readthedocs.org/en/latest/index.html

    Here is some discussion if / how it could be used:
    https://github.com/radio-astro-tools/spectral-cube/issues/110

    For now we re-implement what we need here.

    The order of the spectral cube axes can be very confusing ... this should help:

    * The ``data`` array axis order is ``(energy, lat, lon)``.
    * The ``wcs`` object axis order is ``(lon, lat, energy)``.
    * Methods use the ``wcs`` order of ``(lon, lat, energy)``,
      but internally when accessing the data often the reverse order is used.
      We use ``(xx, yy, zz)`` as pixel coordinates for ``(lon, lat, energy)``,
      as that matches the common definition of ``x`` and ``y`` in image viewers.

    Parameters
    ----------
    data : array_like
        Data array (3-dim)
    wcs : `~astropy.wcs.WCS`
        Word coordinate system transformation
    energy : `~astropy.units.Quantity`
        Energy array

    Notes
    -----
    Diffuse model files in this format are distributed with the Fermi Science tools
    software and can also be downloaded at
    http://fermi.gsfc.nasa.gov/ssc/data/access/lat/BackgroundModels.html

    E.g. the 2-year diffuse model that was used in the 2FGL catalog production is at
    http://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/gal_2yearp7v6_v0.fits
    """
    def __init__(self, data, wcs, energy):
        # TODO: check validity of inputs
        self.data = data
        self.wcs = wcs

        # TODO: decide whether we want to use an EnergyAxis object or just use the array directly.
        self.energy = energy
        self.energy_axis = LogEnergyAxis(energy)

        self._interpolate_cache = None

    @property
    def _interpolate(self):
        if self._interpolate_cache == None:
            # Initialise the interpolator
            # This doesn't do any computations ... I'm not sure if it allocates extra arrays.
            from scipy.interpolate import RegularGridInterpolator
            points = list(map(np.arange, self.data.shape))
            self._interpolate_cache = RegularGridInterpolator(points, self.data.value,
                                                              fill_value=None, bounds_error=False)

        return self._interpolate_cache

    @staticmethod
    def read_hdu(hdu_list):
        """Read spectral cube from HDU.

        Parameters
        ----------
        object_hdu : `astropy.io.fits.ImageHDU`
            Image HDU object to be read
        energy_table_hdu : `astropy.io.fits.TableHDU`
            Table HDU object giving energies of each slice
            of the Image HDU object_hdu

        Returns
        -------
        spectral_cube : `GammaSpectralCube`
            Spectral cube
        """
        object_hdu = hdu_list[0]
        energy_table_hdu = hdu_list[1]
        data = object_hdu.data
        data = Quantity(data, '1 / (cm2 MeV s sr)')
        # Note: the energy axis of the FITS cube is unusable.
        # We only use proj for LON, LAT and call ENERGY from the table extension
        header = object_hdu.header
        wcs = WCS(header)
        energy = energy_table_hdu.data['Energy']
        energy = Quantity(energy, 'MeV')
        return GammaSpectralCube(data, wcs, energy)

    @staticmethod
    def read(filename):
        """Read spectral cube from FITS file.

        Parameters
        ----------
        filename : str
            File name

        Returns
        -------
        spectral_cube : `GammaSpectralCube`
            Spectral cube
        """
        data = fits.getdata(filename)
        data = Quantity(data, '1 / (cm2 MeV s sr)')
        # Note: the energy axis of the FITS cube is unusable.
        # We only use proj for LON, LAT and do ENERGY ourselves
        header = fits.getheader(filename)
        wcs = WCS(header)
        energy = Table.read(filename, 'ENERGIES')['Energy']
        energy = Quantity(energy, 'MeV')

        return GammaSpectralCube(data, wcs, energy)

    def world2pix(self, lon, lat, energy, combine=False):
        """Convert world to pixel coordinates.

        Parameters
        ----------
        lon, lat, energy

        Returns
        -------
        x, y, z or array with (x, y, z) as columns
        """
        lon = lon.to('deg').value
        lat = lat.to('deg').value
        x, y, _ = self.wcs.wcs_world2pix(lon, lat, 0, 0)

        z = self.energy_axis.world2pix(energy)

        shape = (x * y * z).shape
        x = x * np.ones(shape)
        y = y * np.ones(shape)
        z = z * np.ones(shape)

        if combine == True:
            x = np.array(x).flat
            y = np.array(y).flat
            z = np.array(z).flat
            return np.column_stack([z, y, x])
        else:
            return x, y, z

    def pix2world(self, x, y, z):
        """Convert world to pixel coordinates.

        Parameters
        ----------
        x, y, z

        Returns
        -------
        lon, lat, energy
        """
        lon, lat, _ = self.wcs.wcs_pix2world(x, y, 0, 0)
        energy = self.energy_axis.pix2world(z)

        lon = Quantity(lon, 'deg')
        lat = Quantity(lat, 'deg')
        energy = Quantity(energy, self.energy.unit)

        return lon, lat, energy

    @property
    def spatial_coordinate_images(self):
        """Spatial coordinate iamges.

        Returns
        -------
        lon, lat : `astropy.units.Quantity`
            Arrays of longitude and latitude pixel coordinates.
        """
        n_lon = self.data.shape[2]
        n_lat = self.data.shape[1]
        i_lat, i_lon = np.indices((n_lat, n_lon))
        lon, lat, _ = self.pix2world(i_lon, i_lat, 0)

        return lon, lat

    @property
    def solid_angle_image(self):
        """Solid angle image.

        TODO: currently uses CDELT1 x CDELT2, which only
              works for cartesian images near the equator.

        Returns
        -------
        solid_angle_image : `~astropy.units.Quantity`
            Solid angle image (steradian)
        """
        cdelt = self.wcs.wcs.cdelt
        solid_angle = np.abs(cdelt[0]) * np.abs(cdelt[1])
        shape = self.data.shape[:2]

        solid_angle = solid_angle * np.ones(shape, dtype=float)
        solid_angle = Quantity(solid_angle, 'deg^2')

        return solid_angle.to('steradian')

    def flux(self, lon, lat, energy):
        """Differential flux (linear interpolation).

        Parameters
        ----------
        lon : `~astropy.units.Quantity` or `~astropy.coordinates.Angle`
            Longitude
        lat : `~astropy.units.Quantity` or `~astropy.coordinates.Angle`
            Latitude
        energy : `~astropy.units.Quantity`
            Energy

        Returns
        -------
        flux : `~astropy.units.Quantity`
            Differential flux (1 / (cm2 MeV s sr))
        """
        # Determine output shape by creating some array via broadcasting
        shape = (lon * lat * energy).shape

        pix_coord = self.world2pix(lon, lat, energy, combine=True)
        values = self._interpolate(pix_coord)
        values = values.reshape(shape)

        return Quantity(values, '1 / (cm2 MeV s sr)')

    def spectral_index(self, lon, lat, energy, dz=1e-3):
        """Power law spectral index.

        A forward finite difference method with step ``dz`` is used along
        the ``z = log10(energy)`` axis.

        Parameters
        ----------
        lon : `~astropy.units.Quantity` or `~astropy.coordinates.Angle`
            Longitude
        lat : `~astropy.units.Quantity` or `~astropy.coordinates.Angle`
            Latitude
        energy : `~astropy.units.Quantity`
            Energy

        Returns
        -------
        spectral_index : `~astropy.units.Quantity`
            Spectral index
        """
        raise NotImplementedError
        # Compute flux at `z = log(energy)`
        pix_coord = self.world2pix(lon, lat, energy, combine=True)
        flux1 = self._interpolate(pix_coord)

        # Compute flux at `z + dz`
        pix_coord[:, 0] += dz
        # pixel_coordinates += np.zeros(pixel_coordinates.shape)
        flux2 = self._interpolate(pix_coord)

        # Power-law spectral index through these two flux points
        # d_log_flux = np.log(flux2 / flux1)
        # spectral_index = d_log_flux / dz
        energy1 = energy
        energy2 = (1. + dz) * energy
        spectral_index = powerlaw.g_from_points(energy1, energy2, flux1, flux2)

        return spectral_index

    def integral_flux_image(self, energy_band, energy_bins=10):
        """Integral flux image for a given energy band.

        A local power-law approximation in small energy bins is
        used to compute the integral.

        Parameters
        ----------
        energy_band : `~astropy.units.Quantity`
            Tuple ``(energy_min, energy_max)``
        energy_bins : int or `~astropy.units.Quantity`
            Energy bin definition.

        Returns
        -------
        image : `~astropy.io.fits.ImageHDU`
            Integral flux image (1 / (m^2 s))
        """
        if isinstance(energy_bins, int):
            energy_bins = energy_bounds_equal_log_spacing(energy_band, energy_bins)

        energy1 = energy_bins[:-1]
        energy2 = energy_bins[1:]

        # Compute differential flux at energy bin edges of all pixels
        xx = np.arange(self.data.shape[2])
        yy = np.arange(self.data.shape[1])
        zz = self.energy_axis.world2pix(energy_bins)

        xx, yy, zz = np.meshgrid(zz, yy, xx, indexing='ij')
        shape = xx.shape

        pix_coords = np.column_stack([xx.flat, yy.flat, zz.flat])
        flux = self._interpolate(pix_coords)
        flux = flux.reshape(shape)

        # Compute integral flux using power-law approximation in each bin
        flux1 = flux[:-1, :, :]
        flux2 = flux[1:, :, :]
        energy1 = energy1[:, np.newaxis, np.newaxis].value
        energy2 = energy2[:, np.newaxis, np.newaxis].value
        integral_flux = powerlaw.I_from_points(energy1, energy2, flux1, flux2)

        integral_flux = integral_flux.sum(axis=0)

        # TODO: check units ... set correctly if not OK.
        header = self.wcs.sub(['longitude', 'latitude']).to_header()
        hdu = fits.ImageHDU(data=integral_flux, header=header, name='integral_flux')

        return hdu

    def __repr__(self):
        # Copied from `spectral-cube` package
        s = "GammaSpectralCube with shape={0}".format(self.data.shape)
        if self.data.unit is u.dimensionless_unscaled:
            s += ":\n"
        else:
            s += " and unit={0}:\n".format(self.data.unit)
        s += " n_x: {0:5d}  type_x: {1:15s}  unit_x: {2}\n".format(self.data.shape[2], self.wcs.wcs.ctype[0], self.wcs.wcs.cunit[0])
        s += " n_y: {0:5d}  type_y: {1:15s}  unit_y: {2}\n".format(self.data.shape[1], self.wcs.wcs.ctype[1], self.wcs.wcs.cunit[1])
        s += " n_s: {0:5d}  type_s: {1:15s}  unit_s: {2}".format(self.data.shape[0], self.wcs.wcs.ctype[2], self.wcs.wcs.cunit[2])
        return s


def compute_npred_cube(flux_cube, exposure_cube, energy_bounds):
    """TODO

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO
    """
    if flux_cube.data.shape[1:] != exposure_cube.data.shape[1:]:
        raise ValueError('flux_cube and exposure cube must have the same shape!\n'
                         'flux_cube: {0}\nexposure_cube: {1}'
                         ''.format(flux_cube.data.shape[1:], exposure_cube.data.shape[1:]))

    energy_centers = energy_bin_centers_log_spacing(energy_bounds)
    wcs = exposure_cube.wcs
    lon, lat = exposure_cube.spatial_coordinate_images
    solid_angle = exposure_cube.solid_angle_image

    npred_cube = np.zeros_like(solid_angle)
    for i in range(len(energy_bounds) - 1):
        energy_bound = energy_bounds[i:i + 2]
        int_flux = flux_cube.integral_flux_image(energy_bound)
        exposure = exposure_cube.flux(lon, lat, energy_centers[i])
        npred_image = int_flux.data * exposure * solid_angle
        npred_cube[i] = npred_image

    npred_cube = GammaSpectralCube(data=npred_cube,
                                   wcs=wcs,
                                   energy=energy_bounds)

    return npred_cube


def convolve_npred_cube(npred_cube, max_offset, resolution=1):
    """TODO

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO
    """
    from scipy.ndimage import convolve

    pixel_size = Angle(resolution, 'deg')
    offset_max = Angle(max_offset, 'deg')
    if energy == 'None':
        fermi_psf = EnergyDependentTablePSF.read(filename)
        # PSF energy band calculation doesn't work, so implemented at log center energy instead
        energy = Quantity(np.sqrt(energy_band[0] * energy_band[1]), 'MeV')
        psf = fermi_psf.table_psf_at_energy(energy)
    else:
        energy = Quantity(energy, 'MeV')
        fermi_psf = EnergyDependentTablePSF.read(filename)
        psf = fermi_psf.table_psf_at_energy(energy=energy)
    psf.normalize()
    kernel = psf.kernel(pixel_size=pixel_size, offset_max=offset_max)
    kernel_image = kernel.value / kernel.value.sum()
    result = convolve(image, kernel_image, mode='constant')

    return result
