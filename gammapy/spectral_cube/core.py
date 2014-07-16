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
from ..spectrum import LogEnergyAxis, energy_bounds_equal_log_spacing
from ..spectrum import powerlaw
from ..image.utils import coordinates

__all__ = ['GammaSpectralCube', 'reproject_cube', 'compute_npred_cube', 'correlate_fermi_psf',
           'interp_flux', 'interp_exposure']


def correlate_fermi_psf(image, max_offset, resolution=0.1, energy = 'None', energy_band=[10, 500]):
    from ..datasets import FermiGalacticCenter
    from ..irf import EnergyDependentTablePSF
    from astropy.coordinates import Angle
    from scipy.ndimage import convolve

    filename = FermiGalacticCenter.filenames()['psf']
    pixel_size = Angle(resolution, 'deg')
    offset_max = Angle(max_offset, 'deg')
    if energy == 'None':
        energy_band = Quantity(energy_band, 'MeV')
        fermi_psf = EnergyDependentTablePSF.read(filename)
        psf = fermi_psf.table_psf_in_energy_band(energy_band=energy_band, spectral_index=2.5)
    else:
        energy = Quantity(energy, 'MeV')
        fermi_psf = EnergyDependentTablePSF.read(filename)
        psf = fermi_psf.table_psf_at_energy(energy=energy)
    psf.normalize()
    kernel = psf.kernel(pixel_size=pixel_size, offset_max=offset_max)
    kernel_image = kernel.value/kernel.value.sum()
    return convolve(image, kernel_image, mode='constant')


def interp_flux(hdu, new_energy):
    cube = GammaSpectralCube.read(hdu)
    lat, lon = coordinates(hdu, world=True, radians=True)

    lat = Quantity(lat, 'rad')
    lon = Quantity(lon, 'rad')

    array = cube.flux(lat, lon, new_energy)

    return array.reshape(lat.shape)


def interp_exposure(hdu, new_energy):
    max_energy = hdu.header['CRPIX3'] + (hdu.header['NAXIS3'] * hdu.header['CDELT3'])
    if new_energy >= max_energy:
        max_index = hdu.header['NAXIS3']
        return hdu.data[max_index]
    else:
        cube = GammaSpectralCube.read(hdu)
        lat, lon = coordinates(hdu, world=True, radians=True)

        lat = Quantity(lat, 'rad')
        lon = Quantity(lon, 'rad')

        array = cube.flux(lat, lon, new_energy)

        return array.reshape(lat.shape)


def reproject_cube(hdu1, hdu2, smooth=False):
    """Reprojects hdu1 to the header of hdu2 and returns as hdu.

    Optionally smooths HDUs to match resolution.
    """
    from FITS_tools.cube_regrid import regrid_cube_hdu
    if smooth == False:
        return regrid_cube_hdu(hdu1, hdu2.header)
    else:
        return regrid_cube_hdu(hdu1, hdu2.header, smooth=True)


def compute_npred_cube(flux_hdu, exposure_hdu, desired_energy=None, convolve='Yes', max_convolution_offset=10):
    """ Computes predicted counts cube from model flux cube at given energy.

    TODO: Update this docstring.

    Interpolates/Extrapolates flux cube and exposure cube to the same energies
    assuming a power law index, multiplying and optionally convolving by the
    Fermi/LAT (or provided) energy dependent PSF and returns either counts cube
    (for specified energy range) or counts image at a specified energy value.

    Parameters
    ----------
    flux_cube : array_like
        Data array (3-dim)
    energy : float (optional)
        Energy at which npred is returned, MeV
    psf : array_like (optional)
        Data array PSF image (2-dim)
    e_bounds : array_like
        Data array for energy bounds (1-dim), MeV
    n_energies : int        
        Number of equal sized logarithmic energy bins between the
        maximum and minimum energies specified in e_bounds

    Returns
    -------
    npred_cube : array_like
        Data array (2-dim or 3-dim). Either a data cube if energy range is
        specified, or an image slice of the data cube if single energy is
        specified

    Notes
    -----
    Returns the cube as a 3-dimensional array rather than a hdu object.
    """
    flux_rep = reproject_cube(flux_hdu, exposure_hdu, smooth=True)

    if desired_energy != None:
        flux = interp_flux(flux_rep, desired_energy)
        exposure = interp_exposure(exposure_hdu, desired_energy)
    else:
        flux = flux_rep.data
        exposure = exposure.data

    flux *= exposure

    if convolve == 'Yes':
        energy_bins = flux.header['NAXIS3']
        energy_min = flux.header['CRPIX3']
        bin_size = flux.header['CDELT3']
        resolution = flux.header['CDELT2']
        indices = np.arange(energy_bins)
        counts = flux.copy()
        for index in indices:
            #energy = Quantity(energy_min + index * bin_size, 'MeV')
            energy = energy_min + index * bin_size
            cube_layer = flux[index]
            eband_min = energy_min + index * bin_size
            eband_max = energy_min + (index + 1) * bin_size
            counts[index] = correlated_fermi_psf(cube_layer, max_convolution_offset,
                                                 resolution, energy, energy_band=[eband_min, eband_max])
    else:
        counts = flux
    return fits.ImageHDU(counts)


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
        values.reshape(shape)

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
