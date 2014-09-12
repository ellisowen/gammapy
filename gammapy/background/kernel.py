# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function, division
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import gc
from ..stats import significance
from ..image import binary_dilation_circle


__all__ = ['GammaImages', 'IterativeKernelBackgroundEstimator']


class GammaImages(object):
    """TODO: implement a more general images container class
    that can be re-used in other places as well.
    """

    def __init__(self, counts, background=None, mask=None):
        self.counts = np.asarray(counts.data, dtype=float)
        self.header = counts.header
        if background == None:
            # Start with a flat background estimate
            self.background = 0.001 * np.ones_like(counts.data, dtype=float)
        else:
            self.background = np.asarray(background.data, dtype=float)

        if mask == None:
            self.mask = np.ones_like(counts.data, dtype=bool)
        else:
            self.mask = np.asarray(mask.data, dtype=bool)
    
    def compute_correlated_maps(self, kernel):
        """Compute significance image for a given kernel.
        """
        from scipy.ndimage import convolve
        print("Calls Images!")
        self.counts_corr = convolve(self.counts, kernel)
        self.background_corr = convolve(self.background, kernel)
        self.significance = significance(self.counts_corr, self.background_corr)
        return self


class IterativeKernelBackgroundEstimator(object):
    """Iteratively estimate a background model.

    Parameters
    ----------
    images : `~gammapy.background.GammaImages`
        GammaImages object containing counts image and (optional) initial
        background estimation.
    source_kernel : `numpy.ndarray`
        Source kernel as a numpy array.
    background_kernel : `numpy.ndarray`
        Background convolution kernel as a numpy array.
    significance_threshold : float
        Significance threshold above which regions are excluded.
    mask_dilation_radius : float
        Amount by which mask is dilated with each iteration.
    delete_intermediate_results : bool
        Specify whether results of intermediate iterations should be deleted.
        (Otherwise, these are held in memory). Default True.
    save_intermediate_results : bool
        Specify whether to save intermediate results as FITS files to disk.
        Default False.
    filebase : str (optional)
        Base of filenames if save_intermediate_results = True. Default 'temp'.

    Returns
    -------
    mask : `~astropy.io.fits.ImageHDU`
        Final exclusion mask, when iterations are complete.
    background : `~astropy.io.fits.ImageHDU`
        Final background estimation, when iterations are complete.
    """

    def __init__(self, images, source_kernel, background_kernel,
                 significance_threshold, mask_dilation_radius,
                 delete_intermediate_results=True,
                 save_intermediate_results=False, filebase='temp'):
        
        # self._data[i] is a GammaImages object representing iteration number `i`.
        self._data = list()
        self._data.append(images)
        
        self.header = images.header
        
        self.source_kernel = source_kernel
        self.background_kernel = background_kernel

        self.significance_threshold = significance_threshold
        self.mask_dilation_radius = mask_dilation_radius
        
        self.delete_intermediate_results = delete_intermediate_results
        self.save_intermediate_results = save_intermediate_results
        # Calculate initial significance image
        self._data[-1].compute_correlated_maps(self.source_kernel)
        gc.collect()
    
    def run_ntimes(self, n_iterations, filebase=None):
        """Run N iterations."""

        if self.save_intermediate_results:
            self.save_files(filebase, index=0)

        for ii in range(1, n_iterations + 1):
            if ii == 1:
                # This is needed to avoid excluding the whole Galactic plane
                # in case the initial background estimate is much too low.
                update_mask = False
            else:
                update_mask = True
            
            self.run_iteration(update_mask)
            
            if self.save_intermediate_results: 
                self.save_files(filebase, index=ii)

            if self.delete_intermediate_results:
                # Remove results from previous iteration
                del self._data[0]
                gc.collect()

            mask = self._data[0].mask
            background = self._data[0].background

        return mask, background


    def run(self, filebase=None):
        """Run until mask is stable."""

        ii = 0
        if self.save_intermediate_results:
            self.save_files(filebase, index=ii)
        
        self.run_iteration(update_mask=False)
            
        ii = ii + 1
        if self.save_intermediate_results: 
            self.save_files(filebase, index=ii)
        
        check_mask = False

        while check_mask == False:

            ii = ii + 1

            self.run_iteration(update_mask=True)
            
            if self.save_intermediate_results: 
                self.save_files(filebase, index=ii)

            new_mask = self._data[-1].mask
            old_mask = self._data[0].mask
            background = self._data[-1].mask
            if self.delete_intermediate_results:
                # Remove results from previous iteration
                del self._data[0]
                gc.collect()

            # Test mask (note will terminate as soon as mask does not change).
            new = new_mask.sum()
            old = old_mask.sum()

            if new == old:
                check_mask = True
            else:
                check_mask = False

        return new_mask, background


    def run_iteration(self, update_mask=True):
        """Run one iteration.

        Parameters
        ----------
        update_mask : bool
            Specify whether to update the exclusion mask stored in the input
            `~gammapy.background.GammaImages` object with the exclusion mask
            newly calculated in this method.
        """

        from scipy.ndimage import convolve
        # Start with images from the last iteration. If not, makes one.
        # Check if initial mask exists:
        if self._data[0].mask.sum() == self._data[-1].mask.size:
            mask = np.where(self._data[0].significance > self.significance_threshold, 0, 1)
        else:
            # Compute new exclusion mask:
            if update_mask:
                # Check first if mask needs updating:
                check = np.where(self._data[0].significance[self._data[-1].mask] > self.significance_threshold, 1, 0)
                if check.sum() == 0:
                    mask = self._data[-1].mask.copy()
                else:
                    # Only update mask if there are further significant regions to exclude
                    mask = np.invert(binary_dilation_circle(np.invert(self._data[-1].mask), radius=self.mask_dilation_radius))
            else:
                mask = self._data[-1].mask.copy()
        
        # Compute new background estimate:
        # Convolve old background estimate with background kernel,
        # excluding sources via the old mask.
        weighted_counts = convolve(mask * self._data[0].counts, self.background_kernel)
        weighted_counts_normalisation = convolve(mask.astype(float), self.background_kernel)
        background = weighted_counts / weighted_counts_normalisation
        
        # Convert new Images to HDUs to store in a GammaImages object
        counts = fits.ImageHDU(data=self._data[-1].counts, header=self._data[-1].header)
        background = fits.ImageHDU(data=background, header=self._data[-1].header)
        mask = fits.ImageHDU(data=mask.astype(int), header=self._data[-1].header)
        images = GammaImages(counts, background, mask)
        images.compute_correlated_maps(self.source_kernel)
        self._data.append(images)
        significance = fits.ImageHDU(data=self._data[-1].significance, header=self._data[-1].header)
    
    def save_files(self, filebase, index):
        """Saves files to fits."""

        header = self.header

        filename = filebase + '{0:02d}_mask'.format(index) + '.fits'
        hdu = fits.ImageHDU(data=self._data[-1].mask.astype(np.uint8), header=header)
        hdu.writeto(filename, clobber=True)

        filename = filebase + '{0:02d}_background'.format(index) + '.fits'
        hdu = fits.ImageHDU(data=self._data[-1].background, header=header)
        hdu.writeto(filename, clobber=True)

        filename = filebase + '{0:02d}_significance'.format(index) + '.fits'
        hdu = fits.ImageHDU(data=self._data[-1].significance, header=header)
        hdu.writeto(filename, clobber=True)    

    @property
    def mask_image_hdu(self):
        """Returns mask as `~astropy.io.fits.ImageHDU`."""

        header = self.header
        #import IPython; IPython.embed() 
        return fits.ImageHDU(data=self._data[-1].mask.astype(np.uint8), header=header)    

    @property
    def background_image_hdu(self):
        """Returns resulting background estimate as `~astropy.io.fits.ImageHDU`."""

        header = self.header

        return fits.ImageHDU(data=self._data[-1].background, header=header)

    @property
    def significance_image_hdu(self):
        """Returns resulting background estimate as `~astropy.io.fits.ImageHDU`."""

        header = self.header

        return fits.ImageHDU(data=self._data[-1].significance, header=header)
