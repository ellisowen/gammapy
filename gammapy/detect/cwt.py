# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function, division
import logging
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS

__all__ = ['CWT']


def gauss_kernel(radius, n_sigmas=8):
    """Normalized 2D gauss kernel array.
    """
    sizex = int(n_sigmas * radius)
    sizey = int(n_sigmas * radius)
    radius = float(radius)
    xc = 0.5 * sizex
    yc = 0.5 * sizey
    y, x = np.mgrid[0:sizey - 1, 0:sizex - 1]
    x = x - xc
    y = y - yc
    x = x / radius
    y = y / radius
    g = np.exp(-0.5 * (x ** 2 + y ** 2))
    return g / (2 * np.pi * radius ** 2)  # g.sum()


def difference_of_gauss_kernel(radius, scale_step, n_sigmas=8):
    """Difference of 2 Gaussians (i.e. Mexican hat) kernel array.
    """
    sizex = int(n_sigmas * scale_step * radius)
    sizey = int(n_sigmas * scale_step * radius)
    radius = float(radius)
    xc = 0.5 * sizex
    yc = 0.5 * sizey
    y, x = np.mgrid[0:sizey - 1, 0:sizex - 1]
    x = x - xc
    y = y - yc
    x1 = x / radius
    y1 = y / radius
    g1 = np.exp(-0.5 * (x1 ** 2 + y1 ** 2))
    g1 = g1 / (2 * np.pi * radius ** 2)  # g1.sum()
    x1 = x1 / scale_step
    y1 = y1 / scale_step
    g2 = np.exp(-0.5 * (x1 ** 2 + y1 ** 2))
    g2 = g2 / (2 * np.pi * radius ** 2 * scale_step ** 2)  # g2.sum()
    return g1 - g2


class CWT(object):
    """Continuous wavelet transform.

    TODO: describe algorithm

    TODO: give references

    Initialization of wavelet family.

    Parameters
    ----------
    min_scale : float
        first scale used
    nscales : int
        number of scales considered
    scale_step : float
        base scaling factor
    """
    def __init__(self, min_scale, nscales, scale_step):
        self.kernbase = dict()
        self.scales = dict()
        self.nscale = nscales
        self.scale_step = scale_step
        for ns in np.arange(0, nscales):
            scale = min_scale * scale_step ** ns
            self.scales[ns] = scale
            self.kernbase[ns] = difference_of_gauss_kernel(scale, scale_step)

        # TODO: do we need self.scales and self.scale_array?
        self.scale_array = (scale_step ** (np.arange(0, nscales))) * min_scale

        max_scale = min_scale * scale_step ** nscales
        self.kern_approx = gauss_kernel(max_scale)

#        self.transform = dict()
#        self.error = dict()
#        self.support = dict()

        self.header = None
        self.wcs = None

    def set_data(self, image, background):
        """Set input images."""
        # TODO: check that image and background are consistent 
        self.image = image - 0.0
        self.nx, self.ny = self.image.shape
        self.filter = np.zeros((self.nx, self.ny)) 
        self.background = background - 0.0  # hack because of some bug with old version of fft in numpy
        self.model = np.zeros((self.nx, self.ny)) 
        self.approx = np.zeros((self.nx, self.ny))

        self.transform = np.zeros((self.nscale, self.nx, self.ny))
        self.error = np.zeros((self.nscale, self.nx, self.ny)) 
        self.support = np.zeros((self.nscale, self.nx, self.ny))

    def set_file(self, filename):
        """Set input images from FITS file"""
        # TODO: check the existence of extensions
        # Open fits files
        hdulist = fits.open(filename)
        # TODO: don't hardcode extension numbers and names here ... pass on from gp-cwt
        self.set_data(hdulist[0].data, hdulist['NormOffMap'].data)
        self.header = hdulist[0].header
        self.wcs = WCS(self.header)

    def do_transform(self):
        """Do the transform itself."""
        # TODO: after unit tests are added switch to astropy fftconvolve here.
        from scipy.signal import fftconvolve
        total_background = self.model + self.background + self.approx
        excess = self.image - total_background
        for key, kern in self.kernbase.items():
            self.transform[key] = fftconvolve(excess, kern, mode='same')
            self.error[key] = np.sqrt(fftconvolve(total_background, kern ** 2, mode='same'))

        self.approx = fftconvolve(self.image - self.model - self.bkg,
                                  self.kern_approx, mode='same')
        self.approx_bkg = fftconvolve(self.bkg, self.kern_approx, mode='same')

    def compute_support_peak(self, nsigma=2.0, nsigmap=4.0, remove_isolated=True):
        """Compute the multiresolution support with hard sigma clipping.

        Imposing a minimum significance on a connex region of significant pixels
        (i.e. source detection)
        """ 
        from scipy.ndimage import label
        # TODO: check that transform has been performed
        sig = self.transform / self.error
        for key in self.scales.keys():
            tmp = sig[key] > nsigma
            # produce a list of connex structures in the support
            l, n = label(tmp)
            for id in range(1, n):
                index = np.where(l == id)
                if remove_isolated:
                    if index[0].size == 1:
                        tmp[index] *= 0.0  # Remove isolated pixels from support
                signif = sig[key][index]
                if signif.max() < nsigmap:  # Remove significant pixels island from support 
                    tmp[index] *= 0.0  # if max does not reach maximal significance

            self.support[key] += tmp
            self.support[key] = self.support[key] > 0.

    def inverse_transform(self):
        """Do the inverse transform (reconstruct the image)."""
        res = np.sum(self.support * self.transform, 0)
        self.filter += res * (res > 0)
        self.model = self.filter
        return res

    def iterative_filter_peak(self, nsigma=3.0, nsigmap=4.0, niter=2, convergence=1e-5):
        """Run iterative filter peak algorithm."""
        var_ratio = 0.0
        for iiter in range(niter):
            self.do_transform()
            self.compute_support_peak(nsigma, nsigmap)
            res = self.inverse_transform()
            residual = self.image - (self.model + self.approx)
            tmp_var = residual.var()
            if iiter > 0:
                var_ratio = abs((self.residual_var - tmp_var) / self.residual_var)
                if var_ratio < convergence:
                    logging.info("Convergence reached at iteration {0}".format(iiter + 1))
                    return res
            self.residual_var = tmp_var
        logging.info("Convergence not formally reached at iteration {0}".format(iiter + 1))
        logging.info("Final convergence parameter {0}. Objective was {1}."
                     "".format(convergence, var_ratio)) 
        return res

    def max_scale_image(self):
        """Compute the maximum scale image."""
        maximum = np.argmax(self.transform, 0)
        return self.scale_array[maximum] * (self.support.sum(0) > 0)

    def save_filter(self, filename, clobber=False):
        """Save filter to file."""
        hdu = fits.PrimaryHDU(self.filter, self.header)
        hdu.writeto(filename, clobber=clobber)
        fits.append(filename, self.approx, self.header)
        fits.append(filename, self.filter + self.approx, self.header)
        fits.append(filename, self.max_scale_image(), self.header)


class GammaImages(object):
    """Container for a set of related images.
    
    Meaning of mask:
    * 1 = background region
    * 0 = source region
    (such that multiplying with the mask zeros out the source regions)

    TODO: document
    """
    def __init__(self, counts, background=None, mask=None):
        self.counts = np.asarray(counts, dtype=float)

        if background == None:
            # Start with a flat background estimate
            self.background = np.ones_like(background, dtype=float)
        else:
            self.background = np.asarray(background, dtype=float)

        if mask == None:
            self.mask = np.ones_like(counts, dtype=bool)
        else:
            self.mask = np.asarray(mask, dtype=bool)
    
    def compute_correlated_maps(self, kernel):
        """Compute significance image for a given kernel.
        """
        #import IPython; IPython.embed()
        self.counts_corr = convolve(self.counts, kernel)
        self.background_corr = convolve(self.background, kernel)
        self.significance = significance(self.counts_corr, self.background_corr)

        return self

    def print_info(self):
        logging.info('Counts sum: {0}'.format(self.counts.sum()))
        logging.info('Background sum: {0}'.format(self.background.sum()))
        background_fraction = 100. * self.background.sum() / self.counts.sum()
        logging.info('Background fraction: {0}'.format(background_fraction))
        excluded_fraction = 100. * (1 - np.mean(self.mask))
        logging.info('Excluded fraction: {0}%'.format(excluded_fraction))
    
    def save(self, filename):
        logging.info('Writing {0}'.format(filename))


class IterativeBackgroundEstimator(object):
    """Iteratively estimate a background model.

    see: SciNeGHE source_diffuse_estimation.py
    TODO: document

    Parameters
    ----------
    image : `GammaImages`
        Gamma images

    See also
    --------
    `gammapy.detect.CWT`
    """
    def __init__(self, images, source_kernel, background_kernel,
                 significance_threshold, mask_dilation_radius,
                 delete_intermediate_results=True):
        
        # self._data[i] is a GammaImages object representing iteration number `i`.
        self._data = list()
        self._data.append(images)
        
        self.source_kernel = source_kernel
        self.background_kernel = background_kernel

        self.significance_threshold = significance_threshold
        self.mask_dilation_radius = mask_dilation_radius
        
        self.delete_intermediate_results = delete_intermediate_results
        
        gc.collect()
    
    def run(self, n_iterations, filebase):
        """Run N iterations."""

        self._data[-1].compute_correlated_maps(self.source_kernel)
        self.save_files(filebase, index=0)

        for ii in range(1, n_iterations + 1):
            logging.info('Running iteration #{0}'.format(ii))
            if ii == 1:
                # This is needed to avoid excluding the whole Galactic plane
                # in case the initial background estimate is much too low.
                update_mask = False
            else:
                update_mask = True
            
            self.run_iteration(update_mask)

            self.save_files(filebase, index=ii)

            if self.delete_intermediate_results:
                # Remove results from previous iteration
                del self._data[0]
                gc.collect()

    def run_iteration(self, update_mask=True):
        """Run one iteration."""
        # Start with images from the last iteration
        images = self._data[-1]
        
        logging.info('*** INPUT IMAGES ***')
        images.print_info()

        # Compute new exclusion mask:
        if update_mask:
            logging.info('Computing new exclusion mask')
            mask = np.where(images.significance > self.significance_threshold, 0, 1)
            #print('===', (mask == 0).sum())
            mask = np.invert(binary_dilation_circle(mask == 0, radius=self.mask_dilation_radius))
            #print('===', (mask == 0).sum())
        else:
            mask = images.mask.copy()
        
        # Compute new background estimate:
        # Convolve old background estimate with background kernel,
        # excluding sources via the old mask.
        weighted_counts = convolve(images.mask * images.counts, self.background_kernel)
        weighted_counts_normalisation = convolve(images.mask.astype(float), self.background_kernel)
        background = weighted_counts / weighted_counts_normalisation
        
        # Store new images
        images = GammaImages(counts, background, mask)
        logging.info('Computing source kernel correlated images.')
        images.compute_correlated_maps(self.source_kernel)

        logging.info('*** OUTPUT IMAGES ***')
        images.print_info()
        self._data.append(images)
    
    def save_files(self, filebase, index):

        # TODO: header should be stored as class member instead
        # This is a hack:
        header = fits.getheader('sources.fits.gz', 1)

        images = self._data[-1]

        filename = filebase + '{0:02d}_mask'.format(index) + '.fits'
        logging.info('Writing {0}'.format(filename))
        hdu = fits.ImageHDU(data=images.mask.astype(np.uint8), header=header)
        hdu.writeto(filename, clobber=True)

        filename = filebase + '{0:02d}_background'.format(index) + '.fits'
        logging.info('Writing {0}'.format(filename))
        hdu = fits.ImageHDU(data=images.background, header=header)
        hdu.writeto(filename, clobber=True)

        filename = filebase + '{0:02d}_significance'.format(index) + '.fits'
        logging.info('Writing {0}'.format(filename))
        hdu = fits.ImageHDU(data=images.significance, header=header)
        hdu.writeto(filename, clobber=True)
