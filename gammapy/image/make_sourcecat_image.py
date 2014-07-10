""" make an image from a source catalog, e.g 1FHL 2FGL etc

In progress... will later add to gammapy
"""
__all__ = ['make_reference_survey_image', 'fetch_fermi_extended_sources', 'correlate_fermi_psf', 'correlate_gauss', 'make_model_image', 'get_table_from_fits']

from astropy.io import fits
import time, sys
import numpy as np
from numpy import sign, sqrt, log
from scipy.ndimage import convolve
import pyfits
from gammapy.image import block_reduce_hdu

def make_reference_survey_image(resolution, latitudes, longitudes, center=[0,0], units='ph/cm2/s/sr'):
    """resolution: float
    latitudes, longitudes: [min, max]
    center: [center_lat, center_lon]
    units: str
    """
    import datetime
    today = datetime.date.today().strftime("%B %d, %Y")

    lat_pix = (latitudes[1]-latitudes[0])/resolution
    lon_pix = (longitudes[1]-longitudes[0])/resolution
    pars = {'SIMPLE': 'T', 'BITPIX':-32, 'NAXIS': 2,
           'NAXIS1': lon_pix,
           'NAXIS2': lat_pix, 'EXTEND': 'T', 'DATE': today , 'COMMENT': 0, 'COMMENT': 0, 'CTYPE1': 'GLON-CAR', 'CTYPE2': 'GLAT-CAR',
           'EQUINOX': 2000.00, 'CDELT1':(-1 * resolution), 'CDELT2': resolution , 'CROTA2': 0, 'CRPIX1': np.floor(0.5*lon_pix)+0.5,
           'CRPIX2': np.floor(0.5*lat_pix)+0.5, 'CRVAL1': center[1], 'CRVAL2': center[0] ,
           'BUNIT': units , 'PEDVAL': 0, 'POWER': 0
           }
    
    header = fits.Header()
    header.update(pars)
                       
    shape = (header['NAXIS2'], header['NAXIS1'])
    data = np.zeros(shape)
    return fits.PrimaryHDU(data, header)

def fetch_fermi_extended_sources(catalog):
    from astropy.utils import data
    import tarfile
    import os
    # Note: not actually used
    """Gets fermi extended sources fits data for a catalog:
    """
    BASE_URL = 'http://fermi.gsfc.nasa.gov/ssc/data/access/lat/'
    if catalog == '2FGL':
        url = BASE_URL + '2yr_catalog/gll_psc_v07_templates.tgz'
    elif catalog == '1FGL':
        raise NotImplementedError
    elif catalog == '1FHL':
        url = BASE_URL + '1FHL/LAT_extended_sources_v12.tar'
    elif catalog == '2PC':
        raise NotImplementedError
    else:
        ss = 'Invalid catalog: {0}\n'.format(catalog)
        ss += 'Available: {0}'.format(', '.join(FERMI_CATALOGS))
        raise ValueError(ss)
    # Downloads the relevant tarball
    directory = data.download_file(url, cache=True)
    # Unzips the tarball
    tar = tarfile.open(directory)
    tar.extractall()
    tar.close()
    # creates a list of fits pathsways in the unzipped directory
    os.system('ls -d -1 $PWD/Templates/** > pathways.txt')
    # Creates a fits hdu list of all fits files in the directory
    inputFn = "pathways.txt"
    lines = file(inputFn, 'r').readlines()
    length = len(lines)
    entries = np.arange(length)
    
    hdu_list = []
    for entry in entries:
        line = lines[entry][:-1]
        hdu = fits.open(line)[0]
        hdu_list.append(hdu)
    headers = fits.HDUList(hdu_list)
    
    return headers

from gammapy.image.utils import disk_correlate

def correlate_fermi_psf(image, max_offset, resolution=0.1, energy = 'None', energy_band=[10, 500]):
    from astropy.coordinates import Angle
    from astropy.units import Quantity
    from gammapy.datasets import FermiGalacticCenter
    from gammapy.irf import EnergyDependentTablePSF

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

def correlate_gauss(image, radius):
    from gammapy.morphology import Gauss2DPDF
    """Correlate image with gaussian of a given radius.
    
    This is also called "gaussian correlation" and it means that
    the value of a given pixel in the output image is the
    integral of all pixels within a given input range convolved with a gaussian kernel.
    """
    radius = int(radius)
    gaussian = Gauss2DPDF(radius)
    y, x = np.mgrid[-radius: radius + 1, -radius: radius + 1]
    structure = gaussian(x, y)
    return convolve(image, structure, mode='constant')

def make_model_image(psf='Fermi', resolution=0.1, center=[0, 0], lat_range=[0, 180], lon_range=[0, 360], catalog='1FHL', total_flux='False', filename='1fhl_fermi_psf.fits'):  
    from gammapy.image import coordinates
    from astropy.convolution import convolve
    from astropy.modeling.models import Gaussian2D, Disk2D
    from gammapy.image import images_to_cube
    from gammapy.image import paste_cutout_into_image
    from gammapy.image.utils import solid_angle
    
    reference = make_reference_survey_image(resolution, lat_range, lon_range, center, units='ph/cm2/s/sr') #Check these units!
    
    lons, lats = coordinates(reference)
     
    source_table = get_table_from_fits(catalog, extended='No')
    print source_table
    sources = np.arange(len(source_table['Flux'])).astype(int)
    
    if psf == 'None':
        # If there is no PSF defined, sources will be modelled as Gaussians
        source_kernel = Gaussian2D(0, 0, 0, 0.1, 0.1)
        new_image = im_1 = im_2 = im_3 = source_kernel(lats, lons)
    else:
        # Otherwise, all of the flux will be placed into the reference pixel to be later PSF convolved with the defined PSF
        # Hence original reference empty image is called
        new_image = im_1 = im_2 = im_3 = reference.data
    print len(sources)
    total_point_image = fits.ImageHDU(header=reference.header, data=new_image)
    from astropy.wcs import WCS
    wcs = WCS(total_point_image.header)
    new_image = np.zeros_like(total_point_image.data, dtype=np.float64)
    for source in sources:
        source_type = 'PointSource'#source_table['Source_Type'][source]
        print source
        if source_type == 'ExtendedSource':
            raise NotImplementedError
            # This needs more work...
            #image = source_table['Image'][source]
            #image.data = (image.data * solid_angle(image).value.mean()) / 1000000000000  # TODO: fix this hack... units???
            #resample_factor1 = np.round(reference.header['CDELT1'] / image.header['CDELT1'])
            #resample_factor2 = np.round(reference.header['CDELT2'] / image.header['CDELT2'])
            #block_factors = np.array([resample_factor1, resample_factor2])  # TODO: fix this approximation... kapteyn image utils reprojectto?
            #resampled_image = block_reduce_hdu(image, block_factors, np.sum)
            #paste_cutout_into_image(total_point_image, resampled_image)
            
        elif source_type == 'PointSource':
            lon = source_table['GLON'][source]
            lat = source_table['GLAT'][source]
            flux = source_table['Flux'][source]
            precise = False
            if precise:
                raise NotImplementedError
            else:
                #print(lon, lat)
                #import IPython; IPython.embed()
                print lon, lat
                x, y = wcs.wcs_world2pix(lon, lat, 0)
                print x, y
                
                # TODO: check if this is 0.5 or 1 pix off
                xi, yi = x.astype(int), y.astype(int)
                
                new_image[yi, xi] += flux
            
            
    total_point_image = fits.ImageHDU(header=reference.header, data=new_image)
    
    # Ensure flux or counts remain the same
    if total_flux == 'True':
        factor = source_table['Flux'].sum()
    else:
        factor = total_flux
        
    if psf == 'None':
        new_image = (new_image/new_image.sum()) * factor
    elif psf == 'Gaussian':
        new_image = correlate_gauss((new_image/new_image.sum()), 2)
    elif psf == 'Disk':
        new_image = disk_correlate((new_image/new_image.sum()), 2)
    elif psf == 'Fermi':
        print "Original Flux"
        print new_image.sum()
        new_image = correlate_fermi_psf(new_image, 5)
        print "PSF Convolved Flux"
        print new_image.sum()
    header = reference.header
    image = fits.ImageHDU(data=new_image, header=header)
    image.writeto(filename, clobber=True)
    
    # TODO: Implement these lines as part of energy binning
    #factor_1 = source_table['Flux10_30GeV'].sum()/im_1.sum()
    #im_1 =  correlate_gauss((im_1) * factor_1, 3)
    #image1 = fits.ImageHDU(data=im_1, header=header)
    #factor_2 = source_table['Flux30_100GeV'].sum()/im_2.sum()
    #im_2 =  correlate_gauss((im_2) * factor_2, 3)
    #image2 = fits.ImageHDU(data=im_2, header=header)
    #factor_3 = source_table['Flux100_500GeV'].sum()/im_3.sum()
    #im_3 =  correlate_gauss((im_3) * factor_3, 3)
    #image3 = fits.ImageHDU(data=im_3, header=header)
    
    # Also for energy binning
    # hdu_list = [image1, image2, image3]
    #cube_hdu = images_to_cube(hdu_list)
    #cube_hdu.writeto('1fhl_cube_gauss_2.fits', clobber=True) 
           
def get_table_from_fits(catalog, extended='No'):
    
    from gammapy.datasets import fetch_fermi_catalog
    from astropy.table import Table
    from scipy.stats import gmean
    
    data = []
    #TODO: implement multiple catalog feature (low priority)
    if catalog == '1FHL':
        point_sources = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')
    elif catalog == 'TeVCat':
        point_sources = fits.open('TeVCat_all.fits')[1].data
    n_sources = len(point_sources)
    sources = np.arange(n_sources)
    print sources
    for source in sources:
        print source
        glon = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['GLON'][source]
        glat = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['GLAT'][source]
        if catalog == '1FHL':
            flux_headers = ['Flux10_30GeV', 'Flux30_100GeV', 'Flux100_500GeV']
            flux_bol = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Flux'][source]
            #NOTE: Implementation of energy bands is not yet complete and doesn't appear in the output
            Flux10_30GeV = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Flux10_30GeV'][source]
            Flux30_100GeV = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Flux30_100GeV'][source]
            Flux100_500GeV = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Flux100_500GeV'][source]
        elif catalog == 'TeVCat':
            flux_headers = ['Flux10_30GeV', 'Flux30_100GeV', 'Flux100_500GeV']
            flux_bol = fits.open('TeVCat_all.fits')[1].data['Flux_CU'][source]

            #NOTE: Implementation of energy bands is not yet complete and doesn't appear in the output
            Flux10_30GeV = 0 #fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Flux10_30GeV'][source]
            Flux30_100GeV = 0# fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Flux30_100GeV'][source]
            Flux100_500GeV = 0# fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Flux100_500GeV'][source]
        elif catalog == '2FGL':
            flux_headers = ['Flux30_100', 'Flux100_300', 'Flux300_1000', 'Flux1000_3000', 'Flux3000_10000', 'Flux10000_100000'] 
            # Update this to reflect the table headers...
            flux_bol = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Flux_Density'][source]
            Flux10_30GeV = 0#fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Flux10_30'][source]
            Flux30_100GeV = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Flux30_100'][source]
            Flux100_500GeV = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Flux100_300'][source]  #!!! TODO: Fix this hack
        source_name = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')['Source_Name'][source]
        row = dict(Source_Type='PointSource', Source_Name=source_name, GLON=glon, GLAT=glat, Flux=flux_bol, Flux10_30GeV=Flux10_30GeV, Flux30_100GeV=Flux30_100GeV, Flux100_500GeV=Flux100_500GeV, Image='None')
        data.append(row)
    if extended=='Yes':
        raise NotImplementedError
        ext_sources = fetch_fermi_catalog(catalog, 'ExtendedSources')  # Change 2fgl to catalog here when fixed TODO:
        n_sources = len(ext_sources)
        sources = np.arange(n_sources)
        for source in sources:
            # glon = fetch_fermi_catalog(catalog, 'ExtendedSources')['GLON'][source]
            # glat = fetch_fermi_catalog(catalog, 'ExtendedSources')['GLAT'][source]
            #if catalog == '1FHL':
            # flux_headers = ['Flux10_30GeV', 'Flux30_100GeV', 'Flux100_500GeV']
            image_name = fetch_fermi_catalog(catalog, 'ExtendedSources')['Spatial_Filename'][source]  # TODO: Fix this hack
            flux_image = fits.open('Templates/{0}'.format(image_name))[0]
            print flux_image.header
            row = dict(Source_Type='ExtendedSource', Source_Name='None', GLON=0, GLAT=0, Flux=0, Flux10_30GeV=0, Flux30_100GeV=0, Flux100_500GeV=0, Image=flux_image)
            data.append(row)
        
    table = Table(data)
        # data required to pass on:
            # point source flux
            # point source lat, lon
        
        # extended_sources = fetch_fermi_catalog(catalog, 'ExtendedSources')
        
        # data required to pass on:
            # extended source peak flux
            # lat, lon of center of the extended source
            # Model Form - better to reconstruct using the astropy models than use the fits files in the tarball
    return table

#if __name__ == '__main__':
#    make_model_image(psf='Fermi', resolution=0.1, center=[0, 0], lat_range=[0, 180], lon_range=[0, 360], catalog='1FHL', total_flux='True', filename='1fhl_fermi_psf.fits')

