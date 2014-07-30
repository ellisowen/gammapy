# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" make an image from a source catalog, or simulated catalog, e.g 1FHL 2FGL etc
"""
__all__ = ['reference_image', 'catalog_image']

from astropy.io import fits
from astropy.modeling.models import Gaussian2D, Disk2D
from ..image.utils import make_header, disk_correlate, fermipsf_correlate, gauss_correlate
from ..datasets.load import fetch_fermi_extended_soruces
from ..image import coordinates
from ..catalog.utils import catalog_table


def reference_image(resolution, latitudes, longitudes, center=[0,0], units='ph/cm2/s/sr'):
    """TODO
    """
    lat_pix = (latitudes[1]-latitudes[0])/resolution
    lon_pix = (longitudes[1]-longitudes[0])/resolution
    header = make_header(lon_pix, lat_pix, resolution, center[0], center[1],
                         'CAR', 'GAL', np.floor(0.5*lon_pix)+0.5, np.floor(0.5*lat_pix)+0.5)
    shape = (header['NAXIS2'], header['NAXIS1'])
    data = np.zeros(shape)
    return fits.PrimaryHDU(data, header)


def catalog_image(psf='Fermi', resolution=0.1, center=[0, 0], lat_range=[0, 180], lon_range=[0, 360], catalog='1FHL', total_flux='False', filename='1fhl_fermi_psf.fits'):  
    """TODO
    """
    reference = reference_image(resolution, lat_range, lon_range, center, units='ph/cm2/s/sr') #Check these units!
    lons, lats = coordinates(reference)
    source_table = catalog_table(catalog, extended='No')
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
        new_image = gauss_correlate((new_image/new_image.sum()), 2)
    elif psf == 'Disk':
        new_image = disk_correlate((new_image/new_image.sum()), 2)
    elif psf == 'Fermi':
        print "Original Flux"
        print new_image.sum()
        new_image = fermipsf_correlate(new_image, 5)
        print "PSF Convolved Flux"
        print new_image.sum()
    header = reference.header
    image = fits.ImageHDU(data=new_image, header=header)
    image.writeto(filename, clobber=True)
    # TODO: Implement these lines as part of energy binning
    #factor_1 = source_table['Flux10_30GeV'].sum()/im_1.sum()
    #im_1 =  gauss_correlate((im_1) * factor_1, 3)
    #image1 = fits.ImageHDU(data=im_1, header=header)
    #factor_2 = source_table['Flux30_100GeV'].sum()/im_2.sum()
    #im_2 =  gauss_correlate((im_2) * factor_2, 3)
    #image2 = fits.ImageHDU(data=im_2, header=header)
    #factor_3 = source_table['Flux100_500GeV'].sum()/im_3.sum()
    #im_3 =  gauss_correlate((im_3) * factor_3, 3)
    #image3 = fits.ImageHDU(data=im_3, header=header)
    # Also for energy binning
    # hdu_list = [image1, image2, image3]
    #cube_hdu = images_to_cube(hdu_list)
    #cube_hdu.writeto('1fhl_cube_gauss_2.fits', clobber=True) 

def catalog_to_table(catalog):
    """ TODO
    """
    from scipy.stats import gmean
    data = []
    if catalog == '1FHL':
        point_sources = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')
    elif catalog == '2FGL':
        point_sources = fetch_fermi_catalog(catalog, 'LAT_Point_Source_Catalog')
    n_sources = len(point_sources)
    sources = np.arange(n_sources)
    for source in sources:
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