Fermi-LAT small data sets
=========================

We add a few small Fermi-LAT data files that we use for unit tests, examples and tutorials.

Parameters
----------

* 5 years of observation time (2008-08-05 to 2013-08-05)
* Event class and IRF: P7REP_CLEAN_V15
* Max zenith angle cut: 105 deg
* 10 GeV < Energy < 500 GeV

Files
-----

* ``fermi_counts.fits.gz`` -- Galactic center counts image
* ``psf.fits`` -- Galactic center PSF
* ``gll_iem_v02_cutout.fits`` -- Galactic center diffuse model cube
* ``fermi_exposure.fits.gz`` -- Galactic center exposure cube


Details
-------

Diffuse Model Cube: `gll_iem_v02_cutout.fits`

Fermi LAT data server query parameters:

* Equatorial coordinates (degrees) (266.405,-28.9362)
* Time range (MET)  (239587200,397353600)
* Time range (Gregorian)  (2008-08-05 00:00:00,2013-08-05 00:00:00)
* Energy range (MeV)   (10000,500000)
* Search radius (degrees) 30

Commands:

I produced the `gll_iem_v02_cutout.fits` file using these commands::

   $ wget http://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/gll_iem_v02.fit
   $ ftcopy 'gll_iem_v02.fit[330:390,170:190,*]' gll_iem_v02_cutout.fits
   $ fchecksum gll_iem_v02_cutout.fits update+ datasum+

   
Exposure Cube: `fermi_exposure.fits.gz`

Fermi LAT Key Parameters:

* Equatorial coordinates (degrees) (266.405,-28.9362)
* Time range (MET)  (239587200,397353600)
* Time range (Gregorian)  (2008-08-05 00:00:00,2013-08-05 00:00:00)
* Energy range (MeV)   (10000,1000000)
* Search radius (degrees) 30

Commands:

I produced the `fermi_exposure.fits.gz` file using these commands from the FSSC Fermi Science Tools::

   $ gtselect
   Input FT1 file[events.txt] 
   Output FT1 file[events.fits] 
   RA for new search center (degrees) (0:360) [266.404947172699] 
   Dec for new search center (degrees) (-90:90) [-28.9362422432238] 
   radius of new search region (degrees) (0:180) [30]
   start time (MET in s) (0:) [239587200] 
   end time (MET in s) (0:) [397353600] 
   lower energy limit (MeV) (0:) [50] 
   upper energy limit (MeV) (0:) [1000000] 
   maximum zenith angle value (degrees) (0:180) [105] 
   $ gtbin
   This is gtbin version ScienceTools-v9r32p5-fssc-20130916
   Type of output file (CCUBE|CMAP|LC|PHA1|PHA2|HEALPIX) [CCUBE] 
   Event data file name[events.fits] 
   Output file name[counts.fits] 
   Spacecraft data file name[spacecraft.fits] 
   Size of the X axis in pixels[61] 
   Size of the Y axis in pixels[21] 
   Image scale (in degrees/pixel)[1] 
   Coordinate system (CEL - celestial, GAL -galactic) (CEL|GAL) [GAL]
   First coordinate of image center in degrees (RA or galactic l)[0] 
   Second coordinate of image center in degrees (DEC or galactic b)[0] 
   Rotation angle of image axis, in degrees[0] 
   Projection method e.g. AIT|ARC|CAR|GLS|MER|NCP|SIN|STG|TAN:[CAR] 
   Algorithm for defining energy bins (FILE|LIN|LOG) [LOG] 
   Start value for first energy bin in MeV[50] 
   Stop value for last energy bin in MeV[1000000]
   Number of logarithmically uniform energy bins[60] 
   $ gtltcube
   Event data file[events.fits] 
   Spacecraft data file[spacecraft.fits] 
   Output file[livetime.fits] 
   Step size in cos(theta) (0.:1.) [0.1] 
   Pixel size (degrees)[1]
   $ gtexpcube2
   Livetime cube file[livetime.fits] 
   Counts map file[counts.fits] 
   Output file name[exposure_cube.fits] 
   Response functions to use[P7REP_CLEAN_V15] 
   