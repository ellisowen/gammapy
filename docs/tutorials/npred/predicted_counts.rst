Fermi-LAT diffuse model predicted counts image
==============================================

The `gammapy.spectral_cube` module allows for image-based analysis in energy
bands. In particular, similar functionality to gtmodel in the Fermi Science
tools [FSSC2013]_ is offered in `gammapy.spectral_cube.compute_npred_cube`
which generates a predicted instrument PSF-convolved counts cube based on a
provided background model. Unlike the science tools, this implementation is
appropriate for use with large regions of the sky. 


Predicting Counts
-----------------

The example script below computes the Fermi PSF-convolved predicted counts map
using `gammapy.spectral_cube`. This is then used to produce a Li & Ma significance
image [LiMa1983]_. The left image shows the significance image produced using the
methods in `gammapy.spectral_cube`, while a comparison against the significance image
produced using the Fermi Science tools is shown on the right. These results are
for the Vela region for energies between 10 and 500 GeV.

.. literalinclude:: npred_general.py

.. plot:: tutorials/npred/npred_convolved_significance.py
	:include-source:
   
   
Checks
------

For small regions, the predicted counts cube and significance images may be
checked against the gtmodel output. The Vela region shown above is taken in
this example in one energy band with the following parameters:

  * 5 years of observation time (2008-08-05 to 2013-08-05)
  * Event class and IRF: P7REP_CLEAN_V15
  * Max zenith angle cut: 105 deg
  * 10 GeV < Energy < 500 GeV
  * Position: GLON = 263.05836967702709; GLAT = 3.9298511274632784
  * Image Size: 5 x 5 degrees (Search center 3 degrees)
  * Pixel size: 0.1 degrees/pixel
  * Image Shape: (50, 50) pixels

Images for the predicted background counts in this region in the Gammapy case
(left) and Fermi Science Tools gtmodel case (center) are shown below, based on
the differential flux contribution of the Fermi diffuse model gll_iem_v05_rev1.fit.
The image on the right shows the ratio.

.. plot:: tutorials/npred/npred_convolved.py
	:include-source:

We may compare these against the true counts observed by Fermi LAT in this region
for the same parameters:

 * True counts: 261
 * Fermi Tools gtmodel predicted counts: 145
 * Gammapy predicted counts: 271
