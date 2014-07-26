Background Models & Spectral Cubes
==================================

The `gammapy.spectral_cube` module allows for image-based analysis in energy bands. In particular, similar functionality to gtmodel in the Fermi Science
tools [FSSC2014]_ is offered in `gammapy.spectral_cube.compute_npred_cube` which generates a predicted instrument PSF-convolved counts cube based on a
provided background model. Unlike the science tools, this implementation is appropriate for use with large regions of the sky. 


Predicting Counts
-----------------

The example script below computes the Fermi PSF-convolved predicted counts map using `gammapy.spectral_cube`.


.. plot:: tutorials/npred_convolved.py
	:include-source:
   
   
Checks
------

For small regions, the predicted counts cube and significance images may be checked against the gtmodel output. The Vela region is taken in this example in one
energy band with the following parameters:

The gtmodel result from 4 years of LAT data can be reproduced with the following commands... TODO...

From the predicted counts and observed counts, a Li & Ma significance image [LiMa1983] and excess image are produced.
 
TODO: include source

