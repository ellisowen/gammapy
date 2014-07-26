Fermi-LAT small data sets
=========================

We add a few small Fermi-LAT data files for the Vela region that we use for unit tests, examples and tutorials.

Parameters
----------

* 5 years of observation time (2008-08-05 to 2013-08-05)
* Event class and IRF: P7REP_CLEAN_V15
* Max zenith angle cut: 105 deg
* 50 MeV < Energy < 500000 MeV

Files
-----

* ``counts_vela.fits.gz`` -- Vela region counts cube 
* ``exposure_vela.fits.gz`` --	Vela region exposure cube


Details
-------

Commands:

I produced the `counts_vela.fits.gz` file using the commands in the executable script ``counts_commands.sh`` with the FSSC Fermi Science Tools.
These were also used to produce the `exposure_vela.fits.gz` file using the commands in the executable script ``exposure_commands.sh``.
   