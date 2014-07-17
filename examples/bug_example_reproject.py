from FITS_tools.cube_regrid import regrid_cube_hdu

diffuse = fits.open(FermiGalacticCenter.filenames()['diffuse_model'])[0]
exposure = fits.open(FermiGalacticCenter.filenames()['exposure_cube'])[0]

regrid_cube_hdu(hdu1, hdu2.header, smooth=True)