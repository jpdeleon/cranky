# k2photometry 2
This code was forked from [V. Van Eylen](https://github.com/vincentvaneylen/k2photometry) used to read, reduce and detrend K2 photometry. 

Changes: 

Edited `pixeltoflux()` and `run_pipeline()`
* simplified get_pixelfluxes() & edited gotoflux()
* removed `campaign` variable
* created example2.py using new TPF
* uses astropy.io.fits instead of pyfits

To do:
* improve sigma clipping and centroiding

Additional requirements:
* [bls](https://github.com/dfm/python-bls)
* [lmfit](https://github.com/lmfit/lmfit-py/)

