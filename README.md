# k2photometry 2
======
This code was forked from [V. Van Eylen](https://github.com/vincentvaneylen/k2photometry) used to read, reduce and detrend K2 photometry. 

## Installation
======
Clone the source code from the [GitHub repository](https://github.com/jpdeleon/k2photometry) and install using the standard python tools:

```bash
    git clone https://github.com/dfm/python-bls.git
    cd python-bls
    python setup.py install
```

## Quick test
======
To run, specify at least the path to the (.fits) data directory, `indir`. Set `showfig` to `False` to quickly loop all data.

```bash
$ cranker --indir 'data/' --showfig False
```

## Notes
====== 
09/12
* created `cranker` script
* added setup.py

09/11
Edited `pixeltoflux()` and `run_pipeline()`
* simplified get_pixelfluxes() & edited gotoflux()
* removed `campaign` variable
* created example2.py using new TPF
* uses astropy.io.fits instead of pyfits

TODO:
* improve sigma clipping and centroiding
* add basic license

Additional requirements:
* [bls](https://github.com/dfm/python-bls)
* [lmfit](https://github.com/lmfit/lmfit-py/)

