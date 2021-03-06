# cranky

This code was forked from [V. Van Eylen](https://github.com/vincentvaneylen/k2photometry). This pipeline is used to read, reduce and detrend K2 photometry. 

## Installation

Clone the source code from the [GitHub repository](https://github.com/jpdeleon/cranky) and install using the standard python tools:

```bash
git clone https://github.com/jpdeleon/cranky.git
cd cranky
python setup.py install
```

## Quick test

To run, specify at least the path to the data directory, `-i`. `showfig` is a flag.

For simple sigma clipping and detrending, 

```bash
$ step1 -i '/path/to/data/*.fits.gz' -showfig
```

Outputs are saved in `output1` folder by default unless specified in `-o`.


`step2` does simple phase-folding and parameter inference using MLE. The input comes from `step1` output: `*detrended_lc_<starname>.txt` which containes (time,flux,x,y).

```bash
$ step2 -i '/output/starname/*detrended_lc_<starname>.txt'
```

## Notes

09/28
* created summary/ quick-look figure
 
09/25
* created scripts:
  * `step1` for simple sigma clipping and detrending
  * `step2` for phase-folding, and parameter inference via MLE
 
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

