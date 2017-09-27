# general python files
from os.path import join
import matplotlib.pyplot as pl
import numpy as np
pl.style.use('seaborn-white')

import pixeltoflux
import centroidfit
import pandas as pd

def makelc(inputpath='',outputpath='', find_transits=True,chunksize=300,cutoff_limit=1.1, showfig=None):#,campaign=1):
  # Takes strings with the EPIC number of the star and input/outputpath. Campaign number is used to complete the correct filename as downloaded from MAST
  # Set makelightcurve or find_transits to False to run only partial

  starname = inputpath.split('/')[-1].split('-')[0][4:-1]

  outputfolder = join(outputpath,str(starname))

  # makes raw light curve from pixel file0
  t,f_t,Xc,Yc = pixeltoflux.gotoflux(inputpath=inputpath,outputpath=outputpath,cutoff_limit=cutoff_limit, showfig=None)#,campaign=campaign)
  # removes outlying data points where thrusters are fired
  t,f_t,Xc,Yc = centroidfit.find_thruster_events(t,f_t,Xc,Yc,starname=starname,outputpath=outputfolder)

  if showfig:
    pl.show() # comment out to keep things running for multiple stars
  pl.close('all')

  return t,f_t,Xc,Yc

from astropy.stats import sigma_clip

def load_df(outputpath,starname,sigma_upper=2,sigma_lower=10):
    fname='lightcurve_raw_'+starname+'.txt'
    txtfile=join(outputpath,starname,fname)
    print(txtfile)
    df = pd.read_csv(txtfile, skiprows=1, delimiter=' ', names=['t','f','x','y'])
    df['f_clip'] = sigma_clip(df.f, sigma_upper=sigma_upper,sigma_lower=sigma_lower)
    df['fmed'] = df.f_clip.apply(lambda x: x/np.nanmedian(df.f_clip))
    df['f_mask'] = df.fmed[df.fmed.apply(lambda x: (x > 0.9) & (x < 1.1))]
    df.head()
    return df

def plot_f_xy(df,outputfolder=''):
    fig, ax = pl.subplots(2,1,figsize=(15,8))
    #ax[0].plot(df.t, df.f_clip,'.', alpha=0.5, color='r', marker='.')
    ax[0].plot(df.t, df.f_mask,'.', alpha=0.5, color='r', marker='.')
    ax[0].set_ylabel('Flux')
    ax[0].set_xlabel('Time (day)')
    #ax[0].set_ylim([np.median(df.fmed)-0.1,np.median(df.fmed)+0.1])
    ax[0].set_xlim([df.t.iloc[0],df.t.iloc[-1]])
    ax[0].set_title('sigma-clipped lc')
    #centroid drift
    ax[1].plot(df.t, df.x.apply(lambda x: x-df.x.iloc[0]),'.', alpha=0.5, color='g', marker='.', label='x')
    ax[1].plot(df.t, df.y.apply(lambda x: x-df.y.iloc[0]),'.', alpha=0.5, color='b', marker='.', label='y')
    ax[1].set_ylim([-1,1])
    ax[1].set_xlim([df.t.iloc[0],df.t.iloc[-1]])
    pl.legend()
    pl.savefig(join(outputfolder,'raw_lc.png'))

    print('clipped datapoints: {}'.format(np.isnan(df.f_mask).sum()))

def astropy_sigma_clip(t,f,x,y,sigma=3):
    fmask= sigma_clip(f, sigma=sigma).mask
    xmask= sigma_clip(x, sigma=sigma).mask
    ymask= sigma_clip(y, sigma=sigma).mask

    zz = np.c_[fmask, xmask, ymask]
    # idx = zz.sum(axis=1) != 0
    # idx.sum()
    mask=zz.any(axis=1)
    t,f,x,y=t[~mask],f[~mask],x[~mask],y[~mask]
    return t,f,x,y

from scipy.interpolate import Rbf

def fit_rbf(t,f, function='quintic', smooth=100,outputfolder='',showfig=True):
    '''
    A class for radial basis function approximation/interpolation of
    n-dimensional scattered data.
    :param function: The radial basis function, based on the radius, r, given by the norm
    :param smooth: Values greater than zero increase the smoothness of the
    approximation.

    100 is chosen by default such that short time-scale variability is not overfitted
    '''
    rbfi = Rbf(t, f, function='quintic', smooth=100)
    if showfig:
        fig, ax = pl.subplots(1,1,figsize=(15,3))
        ax.plot(t, f, '.')
        ax.plot(t, rbfi(t), 'r-')
        pl.savefig(join(outputfolder,'rbf_fit.png'))
    return t,rbfi(t)

def plot_rbf(t,f,f_rbf,outputfolder=''):
    fig, ax = pl.subplots(1,1,figsize=(15,3))
    f_cd = f / f_rbf
    pl.plot(t,f_cd,'.-', alpha=0.5)
    pl.savefig(join(outputfolder,'rbf_detrended.png'))
    return t, f_cd
