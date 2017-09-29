import numpy as np
import matplotlib.pyplot as pl
import gatspy
from os.path import exists, join

from gatspy.periodic import LombScargleFast

from periodfinder import *
from auxiliaries import *

import bls


def plot_raw_lc(t,f,outputfolder='',starname=''):
    fig, ax = pl.subplots(1,1,figsize=(15,3))
    ax.plot(t, f, '.')
    ax.set_ylabel('Flux')
    ax.set_xlabel('Time (day)')
    ax.set_title('working lc for planet finding: {}'.format(starname))
    pl.savefig(join(outputfolder, 'working_lc_' + str(starname) + '.png'))
    return None


def find_period(t,f,showfig=True,outputfolder='',starname=''):
    '''
    uses Lomb-Scargle to find planet's orbital period
    '''
    model = LombScargleFast().fit(t, f)
    periods, power = model.periodogram_auto(nyquist_factor=100)
    idx1 = periods > 1
    if showfig:
        idx2 = np.argmax(power[idx1])
        period = periods[idx1][idx2]

        fig, ax = pl.subplots(1,1,figsize=(15,5))
        ax.plot(periods, power, 'k-')
        ax.set(xlim=(0.5, 5),
        #         , ylim=(0, 0.01),
           xlabel='period (days)',
           ylabel='Lomb-Scargle Power')
        ax.vlines(period, *ax.get_ylim(), linestyles='dotted', colors='r')
        ax.set_title('best period: {0:.3f}'.format(period))
    pl.savefig(join(outputfolder, 'LombScargle_' + str(starname) + '.png'))
    return period, power

def estimate_k(t,f,p,showfig=False):
    '''
    k=Rp/Rs estimate transit depth assumed to be within 0.01 percentile
    '''
    baseline,minimum=np.percentile(f[np.array(t%p).argsort()], [50,1])
    if showfig:
        pl.hist(f[np.array(t%p).argsort()]);
    #print(baseline,minimum)
    return np.sqrt(baseline-minimum)


def estimate_t0(t,f,tlower=None,tupper=None, showfig=False,outputfolder='',starname=''):
    '''
    very basic search for t0 given interval tlower and tupper
    '''
    if tlower is None and tupper is None:
        tlower, tupper=t[0],t[200]
    idx = (t > tlower) & (t < tupper)
    tsub, fsub = t[idx], f[idx]
    idx = fsub < np.median(fsub) - 0.5 * np.nanstd(fsub)
    t0 = np.median(tsub[idx])
    if showfig:
        fig, ax = pl.subplots(1,1,figsize=(15,3))
        ax.plot(tsub, fsub, '.')
        ax.vlines(t0, *ax.get_ylim())
        #ax.set_title(starname)
        ax.set_title('t0= {0:.3f}'.format(t0))
    pl.savefig(join(outputfolder, str(starname) +'_t0.png'))
    return t0


def get_tns(t, p, t0):
    '''
    determine transit conjuctions
    '''

    idx = t != 0
    t = t[idx]

    while t0-p > t.min():
        t0 -= p
    if t0 < t.min():
        t0 += p

    tns = [t0+p*i for i in range(int((t.max()-t0)/p+1))]

    while tns[-1] > t.max():
        tns.pop()

    while tns[0] < t.min():
        tns = tns[1:]

    return tns


def fold(t, f, p, t0, width=0.4, clip=False, bl=False, t14=0.1):
    '''
    fold lc; very sensitive to p & t0
    '''
    tns = get_tns(t, p, t0)
    tf, ff = np.empty(0), np.empty(0)
    for i,tn in enumerate(tns):
        idx = (t > tn - width/2.) & (t < tn + width/2.)
        ti = t[idx]-tn
        fi = f[idx]
        fi /= np.nanmedian(fi)
        tf = np.append(tf, ti)
        ff = np.append(ff, fi / np.nanmedian(fi))
    idx = np.argsort(tf)
    tf = tf[idx]
    ff = ff[idx]
    return tf, ff

def fold_with_eclipse(t, f, p, t0, width=0.4, t14=0.1, odd_even=False):
    '''
    fold lc; very sensitive to p & t0
    '''
    tns = get_tns(t, p, t0)
    tf, ff = np.empty(0), np.empty(0)

    #for odd_even
    tf1, ff1 = np.empty(0), np.empty(0)
    tf2, ff2 = np.empty(0), np.empty(0)

    if odd_even:
        for i,tn in enumerate(tns[1::2]): #odd
            idx_odd = (t > tn - width/2.) & (t < tn + width/2.)
            ti1 = t[idx_odd]-tn
            fi1 = f[idx_odd]
            fi1 /= np.nanmedian(fi1)
            tf1 = np.append(tf1, ti1)
            ff1 = np.append(ff1, fi1 / np.nanmedian(fi1))

        idx_odd = np.argsort(tf1)
        tf1 = tf1[idx_odd]
        ff1 = ff1[idx_odd]

        for i,tn in enumerate(tns[::2]): #even
            idx_even = (t > tn - width/2.) & (t < tn + width/2.)
            ti2 = t[idx_even]-tn
            fi2 = f[idx_even]
            fi2 /= np.nanmedian(fi2)
            tf2 = np.append(tf2, ti2)
            ff2 = np.append(ff2, fi2 / np.nanmedian(fi2))

        idx_even = np.argsort(tf2)
        tf2 = tf2[idx_even]
        ff2 = ff2[idx_even]

        return tf1, ff1, tf2,ff2
    else:
        for i,tn in enumerate(tns):
            idx = (t > tn - width/2.) & (t < tn + width/2.)
            ti = t[idx]-tn
            fi = f[idx]
            fi /= np.nanmedian(fi)
            tf = np.append(tf, ti)
            ff = np.append(ff, fi / np.nanmedian(fi))

        idx = np.argsort(tf)
        tf = tf[idx]
        ff = ff[idx]
        return tf, ff

def fold_data(t,y,period):
    '''
    alternative to fold function above
    '''
    # simple module to fold data based on period
    folded = t % period
    inds = np.array(folded).argsort()
    t_folded = folded[inds]
    y_folded = y[inds]

    return t_folded,y_folded


def plot_tns(tns,t,f, f1=0.975, f2=0.98,outputfolder='',starname=''):
    fig, ax = pl.subplots(1,1,figsize=(15,5))
    ax.plot(t, f, '.')
    ax.vlines(tns, f1,f2)
    ax.set_title(starname)
    pl.savefig(join(outputfolder, 'tns_' + str(starname) + '.png'))
    return None

def plot_folded_lc(t,f,period,t0,outputfolder='',starname=''):
    tf, ff = fold(t, f, period, t0)

    fig, ax = pl.subplots(1,1,figsize=(15,5))
    ax.plot(tf, ff, '.')
    ax.set_title('Phase-folded lc before optimization'.format(starname))
    ax.set_xlabel('Phase')
    ax.set_ylabel('Normalized Flux')
    pl.savefig(join(outputfolder, 'folded_lc_' + str(starname) + '.png'))
    return None

def estimate_t14(value=0.01):
    if value:
        return value


def t14_circ(p, a, k, b):
    """
    Winn 2014 ("Transits and Occultations"), eq. 14
    """
    i = inclination(a, b)
    alpha = np.sqrt( (1 + k)**2 - b**2 )
    return (p / np.pi) * np.arcsin( alpha / np.sin(i) / a )


def tshape_approx(a, k, b):
    """
    Seager & Mallen-Ornelas 2003, eq. 15
    """
    i = inclination(a, b)
    alpha = (1 - k)**2 - b**2
    beta = (1 + k)**2 - b**2
    return np.sqrt( alpha / beta )


def max_k(tshape):
    """
    Seager & Mallen-Ornelas 2003, eq. 21
    """
    return (1 - tshape) / (1 + tshape)


def scaled_a(p, t14, k, i=np.pi/2, b=0):
    numer = np.sqrt( (k + 1)**2 - b**2 )
    denom = np.sin(i) * np.sin(t14 * np.pi / p)
    return float(numer / denom)


def compute_a(a_scaled,Rstar):
    return Rstar/a_scaled


######################################
from pytransit import MandelAgol

def model_q(theta, t):
    MA = MandelAgol()
    k,tc,p,a,b,q1,q2,_ = theta
    i     = inclination(a, b)
    u1,u2 = q_to_u(q1, q2)
    model = MA.evaluate(t, k, (u1,u2), tc, p, a, i)

    return model


def lnlike(theta, t, f):
    k,t0,p,a,b,q1,q2,sig = theta
    m = model_q(theta, t)
    resid = f - m
    inv_sigma2 = 1.0/(sig**2)

    return -0.5*(np.sum((resid)**2*inv_sigma2 - np.log(inv_sigma2)))

from scipy import stats

def lnprob(theta, t, f):
    k,t0,p,a,b,q1,q2,sig = theta
    #logprior
    # -1<t0<1 if t,f inputs are folded
    # if t0 < -1 or t0 > 1:
    #     return -np.inf
    if q1 < 0 or q1 > 1 or q2 < 0 or q2 > 1 or b < 0 or b > 1 or k < 0 or k > 1:
        return -np.inf

    u1, u2 = q_to_u(q1, q2)

    lp = 0
    # if up is not None:
    #     lp += np.log(stats.norm.pdf(u1, loc=up[0], scale=up[1]))
    #     lp += np.log(stats.norm.pdf(u2, loc=up[2], scale=up[3]))

    #loglike
    ll = lnlike(theta, t, f)

    if np.isnan(ll).any():
        return -np.inf
    return lp + ll


import scipy.optimize as op

def fit_folded_lc(initial, args, method='powell',verbose=True):
    '''
    `Powell` method minimises the function by a bi-directional search along each search vector, in turn
    '''
    nlp = lambda *args: -lnprob(*args)

    opt = op.minimize(nlp, initial, args=args, method=method)
    labels='k,t0,p,a,b,q1,q2,sig'.split(',')

    if verbose:
        print('converged: {}'.format(opt.success))
        for i,j in zip(labels,opt.x):
            print('{0}: {1}'.format(i,j))
    return opt


def plot_fit(theta,t,f,show_model=True,outputfolder='',starname='', folded=False):
    t0, p = theta[1], theta[2]
    #f /= np.median(f)
    fig, ax = pl.subplots(1,1,figsize=(15,5))

    if folded:
        #t,f and t0 are folded
        ax.plot(t, f, '.', label='data')
        #note that t0 in theta is t0_fold
        fmod=model_q(theta, t)
        ax.plot(t,fmod,'r.-',label='model')
        residual = f - fmod

    else:
        #t,f must be folded
        tf, ff = fold(t, f, p, t0)
        ax.plot(tf, ff, '.', label='data')
        fmod=model_q(theta, t)
        ttmod,ffmod=fold(t,fmod,p,t0)
        residual = ff - ffmod
        ax.plot(ttmod,ffmod,'r.-',label='model')

    rms = np.sqrt(np.mean(residual**2))
    ax.set_xlabel('Phase')
    ax.set_ylabel('Normalized Flux')
    ax.set_title('Phase-folded with model fit (MLE): {}'.format(starname))
    ax.legend()
    pl.savefig(join(outputfolder, 'folded_model_optimized_fit' + str(starname) + '.png'))
    return rms


def model_s(theta, x, y, t):
    '''
    theta contains auxiliary coefficients
    linear combination of auxiliary vectors
    '''
    offset = np.ones(len(t))
    s = (np.array(theta)*np.c_[x, y, x*y, x**2, y**2, offset, t]).sum(axis=1)
    return s

def least_squares_model(f,x,y):
    '''
    Least-squares estimation to solve for model weights
    y*X.T=w
    -> w=inv(X*X.T)*(y*X.T)
    '''
    #design matrix
    X = np.c_[x,y,x*y,x**2,y**2,np.ones_like(x)]
    w = np.dot(np.dot(np.linalg.inv(np.dot(X.T,X)),X.T),f)
    model = np.dot(X,w)

    return model

def find_period_combo(t,f,showfig=True,outputfolder='',starname=''):
    '''
    uses Lomb-Scargle to find planet's orbital period
    '''
    model = LombScargleFast().fit(t, f)
    periods, power = model.periodogram_auto(nyquist_factor=100)
    idx1 = periods > 1
    idx2 = np.argmax(power[idx1])
    period = periods[idx1][idx2]

    if showfig:
        fig, ax = pl.subplots(1,1,figsize=(15,5))
        ax.plot(periods, power, 'k-')
        ax.set(xlim=(0.5, 5),
        #         , ylim=(0, 0.01),
           xlabel='period (days)',
           ylabel='Lomb-Scargle Power')
        ax.vlines(period, *ax.get_ylim(), linestyles='dotted', colors='r')
        ax.set_title('best period: {0:.3f}'.format(period))
    return period, periods, power

def summarize(initial,args,freqlist=None,powers=None,minimum=3,outputfolder='',starname=''):
    '''
    Create a summary consisting of the ff plots:
    1. raw lc
    2. LS periodogram
    3. BLS periodogram
    4. folded lc for primary transit with model
    5. folded lc for secondary eclipse to check for FP
    6. folded on odd periods  to check for FP
    7. folded on even periods to check for FP

    :param minimum: percentile of flux to compute for minimum transit depth
    '''
    import matplotlib as mpl
    mpl.rcParams['font.size'] = 20

    fig = pl.figure('Quick Look',figsize=(35.,25.))

    rows,cols= 8,2
    ax1 = pl.subplot2grid((rows,cols), (0,0), colspan=3,rowspan=2)
    ax2 = pl.subplot2grid((rows,cols), (2,0), colspan=3)
    ax3 = pl.subplot2grid((rows,cols), (3,0), colspan=3)
    ax4 = pl.subplot2grid((rows,cols), (4,0), rowspan=2)
    ax5 = pl.subplot2grid((rows,cols), (4,1), rowspan=2)
    ax6 = pl.subplot2grid((rows,cols), (6,0), rowspan=2)
    ax7 = pl.subplot2grid((rows,cols), (6,1), rowspan=2)

    #for modeling

    k,t0,p,a_au,b,q1,q2,sig = initial
    t, f = args

    tt,ff = fold(t,f,p,t0)
    #adjust t0 to search near phase = 0
    t0_fold=estimate_t0(tt,ff, tlower=tt[len(tt)/2-30], tupper=tt[len(tt)/2+30],
                        outputfolder=outputfolder,starname=starname, showfig=True)
    initial_fold = k,t0_fold,p,a_au,b,q1,q2,sig
    args_fold = (tt, ff)

    #--------------------raw lightcurve--------------------#
    print('raw light curve done')
    ax1.plot(t, f, '.')
    ax1.set_ylabel('Flux')
    ax1.set_xlabel('Time (day)')
    ax1.set_title('working lc ({})'.format(starname))

    #--------------------BLS periodogram--------------------#
    print('bls periodogram done')
    if freqlist is None or powers is None: #use_BLS:
        folded,t_folded,period,freqlist,powers = get_period(t,f,[],
                                                    #showfig=True,
                                                    starname=starname,
                                                    get_mandelagolmodel=False)
    else:
        pass
    ax2.plot(1./freqlist, powers, 'k-')
    ax2.set(xlim=(0.5, 5), #, ylim=(0, 0.01),
            xlabel='Period (days)',
            ylabel='BLS Power')
    snr_bls = np.max(powers)/np.median(powers)
    ax2.legend(['period={0:.4f} d\nsnr={1:.3f}'.format(period,snr_bls)],loc=1)
    ax2.vlines(period, *ax2.get_ylim(), linestyles='dotted', colors='r')
    #ax2.set_xticks([])
    #--------------------LS periodogram--------------------#
    print('lomb-scargle periodogram done')
    period, periods, powers = find_period_combo(t,f,showfig=False)

    ax3.plot(periods, powers, 'k-')
    ax3.set(xlim=(0.5, 5), #, ylim=(0, 0.01),
            xlabel='Period (days)',
            ylabel='LS Power')
    snr_ls = np.max(powers)/np.median(powers)
    ax3.legend(['period={0:.4f} d\nsnr={1:.3f}'.format(period,snr_ls)],loc=1)
    ax3.vlines(period, *ax3.get_ylim(), linestyles='dotted', colors='r')
    #ax3.set_xticks([])

    #--------------------model fit--------------------#
    print('primary transit model fitting done')

    method='powell'
    opt=fit_folded_lc(initial, args, method=method, verbose=False)
    t0, p = opt.x[1], opt.x[2]

    tf, ff = fold(t, f, p, t0)

    ax4.plot(tf, ff, '.', label='data')
    ax4.set_xlabel('Phase')
    ax4.set_ylabel('Normalized Flux')
    fmod=model_q(opt.x, t)
    ttmod,ffmod=fold(t,fmod,p,t0)
    ax4.plot(ttmod,ffmod,'r-',label='model')
    ax4.set_title('Phase-folded on primary transit with model fit (MLE)')
    #ax4.set_xticks([])
    ax4.legend()

    #--------------------secondary eclipse--------------------#
    print('secondary eclipse model fitting done')
    tf2,ff2=fold(t, f, period, t0+period/2, width=1, t14=0.1)

    ax5.plot(tf2, ff2, '.', label='data')
    ax5.set_title('Phase-folded on secondary eclipse')# with model fit (MLE)')
    ax5.set_xlabel('Phase')
    ax5.set_ylabel('Normalized Flux')
    ax5.set_ylim(*ax4.get_ylim())
    #ax5.set_xticks([])

    #--------------------odd--------------------#
    print('odd period phase-folding done')

    tf_odd,ff_odd, tf_even,ff_even = fold_with_eclipse(t, f, p, t0, width=0.4, t14=0.1, odd_even=True)
    #MLE fit
    opt_odd  = fit_folded_lc(initial_fold, args=(tf_odd,ff_odd), method=method, verbose=False)
    opt_even = fit_folded_lc(initial_fold, args=(tf_even,ff_even), method=method, verbose=False)

    ax6.plot(tf_odd, ff_odd, '.', label='data')
    #model
    fmod_odd=model_q(opt_odd.x, tf_odd)
    ax6.plot(tf_odd,fmod_odd,'r-',label='model')
    depth_odd=np.percentile(ff_odd, minimum) #np.min(ff_odd)
    ax6.hlines(depth_odd, *ax6.set_xlim(), linestyles='dashed', colors='k')

    ax6.set_title('Phase-folded on ODD periods with model fit (MLE)')
    ax6.set_xlabel('Phase')
    ax6.set_ylabel('Normalized Flux')

    #--------------------even--------------------#
    print('even period phase-folding done')
    ax7.plot(tf_even, ff_even, '.', label='data')
    #model
    fmod_even=model_q(opt_even.x, tf_even)
    ax7.plot(tf_even,fmod_even,'r-',label='model')
    depth_even=np.percentile(ff_even, minimum) #np.min(ff_even)
    ax7.hlines(depth_even, *ax7.set_xlim(), linestyles='dashed', colors='k')

    ax7.set_title('Phase-folded on EVEN periods with model fit (MLE)')
    ax7.set_xlabel('Phase')
    ax7.set_ylabel('Normalized Flux')
    ax7.set_ylim(*ax6.get_ylim())

    fig.tight_layout()
    fig.savefig(join(outputfolder, 'quick_look_' + str(starname) + '.png'))

def baseline(theta, t):
    ti = t - t.mean()
    c0,c1,c2,c3 = theta[-4:]
    return c0 + c1 * ti + c2 * ti**2 + c3 * ti**3
