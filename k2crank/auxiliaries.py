# general python files
import numpy as np

def inclination(a, b, e=None, w=None):
    """
    Winn 2014 ("Transits and Occultations"), eq. 7
    """
    # if e is None and w is None:
    #     return np.arccos(b / a)
    # elif e is not None and w is not None:
    #     return np.arccos(b / a * (1 + e * np.sin(w)) / (1 - e**2))
    # else:
    #     return np.nan
    return np.arccos(b / a)

def q_to_u(q1, q2):
    u1 = 2 * np.sqrt(q1) * q2
    u2 = np.sqrt(q1) * (1 - 2*q2)
    return u1, u2


def u_to_q(u1, u2):
    q1 = (u1 + u2)**2
    q2 = u1 / (2 * (u1 + u2))
    return q1, q2


def impact_param(a,Rstar):
    '''
    a is semi-major axis, not scaled_a = R_star/a
    '''
    return a * np.cos(i)/Rstar

def rebin_dataset(x,y,bin_size,medianbin=False):   ## FIXME calculate proper stdv in case of medianbinning
  #
  # Auxiliary function to do a binning, either mean or median
  #
  inds = np.array(x).argsort()
  x = x[inds]
  y = y[inds]


  x_rebinned = []
  y_rebinned = []
  stdv_bins = []
  i = 0
  while i < len(y):
    j = 0
    average = 0
    x_rebinned.append(x[i])
    bin_values = []
    if not medianbin:
      while (j < bin_size) and (i < len(y)):
	bin_values.append(y[i])
	average += y[i]
	j = j + 1
	i = i + 1
      average = average/j
      stdv=average
    # find stdv
    else:
      while (j < bin_size) and (i < len(y)):
	bin_values.append(y[i])
	j = j + 1
	i = i + 1
      average = np.median(bin_values)
      mad = 1.4826 * np.median( np.abs(average - np.array(bin_values)) ) / ((len(bin_values) -1)**0.5)
      #stdv = np.std(bin_values) / (len(bin_values)**0.5)
      stdv = mad

    stdv_bins.append(stdv)
    y_rebinned.append(average)
  return [x_rebinned,y_rebinned,stdv_bins]


def running_sigma_clip(data,sigma=3,binsize=10,dependent_var=None):
  #
  # Sigma clipping (running): find local outliers
  #
  if dependent_var is not None:
    dep_var_clipped = []
  data_clipped = []
  upperlist = []
  lowerlist = []
  i = 0
  while i < len(data):
    bin_begin = max(0, (i - binsize/2))
    bin_end = min(len(data),(i+binsize/2))
    the_bin = data[bin_begin:bin_end]

    std = np.nanstd(np.sort(the_bin)[1:])
    median = np.median(the_bin)
    upperbound = (median + (sigma*std))
    lowerbound = (median - (sigma*std))
    upperlist.append(upperbound)
    lowerlist.append(lowerbound)
    if (data[i] < upperbound) and (data[i] > lowerbound):
      data_clipped.append(data[i])
      if dependent_var is not None:
	dep_var_clipped.append(dependent_var[i])

    i = i + 1
  if dependent_var is not None:
    return [data_clipped,dep_var_clipped,upperlist,lowerlist]
  else:
    return data_clipped


def sigma_clip(data,sigma=3,dependent_var=None,iterative=False,top_only=False):
  #
  # Auxiliary to sigma clip, option to run more than once, and to include only top data points (lower outliers may be transit events so we have to be careful what we clip)
  #

  repeat=True
  while repeat:
    # sigma clip data for outliers beyond sigma
    data = np.array(data)
    mean = np.nanmean(data)
    std = np.nanstd(data)
    if top_only:
      unclipped = (data < (mean + sigma*std))
    else:
      unclipped = (data < (mean + sigma*std))*(data > (mean - sigma*std)) # array of true and false

    data = data[unclipped]

    if dependent_var is not None:
      dependent_var = np.array(dependent_var)
      dependent_var = dependent_var[unclipped]

    if (iterative and (np.sum(~unclipped) > 0)):
      #	print 'Data points removed:'
      #print np.sum(~unclipped)
      #print 'Repeating..'
      repeat = True

    else:
      #print 'End of sigma clipping, data points removed:'
      #print np.sum(~unclipped)
      repeat = False


  if dependent_var is not None:
    return [data,dependent_var]
  else:
    return data


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

import sys

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")
