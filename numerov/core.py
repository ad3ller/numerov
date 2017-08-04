# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 16:20:34 2017

@author: Adam
"""
from math import exp, ceil, log
import numpy as np

def wf_numerov(n, l, nmax, step=0.005, rmin=0.65):
    """ Use the Numerov method to find the wavefunction for state n*, l, where
        n* = n - delta.
        nmax ensures that wavefunctions from different values of n can be aligned.
    """
    w1 = -0.5 * n**-2.0
    w2 = (l + 0.5)**2.0
    rmax = 2 * nmax * (nmax + 15)
    r_in = n**2.0 - n * (n**2.0 - l*(l + 1.0))**0.5
    step_sq = step**2.0
    # ensure wf arrays will align using nmax
    if n == nmax:
        i = 0
        r_sub2 = rmax
    else:
        i = int(ceil(log(rmax / (2 * n * (n + 15))) / step))
        r_sub2 = rmax * exp(-i*step)
    i += 1

    # initialise
    r_sub1 = rmax * exp(-i*step)
    rvals = [r_sub2, r_sub1]
    g_sub2 = 2.0 * r_sub2**2.0 * (-1.0 / r_sub2 - w1) + w2
    g_sub1 = 2.0 * r_sub1**2.0 * (-1.0 / r_sub1 - w1) + w2
    y_sub2 = 1e-10
    y_sub1 = y_sub2 * (1.0 + step * g_sub2**0.5)
    yvals = [y_sub2, y_sub1]

    # Numerov method
    i += 1
    r = r_sub1
    drr = exp(-step)**(-l - 1) - 1.0
    while r >= rmin:
        ## next step
        r = rmax * exp(-i*step)
        g = 2.0 * r**2.0 * (-1.0 / r - w1) + w2
        y = (y_sub2 * (g_sub2 - (12.0 / step_sq)) + y_sub1 * \
            (10.0 * g_sub1 + (24.0 / step_sq))) / ((12.0 / step_sq) - g)

        ## check for divergence
        if r < r_in:
            dyy = abs((y - y_sub1) / y_sub1)
            if dyy > drr:
                break

        ## store vals
        rvals.append(r)
        yvals.append(y)

        ## next iteration
        g_sub2 = g_sub1
        g_sub1 = g
        y_sub2 = y_sub1
        y_sub1 = y
        i += 1

    rvals = np.array(rvals)
    yvals = np.array(yvals)
    # normalisation
    yvals = yvals * (rvals * np.sum(yvals**2.0 * rvals**2.0))**-0.5 / step**0.5
    return rvals, yvals
