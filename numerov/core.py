#! python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 16:20:34 2017

@author: Adam
"""
from math import exp
import numpy as np

def radial_wf(n, l, step=0.005, rmin=0.65):
    """ Use the Numerov method to find R_nl """
    w1 = -0.5 * n**-2.0
    w2 = (l + 0.5)**2.0
    rmax = 2 * n * (n + 15)
    r_in = n**2.0 - n * (n**2.0 - l*(l + 1.0))**0.5
    step_sq = step**2.0

    # initialise
    i = 1
    r_sub2 = rmax
    r = r_sub1 = rmax * exp(-i*step)
    rvals = [r_sub2, r_sub1]
    g_sub2 = 2.0 * r_sub2**2.0 * (-1.0 / r_sub2 - w1) + w2
    g_sub1 = 2.0 * r_sub1**2.0 * (-1.0 / r_sub1 - w1) + w2
    y_sub2 = 1e-10
    y_sub1 = y_sub2 * (1.0 + step * g_sub2**0.5)
    yvals = [y_sub2, y_sub1]

    # Numerov method
    i += 1
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

def radial_integral(n1, l1, n2, l2, step=0.005, rmin=0.65, p=1):
    """ Use the Numerov method to calculate:

        integrate(R_nl(n1, l1) * r^(2 + p) * R_nl(n2, l2), (r, 0, oo))
    """
    w11 = -0.5 * n1**-2.0
    w12 = (l1 + 0.5)**2.0
    w21 = -0.5 * n2**-2.0
    w22 = (l2 + 0.5)**2.0
    nmax = max(n1, n2)
    rmax = 2 * nmax * (nmax + 15)
    r_in1 = n1**2.0 - n1 * (n1**2.0 - l1*(l1 + 1.0))**0.5
    r_in2 = n2**2.0 - n2 * (n2**2.0 - l2*(l2 + 1.0))**0.5
    step_sq = step**2.0

    # initialise
    r_sub2 = rmax
    r_sub1 = rmax * exp(-step)

    g1_sub2 = 2.0 * r_sub2**2.0 * (-1.0 / r_sub2 - w11) + w12
    g1_sub1 = 2.0 * r_sub1**2.0 * (-1.0 / r_sub1 - w11) + w12
    g2_sub2 = 2.0 * r_sub2**2.0 * (-1.0 / r_sub2 - w21) + w22
    g2_sub1 = 2.0 * r_sub1**2.0 * (-1.0 / r_sub1 - w21) + w22

    y1_sub2 = 1e-10
    y1_sub1 = y1_sub2 * (1.0 + step * g1_sub2**0.5)
    y2_sub2 = 1e-10
    y2_sub1 = y2_sub2 * (1.0 + step * g2_sub2**0.5)

    norm1 = (y1_sub2**2.0 * r_sub2**2.0
             + y1_sub1**2.0 * r_sub1**2.0)
    norm2 = (y2_sub2**2.0 * r_sub2**2.0
             + y2_sub1**2.0 * r_sub1**2.0)

    integral = (y1_sub2 * y2_sub2 * r_sub2**(2.0 + p)
                + y1_sub1 * y2_sub1 * r_sub1**(2.0 + p))

    i = 2
    r = r_sub1
    dr1 = exp(-step)**(-l1 - 1) - 1.0
    dr2 = exp(-step)**(-l2 - 1) - 1.0
    while r >= rmin:
        # Numerov method
        r = rmax * exp(-i*step)
        g1 = 2.0 * r**2.0 * (-1.0 / r - w11) + w12
        g2 = 2.0 * r**2.0 * (-1.0 / r - w21) + w22
        y1 = ((y1_sub2 * (g1_sub2 - (12.0 / step_sq))
               + y1_sub1 * (10.0 * g1_sub1 + (24.0 / step_sq)))
              / ((12.0 / step_sq) - g1))
        y2 = ((y2_sub2 * (g2_sub2 - (12.0 / step_sq))
               + y2_sub1 * (10.0 * g2_sub1 + (24.0 / step_sq)))
              / ((12.0 / step_sq) - g2))

        # check for divergence
        if r < r_in1:
            dy1 = abs((y1 - y1_sub1) / y1_sub1)
            if dy1 > dr1:
                break

        if r < r_in2:
            dy2 = abs((y2 - y2_sub1) / y2_sub1)
            if dy2 > dr2:
                break

        # store vals
        norm1 += y1**2.0 * r**2.0
        norm2 += y2**2.0 * r**2.0
        integral += y1 * y2 * r**(2.0 + p)

        # next iteration
        r_sub1 = r

        g1_sub2 = g1_sub1
        g1_sub1 = g1
        g2_sub2 = g2_sub1
        g2_sub1 = g2

        y1_sub2 = y1_sub1
        y1_sub1 = y1
        y2_sub2 = y2_sub1
        y2_sub1 = y2
        i += 1
    return integral * (norm1 * norm2)**-0.5
