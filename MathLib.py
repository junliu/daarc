#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2015-2020 Jun Liu <jliu@mpifr-bonn.mpg.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function

__author__ = 'Jun LIU'
__copyright__ = 'Copyright (c) 2015-2020 Jun Liu <liuj@xao.ac.cn, jliu@mpifr-bonn.mpg.de>'
__license__ = 'GPL v3'
__version__ = '2.1'


import os
import subprocess
import inspect
import numpy as np
from astropy.io import fits as pyfits
from scipy.optimize import leastsq, curve_fit, minimize
from scipy import interpolate
import scipy.stats as stats
import ephem
import copy
import re



class Gauss:

    def __init__(self, npeak, bo):
        self.npeak = npeak
        self.bo = bo

    def fGauss(self, p, x):

        val = 0
        for i in range(self.npeak):
            val += p[i*3]*np.exp(-(x-p[i*3+1])**2/p[i*3+2]**2*4*np.log(2))
        val += np.polyval(p[3*self.npeak:3*self.npeak+self.bo+1], x)

        return val


    def fResidual(self, p, x, y):
        return y - self.fGauss(p, x)

    def fResidual_wt(self, p, x, y, wt):

        factor = len(wt)/np.sum(wt)
        return (y - self.fGauss(p, x))*wt*factor


    def fJaco(self, p, x, y):

        jac = []
        for i in range(self.npeak):
            a, x0, w = p[i*3], p[i*3+1], p[i*3+2]
            F = a*np.exp(-4.0*np.log(2.0)*(x-x0)**2.0/w**2.0)
            jac.append(F/a)
            jac.append(F*8*np.log(2.0)*(x-x0)/w**2.0)
            jac.append(F*8*np.log(2.0)*(x-x0)**2.0/w**3.0)

        for i in range(self.bo, -1, -1):
            jac.append(x**i)

        return -np.array(jac)


    def fit_basic(self, x, y, p0, despike):

        dof = len(x) - self.npeak*3 - self.bo -2
        scale = min(max(x)-min(x), max(y)-min(y))
        if scale >= 1: scale = 1
        tol = scale*1E-9
        if not despike:
            pars, cov_x, infodict, errmsg, ier = leastsq(self.fResidual, \
                p0, args=(x, y), xtol=tol, ftol=tol, \
                #col_deriv=1, full_output=1)
                Dfun=self.fJaco, col_deriv=1, full_output=1)
        else:
            idx = x.argsort()
            x = x[idx]
            y = y[idx]
            k = (y[1:]-y[:-1])/(x[1:]-x[:-1])
            lin = np.polyfit(x[1:], k, 1)
            dk = np.fabs(k - np.polyval(lin, x[1:]))
            wt = dk**-0.1
            pars, cov_x, infodict, errmsg, ier = leastsq(self.fResidual_wt, \
                    p0, args=(x[1:], y[1:], wt), xtol=tol, ftol=tol, \
                    col_deriv=1, full_output=1)
            dfit = self.fGauss(pars, x[1:]) - y[1:]
            flt = dfit > 0
            wt = dk**-1
            dfit_mean = np.average(dfit[flt], weights=wt[flt])
            dfit_rms = np.std(dfit[flt])

            dist = np.fabs(dfit - dfit_mean)
            dist[dist<despike*dfit_rms] = dfit_rms
            wt = dk**-0.1*dist**-1.0
            pars, cov_x, infodict, errmsg, ier = leastsq(self.fResidual_wt, \
                    pars, args=(x[1:], y[1:], wt), xtol=tol, ftol=tol, \
                    col_deriv=1, full_output=1)

        for i in range(self.npeak):
            pars[3*i+2] = np.fabs(pars[3*i+2])

        # estimate the convariance
        try:
            if ier in [1, 2, 3, 4]:
                fvec = infodict['fvec']
                red_rss = np.sum(fvec**2)/dof
                err = (np.diag(cov_x)*red_rss)**0.5
            else:
                print('poor fitting!!!')
                err = -np.ones(3*self.npeak+self.bo+1)

            return pars, err

        except ValueError:

            p0 = -np.ones(len(p0))
            p0[-1] = min(y)
            return p0, -np.ones(len(p0))


    def fit(self, x, y, p0, despike):

        span = max(x) - min(x)
        if p0 is None:
            bp = fit_poly(x, y, 0, self.bo)[0]
            rp0 = np.zeros(4+self.bo)
            rp0[3+self.bo] = min(y) + bp[self.bo]

            for i in range(self.npeak):
                gs = Gauss((len(rp0)-2)/3, self.bo)
                rp0[i*3+2] = span/4.5
                yclean = y - gs.fGauss(rp0, x)
                rp0[i*3] = max(yclean) - min(yclean)
                rp0[i*3+1] = x[yclean.argmax()]
                rp0, rpe = gs.fit_basic(x, y, rp0, despike)
                if i == self.npeak - 1:
                    break
                rp0 = np.insert(rp0, i*3+3, np.zeros(3))

        elif len(p0)-self.bo-1<3*self.npeak:
            rp0 = np.array(p0)*1.0
            for i in range(self.npeak-(len(p0)-2)/3+1):
                n = (len(rp0)-2)/3
                gs = Gauss(n, self.bo)
                rp0, rpe = gs.FitBasic(x, y, rp0)
                if n == self.npeak:
                    break
                yclean = y - gs.FGauss(rp0, x)
                newp = np.zeros(3)
                newp[0] = max(yclean) - min(yclean)
                newp[1] = x[yclean.argmax()]
                newp[2] = span/4.5
                rp0 = np.insert(rp0, 3*n, newp)

        elif len(p0)-self.bo-1 == 3*self.npeak:
            rp0 = np.array(p0)*1.0
            rp0, rpe = self.fit_basic(x, y, rp0, despike)

        else:
            print('Number of peaks in SMODEL missmatch!!!')
            rp0 = -np.ones(self.npeak+self.bo+1)
            rp0[-1] = min(y)
            rpe = rp0*1.0

        return rp0, rpe


    def fit_window(self, x, y, p0, width, window, despike):

        rp0, rpe = self.fit(x, y, p0, despike)

        if rp0[0] == 0:
            return rp0, rpe

        # taking the 'best' peak
        wlist = np.array([(np.fabs(rp0[3*i+1])+width/20.0)\
                *np.fabs(rp0[3*i+2]-width) for i in range(self.npeak)])
        idx = wlist.argmin()

        x0, w = rp0[idx*3+1], rp0[idx*3+2]
        flt = x >= x0 - w*window[0]
        flt *= x <= x0 + w*window[1]

        if len([e for e in flt if e == True]) < 10:

            return rp0, rpe

        else:

            xfit = x[flt]
            yfit = y[flt]

            rp0_new = np.append(rp0[3*idx:3*idx+3], \
                    rp0[3*self.npeak:3*self.npeak+self.bo+1])
            gs = Gauss(1, self.bo)
            rp0_new, rpe_new = gs.fit(xfit, yfit, rp0_new, despike)

            if np.fabs(rp0_new[1]) >= window[0]*width or \
                np.fabs(rp0_new[2]/width-1)>0.6:
                return rp0, rpe

            else:
                return rp0_new, rpe_new


def fit_poly(x, y, yerr, order):

    N = len(x)
    NP = order + 1
    DOF = N - NP - 1

    if yerr == 0 or yerr == 1:
        yerr = np.ones(N)
    wt = N*yerr**-2.0/np.sum(yerr**-2.0)

    def f_residual(p, x, y):
        return wt*y - np.polyval(p, x)

    scale = min(max(x)-min(x), max(y)-min(y))
    if scale >= 1: scale = 1
    tol = scale*1E-9
    p0 = np.zeros(NP)
    pars, cov_x, infodict, errmsg, ier = leastsq(f_residual, p0, args=(x, y), \
    xtol=tol, ftol=tol, col_deriv=1, full_output=1)

    # estimate the convariance
    fvec = infodict['fvec']
    red_rss = np.sum(fvec**2)/DOF
    err = (np.diag(cov_x)*red_rss)**0.5

    return pars, err, red_rss**0.5

#     def lin2(x, a, b, c):
#         return a*x**2 + b*x + c

#     def lin3(x, a, b, c, d):
#         return a*x**3 + b*x**2 + c*x + d

#     if order == 3:
#         pars, cov_x, infodict, errmsg, ier = curve_fit(lin3, x, y,\
#         xtol=tol, ftol=tol, col_deriv=1, full_output=1)
#     else:
#         bounds = ((-np.inf, 0, -np.inf), (0, np.inf, np.inf))
# #        pars, cov_x, infodict, errmsg, ier = curve_fit(lin2, x, y, bounds=bounds,\
# #        xtol=tol, ftol=tol, col_deriv=1, full_output=1)
#         pars, pcov = curve_fit(lin2, x, y, bounds=bounds, xtol=tol, ftol=tol)

#         err = np.sqrt(np.diag(pcov))
#         red_rss = np.std(y - lin2(x, *pars))

#         return pars, err, red_rss




def RunningMedian(data, n):

    idx = np.arange(n) + np.arange(len(data)-n+1)[:, None]
    b = [e for e in data[idx]]
    return np.array(list(map(np.median, b)))


def DeSpike(data, window=5, sigma=2):

    N = len(data)
    median = RunningMedian(data, window)
    diff = np.fabs(data[window/2:N-window/2]-median)
    rms = np.std(diff)
    nn = len(diff)

    idx_outlier = diff >= rms*sigma
    idx_outlier = np.arange(nn)[idx_outlier]
    for i in idx_outlier:
        data[window/2+i] = median[i]

    return data


def lower_envelope(x, y):

#    def func(x, a, b):
#        return a*x + b
#
#    def get_flipped(ydata, ymodel):
#        flipped = ymodel - ydata
#        flipped[flipped>0] = 0
#        return flipped
#
#    def flipped_resid(pars, x, y):
#        ymodel = func(x, *pars)
#        flipped = get_flipped(y, ymodel)
#        resid = (y + flipped - ymodel)**0.5
#
#        return np.nan_to_num(resid)
#
#    guess = np.polyfit(x, y, 1)
#    pars, flag = leastsq(func=flipped_resid, x0=guess, args=(x, y))
#    flt = func(x, *pars) > y
#    return pars
#    dx = x.max() - x.min()
#    mx = np.mean(x)
#    x0, x1 = mx-0.005*dx, mx+0.005*dx
#    flt = (x>x0)*(x<x1)
#    _y = y[flt]
#    dy = _y.max() - _y.min()
#    print dy
#
#    def func(pars, x, y):
#        a0, b0, a1 = pars
#        d1 = np.fabs(y - (a0*x+b0))**0.4
#        d2 = np.fabs((a1*x+b0+dy) - y)**0.4
#        return d1+d2
#
#    guess = np.polyfit(x, y, 1)
#    a0, b0 = guess
#    a1 = a0
#    pars, flag = leastsq(func=func, x0=(a0, b0, a1), args=(x, y))
#    a = min(pars[0], pars[2])
#    b = pars[1]
#
#    return [a, b+dy-0.6]
#    return [a, b+0.5]

#    N = len(x)
#    filter = np.array([True]*N)
#    nx, ny = x[filter], y[filter]
#    while 1:
#        nn = len(nx[filter])
#        par = np.polyfit(nx[filter], ny[filter] ,1)
#        res = ny - np.polyval(par, nx)
#        filter = res <= np.std(res)*0
#        NN = len(nx[filter])
#        if NN < 0.05*N or NN < 5 or nn == NN: break
#    return nx[filter], ny[filter]

    dx = x.max() - x.min()
    flt0 = (x>0.05*dx)
    dx = x[flt0].max() - x[flt0].min()
    minx = x[flt0].min()
    nx, ny = [], []
    n = 8
    for i in range(n):
        flt = (x>minx+dx*i/n)*(x<=minx+dx*(i+1)/n)
        if True in flt:
            idx = np.argmin(y[flt])
            nx.append(x[flt][idx])
            ny.append(y[flt][idx])

    pars = np.polyfit(nx, ny, 1)
    return pars






