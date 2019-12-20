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
import numpy as np
from scipy.optimize import leastsq
from astropy.io import ascii as ac2
import MathLib as ml

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams as rc
from matplotlib.pyplot import style

style.use('seaborn-bright')
rc['text.usetex'] = True
rc['xtick.direction'] = 'in'
rc['ytick.direction'] = 'in'
rc['xtick.major.size']=4*3
rc['xtick.minor.size']=2*3
rc['ytick.major.size']=3*3
rc['ytick.minor.size']=1.5*3
rc['xtick.minor.visible']=True
rc['ytick.minor.visible']=True
rc['xtick.top'] = True
rc['ytick.right'] = True
rc['savefig.bbox'] = 'tight'
rc['font.size'] = 10
rc['lines.linewidth'] = 1.5
rc['axes.labelsize'] = 12



def fgauss(p, x):
    n_peak = int(len(p)/3)
    # A[i] = p[3*i]
    # xc[i] = p[3*i+1]
    # w[i] = p[3*i+2]
    result = 0
    for i in range(n_peak):
        result += p[3*i]*np.exp(-(x-p[3*i+1])*(x-p[3*i+1])*8*np.log(2)/p[3*i+2]/p[3*i+2]/2)
    return result

def gsres(p, y, x):
    return y - fgauss(p, x)

def gsfit(p0, x, y):
    par = leastsq(gsres, p0, args=(y, x), full_output=1)
    n_par = len(p0)
    dof = len(x) - n_par

    fvec = par[2]["fvec"]
    red_RSS = np.sum(fvec*fvec)/dof

    err = np.zeros(n_par)
    for i in range(n_par):
        err[i] = (red_RSS*par[1][i][i])**0.5

    return par[0], err, red_RSS


def prepare():

    # re-generate LIST.fit.pfm with a larger tolerance

    with open('config', 'r') as f:
        lines = f.read()[:-1].split('\n')
    idx = lines.index('# BEGIN FLAG OPTIONS') + 1
    bw = float(lines[idx].split()[-1])
    doffset = float(lines[idx+1].split()[-1])*2
    dwidth = float(lines[idx+2].split()[-1])*2
    rela_err = float(lines[idx+3].split()[-1])*2
    dsym = float(lines[idx+4].split()[-1])*2
    davg = float(lines[idx+5].split()[-1])*2
    del lines

    DataFile = os.listdir('fit')
    DataFile = [e for e in DataFile if e.endswith('.fit')]
    DataFile.sort()

    outf = open('performance/LIST.fit.pfm', 'w')

    for file in DataFile:
        Amp, Eamp = np.array([]), np.array([])
        Off, Eoff = np.array([]), np.array([])
        Hpbw, Ehpbw = np.array([]), np.array([])
        f = open('fit/'+file, "r")
        trs = f.readlines()
        trs = [e[:-1] for e in trs if e[0] != "#"]
        f.close()
        for i in range(len(trs)):
            tr = trs[i].split()
            Amp = np.append(Amp, float(tr[7]))
            Eamp = np.append(Eamp, float(tr[8]))
            Off = np.append(Off, float(tr[9]))
            Eoff = np.append(Eoff, float(tr[10]))
            Hpbw = np.append(Hpbw, float(tr[11]))
            Epbw = np.append(Ehpbw, float(tr[12]))
            try: reject = (Amp[i]>0)*(Eamp[i]/Amp[i]<rela_err)*\
                          (np.fabs(Off[i])<doffset)*\
                          (np.fabs(Hpbw[i]/bw-1)<dwidth)
            except: reject = False
            if reject == False:
                trs[i] = ''

        trs = [e for e in trs if e != '']
        trs = ("\n").join(trs)
        if trs:
            print(trs, file=outf)
    outf.close()

    #os.popen('grep -v "^#" LIST.fit > LIST.fit.pfm').read()
    scannum = os.popen("awk '{print $2}' performance/LIST.fit.pfm").read().split('\n')
    scannum = scannum[:-1]
    scan = list(set(scannum))
    scan.sort(key = scannum.index)
    scannum = np.array(scannum)

    bad = 1

    for s in scan:
        fscan = scannum == s
        rscan = scannum[fscan]
        n = len(rscan)

        scandir = os.popen("sed -n '%d,%dp' performance/LIST.fit.pfm | awk '{print $3}'"
                       %(bad, bad+n-1)).read().split('\n')[0:n]
        nalon = len([e for e in scandir if e == 'ALON'])
        nalat = n - nalon
        dn = nalon - nalat

        if dn > 0:
            os.popen('sed -i "%d,%dd" performance/LIST.fit.pfm' %(bad, bad+dn-1)).read()
        elif dn < 0:
            os.popen('sed -i "%d,%dd" performance/LIST.fit.pfm' %(bad+n+dn, bad+n-1)).read()

        bad += n - abs(dn)


def rms():
    with open('LIST.fit') as flf:
        lines = flf.readlines()
    lines = [e[:-1] for e in lines[3:] if not e[:-1][-1].endswith('#')]
    mark = np.array([e[0] for e in lines])
    qlt = np.array([(e==' ') and True or False for e in mark])
    lines = [e[1:].split() for e in lines]
    srcs = np.array([e[0] for e in lines])
    scans = np.array([e[1] for e in lines])
    direct = np.array([e[2] for e in lines])
    mjd = np.array([float(e[4]) for e in lines])
    azi = np.array([float(e[5]) for e in lines])
    elv = np.array([float(e[6]) for e in lines])
    amp = np.array([float(e[7]) for e in lines])
    err = np.array([float(e[8]) for e in lines])
    off = np.array([float(e[9]) for e in lines])
    tsys = np.array([float(e[14]) for e in lines])

    snr = amp/err

    with open('performance/SNR.ALL', 'w') as psnr:
        print('# %7s %8s %7s %7s %12s %9s %11s %12s %6s %7s' \
                %('Source', 'Scan', 'direct', 'MJD', 'Azi',
                'Elv', 'Amp', 'Offset', 'Tsys', 'SNR'), file=psnr)
        for i in range(len(scans)):
            print('%s %-11s %4s %6s %15.8f %9.3f %8.3f %12.5e ' \
                    '%7.2f %8.3f %6.2f' \
            %(mark[i], srcs[i], scans[i], direct[i], mjd[i], azi[i], elv[i],
            amp[i], off[i], tsys[i], snr[i]), file=psnr)

    # get and plot snr for primary calibrators
    pcal = ['3C286', '3C48', 'NGC7027']
    ff = open('performance/SNR.CAL', 'w')
    print('# %6s %10s %7s %12s %9s %11s %10s ' \
            '%13s %5s %8s %8s %6s %6s' \
            %('Source', 'direct', 'MJD', 'Azi','Elv', 'Amp', 'err',
            'Offset', 'err', 'Tsys', 'err', 'SNR', 'err'), file=ff)
    for cal in pcal:
        flt1 = (srcs==cal)*(direct=='ALON')*qlt
        flt2 = (srcs==cal)*(direct=='ALAT')*qlt
        if (True in flt1) and (True in flt2):
            print(' %-9s %6s %15.8f %9.3f %8.3f %12.5e %12.5e '\
                    '%7.2f %7.2f %8.3f %8.3f %6.2f %6.2f' \
                    %(cal, 'ALON', np.mean(mjd[flt1]), np.mean(azi[flt1]),
                            np.mean(elv[flt1]), np.mean(amp[flt1]),
                            np.std(amp[flt1]), np.mean(off[flt1]),
                            np.std(off[flt1]), np.mean(tsys[flt1]),
                            np.std(tsys[flt1]), np.mean(snr[flt1]),
                            np.std(snr[flt1])), file=ff)
            print(' %-9s %6s %15.8f %9.3f %8.3f %12.5e %12.5e '\
                    '%7.2f %7.2f %8.3f %8.3f %6.2f %6.2f' \
                    %(cal, 'ALAT', np.mean(mjd[flt2]), np.mean(azi[flt2]),
                            np.mean(elv[flt2]), np.mean(amp[flt2]),
                            np.std(amp[flt2]), np.mean(off[flt2]),
                            np.std(off[flt2]), np.mean(tsys[flt2]),
                            np.std(tsys[flt2]), np.mean(snr[flt2]),
                            np.std(snr[flt2])), file=ff)
    ff.close()


def offset(show=0):
    d = ac2.read('performance/LIST.fit.pfm', names = [str(e) for e in range(17)])
    filter1 = d['2'] == 'ALON'
    filter2 = d['2'] == 'ALAT'
    dd1 = d['9'][filter1]
    dd2 = d['9'][filter2]
    del filter1, filter2, d
    n1, n2 = len(dd1), len(dd2)

    x1lim0, x1lim1 = -3*np.std(dd1)+np.mean(dd1), 3*np.std(dd1)+np.mean(dd1)
    x2lim0, x2lim1 = -3*np.std(dd2)+np.mean(dd2), 3*np.std(dd2)+np.mean(dd2)
    binsize = int(max(x1lim1-x1lim0, x2lim1-x2lim0)/19)

    fig = plt.figure(figsize=(9.0, 3.6))
    gs = gridspec.GridSpec(1,2,wspace=0.1)
    ax0 = plt.subplot(gs[0])

    ax0.set_xlim(x1lim0, x1lim1)
    counts, bins, patches = plt.hist(dd1,
                    bins = np.arange(x1lim0, x1lim1+binsize, binsize),
                    weights=np.ones(n1)/n1*100,
                    density=0, histtype='step', linewidth=2)
    x = (bins[1:]+bins[:-1])/2.0

    p0 = [max(counts), 0, (8*np.log(2))**0.5*np.std(dd1)]
    try:
        par, err, rss = gsfit(p0, x, counts)
    except:
        par, err, rss = None, None, None

    ax0.set_ylim(0, 1.1*max(counts))
    ax0.set_xlabel(r'\rm ALON offset (arcsec)')
    ax0.set_ylabel(r'\rm fraction (\%)')

    plotx = np.linspace(x1lim0, x1lim1, 300)
    # rectangle for +- 3sigma
    if par is not None:
        sigma0 = par[2]/(8*np.log(2))**0.5

        plt.plot(plotx, fgauss(par, plotx), color='r', ls='--')

    # plot parameters
        ax0.text(0.05, 0.88, r'$\rm mean=%0.1f$' %par[1], transform=ax0.transAxes, fontsize=15)
        ax0.text(0.05, 0.78, r'$\sigma=%0.1f$' %sigma0, transform=ax0.transAxes, fontsize=15)

    ax1 = plt.subplot(gs[1])
    ax1.set_xlim(x2lim0, x2lim1)
    counts, bins, patches = plt.hist(dd2,
                    bins = np.arange(x2lim0, x2lim1+binsize, binsize),
                    weights=np.ones(n2)/n2*100,
                    density=0, histtype='step', linewidth=2)
    x = (bins[1:]+bins[:-1])/2.0

    p0 = [max(counts), 0, (8*np.log(2))**0.5*np.std(dd2)]
    try:
        par, err, rss = gsfit(p0, x, counts)
    except:
        par, err, rss = None, None, None

    ax1.set_ylim(0, 1.1*max(counts))
    ax1.set_xlabel(r'\rm ALAT offset (arcsec)')

    plotx = np.linspace(x2lim0, x2lim1, 300)
    # rectangle for +- 3sigma
    if par is not None:
        sigma0 = par[2]/(8*np.log(2))**0.5
        plt.plot(plotx, fgauss(par, plotx), color='r', ls='--')

    # plot parameters
        ax1.text(0.05, 0.88, r'$\rm mean=%0.1f$' %par[1], transform=ax1.transAxes, fontsize=15)
        ax1.text(0.05, 0.78, r'$\sigma=%0.1f$' %sigma0, transform=ax1.transAxes, fontsize=15)

    plt.savefig('performance/Offset.png', bbox_inches='tight', dpi=300)
    if show:
        plt.show()




def pointing(show=0):

    figname = 'performance/Pointing.png'
    d = ac2.read('performance/LIST.fit.pfm', names = [str(e) for e in range(17)])
    filter1 = d['2'] == 'ALON'
    filter2 = d['2'] == 'ALAT'
    az = d['5'][filter1]
    el = d['6'][filter1]
    off1 = d['9'][filter1]
    off2 = d['9'][filter2]
    del d
    az_mean, az_rms = np.mean(off1), np.std(off1)
    el_mean, el_rms = np.mean(off2), np.std(off2)
    az_ylim0, az_ylim1 = az_mean - 3*az_rms, az_mean + 3*az_rms
    el_ylim0, el_ylim1 = el_mean - 3*el_rms, el_mean + 3*el_rms

    fig = plt.figure(figsize=(7.2, 5.4))
    gs = gridspec.GridSpec(2,2,hspace=0.0,wspace=0.0)

    ax0 = plt.subplot(gs[0])
    ax0.set_xlim(0, 360)
    ax0.set_ylim(az_ylim0, az_ylim1)
    ax0.xaxis.set_ticks(np.linspace(0, 360, 7))
#    ax0.yaxis.set_ticks(np.arange(az_ylim0, az_ylim1, az_rms))
    ax0.set_xticklabels([])
    ax0.set_ylabel(r'$\rm \Delta\ Azi\ (arcsec.)$')
    ax0.scatter(az, off1, marker='o', s=25, edgecolors='none', alpha=0.5)
#    ax0.plot([0, 360], [0, 0], color='r', ls='--')

    ax1 = plt.subplot(gs[1])
    ax1.set_xlim(0, 90)
    ax1.set_ylim(az_ylim0, az_ylim1)
    ax1.xaxis.set_ticks(np.linspace(0, 90, 7))
#    ax1.yaxis.set_ticks(np.arange(az_ylim0, az_ylim1, az_rms))
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.scatter(el, off1, marker='o', s=25, edgecolors='none', alpha=0.5)
#    ax1.plot([0, 90], [0, 0], color='r', ls='--')

    ax2 = plt.subplot(gs[2])
    ax2.set_xlim(0, 360)
    ax2.set_ylim(el_ylim0, el_ylim1)
    ax2.xaxis.set_ticks(np.linspace(0, 360, 7))
#    ax2.yaxis.set_ticks(np.arange(el_ylim0, el_ylim1, 15))
    ax2.set_xlabel(r'$\rm Azimuth\ (deg.)$')
    ax2.set_ylabel(r'$\rm \Delta\ Elv\ (arcsec.)$')
    ax2.scatter(az, off2, marker='o', s=25, edgecolors='none', alpha=0.5)
#    ax2.plot([0, 360], [0, 0], color='r', ls='--')

    ax3 = plt.subplot(gs[3])
    ax3.set_xlim(0, 90)
    ax3.set_ylim(el_ylim0, el_ylim1)
    ax3.set_yticklabels([])
    ax3.xaxis.set_ticks(np.linspace(0, 90, 7))
#    ax3.yaxis.set_ticks(np.arange(el_ylim0, el_ylim1, 15))
    ax3.set_xlabel(r'$\rm Elevation\ (deg.)$')
    ax3.scatter(el, off2, marker='o', s=25, edgecolors='none', alpha=0.5)
#    ax3.plot([0, 90], [0, 0], color='r', ls='--')

    plt.savefig(figname)
    if show:
        plt.show()



def off_vs_hpbw(show=0):
    d = ac2.read('performance/LIST.fit.pfm', names = [str(e) for e in range(17)])
    d1 = d['9']
    d2 = d['11']
    n = len(d2)
    xlim0, xlim1 = -3*np.std(d1), 3*np.std(d1)
    ylim0, ylim1 = -3*np.std(d2) + np.mean(d2), 3*np.std(d2) + np.mean(d2)

    fig = plt.figure(figsize=(6.8, 4.2))
    gs = gridspec.GridSpec(2, 2, height_ratios=[2,1.5], hspace=0.0,
                                width_ratios=[2,1], wspace=0.0)
    ax = list(range(4))
    ax[0] = plt.subplot(gs[0])
    ax[0].set_xlim(xlim0, xlim1)
    ax[0].set_ylim(ylim0, ylim1)
    ax[0].set_xticklabels([])
#    ax[0].yaxis.set_ticks(np.linspace(ylim0, ylim1, 10))
    ax[0].set_ylabel(r'$\rm HPBW\ (arcsec)$')
    ax[0].scatter(d1, d2, marker='s', s=15, edgecolors='none', alpha=0.2)

    ax[1] = plt.subplot(gs[1])
    ax[1].set_ylim(ylim0, ylim1)
    ax[1].set_yticklabels([])

    counts, bins, patches = plt.hist(d2,
                    bins = np.arange(ylim0, ylim1, (ylim1-ylim0)/20),
                    weights=np.ones(n)/n*100,
                    density=0, histtype='step', linewidth=2,
                    orientation='horizontal')


    ax[2] = plt.subplot(gs[2])
    ax[2].set_xlim(xlim0, xlim1)
    ax[2].set_xlabel(r'$\rm offset\ (arcsec)$')
    counts, bins, patches = plt.hist(d1,
                    bins = np.arange(xlim0, xlim1+4, (xlim1-xlim0)/20),
                    weights=np.ones(n)/n*100,
                    density=0, histtype='step', linewidth=2)

    plt.savefig('performance/Off_HPBW.png', bbox_inches='tight', dpi=300)
    if show:
        plt.show()


def hpbw(show=0):

    d = ac2.read('performance/LIST.fit.pfm', names = [str(e) for e in range(17)])
    filter1 = d['2'] == 'ALON'
    filter2 = d['2'] == 'ALAT'
    dd1 = d['11'][filter1]
    dd2 = d['11'][filter2]
    del filter1, filter2, d
    n1, n2 = len(dd1), len(dd2)

    x1lim0, x1lim1 = -3*np.std(dd1)+np.mean(dd1), 3*np.std(dd1)+np.mean(dd1)
    x2lim0, x2lim1 = -3*np.std(dd2)+np.mean(dd2), 3*np.std(dd2)+np.mean(dd2)
    binsize = int(max(x1lim1-x1lim0, x2lim1-x2lim0)/24)
    if binsize == 0: binsize = 1

    fig = plt.figure(figsize=(9.0, 3.6))
    gs = gridspec.GridSpec(1,2,wspace=0.1)
    ax0 = plt.subplot(gs[0])

    ax0.set_xlim(x1lim0, x1lim1)
    counts, bins, patches = plt.hist(dd1,
                    bins = np.arange(x1lim0, x1lim1+binsize, binsize),
                    weights=np.ones(n1)/n1*100,
                    density=0, histtype='step', linewidth=2)
    x = (bins[1:]+bins[:-1])/2.0

    p0 = [max(counts), np.mean(dd1), (8*np.log(2))**0.5*np.std(dd1)]
    try:
        par, err, rss = gsfit(p0, x, counts)
    except:
        par, err, rss = None, None, None

    ax0.set_ylim(0, 1.1*max(counts))
    ax0.set_xlabel(r'\rm ALON HPBW (arcsec)')
    ax0.set_ylabel(r'\rm fraction (\%)')

    plotx = np.linspace(x1lim0, x1lim1, 300)
    # rectangle for +- 3sigma
    if par is not None:
        sigma0 = par[2]/(8*np.log(2))**0.5

        plt.plot(plotx, fgauss(par, plotx), color='r', ls='--')

    # plot parameters
        ax0.text(0.05, 0.88, r'$\rm mean=%0.1f$' %par[1], transform=ax0.transAxes, fontsize=15)
        ax0.text(0.05, 0.78, r'$\sigma=%0.1f$' %sigma0, transform=ax0.transAxes, fontsize=15)

    ax1 = plt.subplot(gs[1])
    ax1.set_xlim(x2lim0, x2lim1)
    counts, bins, patches = plt.hist(dd2,
                    bins = np.arange(x2lim0, x2lim1+binsize, binsize),
                    weights=np.ones(n2)/n2*100,
                    density=0, histtype='step', linewidth=2)
    x = (bins[1:]+bins[:-1])/2.0

    p0 = [max(counts), np.mean(dd2), (8*np.log(2))**0.5*np.std(dd2)]
    try:
        par, err, rss = gsfit(p0, x, counts)
    except:
        par, err, rss = None, None, None

    ax1.set_ylim(0, 1.1*max(counts))
    ax1.set_xlabel(r'\rm ALAT HPBW (arcsec)')

    plotx = np.linspace(x2lim0, x2lim1, 300)
    # rectangle for +- 3sigma
    if par is not None:
        sigma0 = par[2]/(8*np.log(2))**0.5
        plt.plot(plotx, fgauss(par, plotx), color='r', ls='--')

    # plot parameters
        ax1.text(0.05, 0.88, r'$\rm mean=%0.1f$' %par[1], transform=ax1.transAxes, fontsize=15)
        ax1.text(0.05, 0.78, r'$\sigma=%0.1f$' %sigma0, transform=ax1.transAxes, fontsize=15)

    plt.savefig('performance/HPBW.png', bbox_inches='tight', dpi=300)
    if show:
        plt.show()


def tsys(tcal, show):

    with open('config', 'r') as f:
        lines = f.readlines()

    cal = [e for e in lines if 'Tcal' in e][0]
    cal = float(cal.split()[1])
    ratio = tcal/cal

    d = ac2.read('performance/LIST.fit.pfm', names = [str(e) for e in range(17)])
    d1 = d['6']
    d2 = d['14']*ratio
    del d
    func = ml.lower_envelope(1/np.sin(np.deg2rad(d1)), d2)

    fig = plt.figure(figsize=(6.8, 3.2))
    gs = gridspec.GridSpec(1,2,wspace=0.0)

    ax0 = plt.subplot(gs[0])
    ax0.xaxis.set_ticks(np.arange(10, 80, 15))
    ax0.set_xlabel(r'$\rm Elevation\ (deg.)$')
    ax0.set_ylabel(r'$\rm T_{sys}\ (K)$')
    ax0.scatter(d1, d2, marker='o', s=15, edgecolors='none', alpha=0.5)

    ax1 = plt.subplot(gs[1])
    ax1.set_xlim(0.8, 4.0)
    ax1.set_ylim(ax0.get_ylim())
    ax1.set_yticklabels([])
    ax1.set_xlabel(r'$\rm Airmass$')
    ax1.scatter(1/np.sin(np.deg2rad(d1)), d2, marker='o', s=15, edgecolors='none', alpha=0.5)
    plotx = np.linspace(0.8, 4.0, 10)
    ax1.plot(plotx, np.polyval(func, plotx), color='r', ls='--')
    ax1.text(0.05, 0.88, r'$\rm T_{sys}=%0.2f\,K+%0.3f*Airmass$' %(func[1], func[0]), transform=ax1.transAxes, fontsize=12)
    ax1.text(0.05, 0.78, r'$\rm T_{zenith}=%0.1f\,K$' %np.polyval(func, 1), transform=ax1.transAxes, fontsize=12)

    plt.savefig('performance/Tsys.png', bbox_inches='tight', dpi=300)

    if show:
        plt.show()

    return func, ratio


def hpbw_value():

    d = ac2.read('performance/LIST.fit.pfm', names = [str(e) for e in range(17)])
    dd1 = d['11']
    n1 = len(dd1)
    del d

    x1lim0, x1lim1 = -3*np.std(dd1)+np.mean(dd1), 3*np.std(dd1)+np.mean(dd1)
    binsize = int((x1lim1-x1lim0)/24)
    if binsize == 0: binsize = 1

    fig = plt.figure(figsize=(9.0, 3.6))
    gs = gridspec.GridSpec(1,2,wspace=0.1)
    ax0 = plt.subplot(gs[0])

    ax0.set_xlim(x1lim0, x1lim1)
    counts, bins, patches = plt.hist(dd1,
                    bins = np.arange(x1lim0, x1lim1+binsize, binsize),
                    weights=np.ones(n1)/n1*100,
                    density=0, histtype='step', linewidth=2)
    x = (bins[1:]+bins[:-1])/2.0

    p0 = [max(counts), np.mean(dd1), (8*np.log(2))**0.5*np.std(dd1)]
    try:
        par, err, rss = gsfit(p0, x, counts)
    except:
        par, err, rss = [None, None, None], None, None

    return par[1], par[2]


def off_inspect(show=0):
    d = ac2.read('performance/LIST.fit.pfm', names = [str(e) for e in range(17)])
    scans = d['1']
    scandir = d['2']
    d1 = d['4']
    d2 = d['9']
    d0 = d['0']
    d7 = d['7']
    scannum = list(set(scans))
    scans = np.array(scans)

    mjd_alon = [np.mean(d1[(scans==s)*(scandir=='ALON')]) for s in scannum]
    mjd_alat = [np.mean(d1[(scans==s)*(scandir=='ALAT')]) for s in scannum]
    off_alon = [np.mean(d2[(scans==s)*(scandir=='ALON')]) for s in scannum]
    off_alat = [np.mean(d2[(scans==s)*(scandir=='ALAT')]) for s in scannum]
    mjd0 = int(min(min(mjd_alon), min(mjd_alat)))

    fig = plt.figure(figsize=(9.0, 3.6))
    gs = gridspec.GridSpec(1,1,wspace=0.0)
    ax = plt.subplot(gs[0])

    ax.plot(np.array(mjd_alon)-mjd0, off_alon, 'ro', label='ALON')
    ax.plot(np.array(mjd_alat)-mjd0, off_alat, 'bs', label='ALAT')
    ax.legend(loc='best')
    ax.set_xlabel(r'MJD-%d' %mjd0)
    ax.set_ylabel(r'arcsec.')
    plt.savefig('performance/Off_MJD.png', bbox_inches='tight', dpi=300)

    cals = ['3C286', '3C48', 'NGC7027']
    fig = plt.figure(figsize=(4.5, 7))
    gs = gridspec.GridSpec(3,1,hspace=0.15)
    for i in range(len(cals)):
        x_lon = [np.mean(d2[(d0==cals[i])*(scandir=='ALON')*(scans==s)]) for s in scannum]
        x_lat = [np.mean(d2[(d0==cals[i])*(scandir=='ALAT')*(scans==s)]) for s in scannum]
        y_lon = [np.mean(d7[(d0==cals[i])*(scandir=='ALON')*(scans==s)]) for s in scannum]
        y_lat = [np.mean(d7[(d0==cals[i])*(scandir=='ALAT')*(scans==s)]) for s in scannum]
        ax = plt.subplot(gs[i])
#    ax.errorbar(x, y, yerr=yerr, fmt='o', capthick=1, label=cals[i])
        ax.plot(x_lat, y_lon, 'ro', label=cals[i]+' ALON')
        ax.plot(x_lon, y_lat, 'bs', label=cals[i]+' ALAT')
        ax.set_ylabel(r'Amplitude [counts]')
        if i == 2:
            ax.set_xlabel(r'MJD')
        ax.legend(loc='best')

    plt.savefig('performance/Off_Amp.png', bbox_inches='tight', dpi=300)



def sys_pfm(D, F, Tcal, Tgr):

    if not os.path.exists('performance'):
        os.mkdir('performance')

    pfile = open('performance/performance', 'w')

    if not os.path.isfile('Performance/LIST.fit.pfm'):
        prepare()


    pointing()
    off_vs_hpbw()
    offset()
    hpbw()
    off_inspect()
    tsys_par, ratio = tsys(Tcal, 0)

    Kbol = 1.3806488E-23
    Tatm = (Tgr+273.15)*1.12-50

    print('Inpur Parameters')
    print('-------------------')
    print('%-6s = %8.1f [meter]' %('D', D))
    print('%-6s = %8.1f [MHz]' %('F', F))
    print('%-6s = %8.3f [K]' %('Tcal', Tcal))
    print('%-6s = %8.1f [deg. cent.]\n' %('Tgr', Tgr))

    print('%-6s = %8.1f [k]' %('Tatm', Tatm))
    print('')

    print('System Performance')
    print('-------------------')

    print('Inpur Parameters', file=pfile)
    print('-------------------', file=pfile)
    print('%-6s = %8.1f [meter]' %('D', D), file=pfile)
    print('%-6s = %8.1f [MHz]' %('F', F), file=pfile)
    print('%-6s = %8.3f [K]' %('Tcal', Tcal), file=pfile)
    print('%-6s = %8.1f [deg. cent.]\n' %('Tgr', Tgr), file=pfile)

    print('%-6s = %8.1f [k]' %('Tatm', Tatm), file=pfile)
    print('', file=pfile)

    print('System Performance', file=pfile)
    print('-------------------', file=pfile)


    Tsys0 = np.polyval(tsys_par, 1)
    print('%-6s = %8.1f [k]' %('Tsys0', Tsys0))
    print('%-6s = %8.4f' %('Tau0', tsys_par[0]/Tatm))
    print('%-6s = %8.1f [k]' %('Tsys0', Tsys0), file=pfile)
    print('%-6s = %8.4f' %('Tau0', tsys_par[0]/Tatm), file=pfile)

    with open('cal.factor') as inputf:
        line = inputf.readlines()[-1]
    DPFU = float(line.split()[4])*ratio
    SEFD = Tsys0/DPFU
    print('%-6s = %8.3f [K/Jy]' %('DPFU', DPFU))
    print('%-6s = %8.3f [Jy]\n' %('SEFD', SEFD))
    print('%-6s = %8.3f [K/Jy]' %('DPFU', DPFU), file=pfile)
    print('%-6s = %8.3f [Jy]\n' %('SEFD', SEFD), file=pfile)

    Ae = 2*Kbol*DPFU*1E26
    EffA = Ae/np.pi/(D/2.0)**2
    print('%-6s = %8.3f [meter^2]' %('Ae', Ae))
    print('%-6s = %8.3f [%%]' %('Eff.A', EffA*100))
    print('%-6s = %8.3f [meter^2]' %('Ae', Ae), file=pfile)
    print('%-6s = %8.3f [%%]' %('Eff.A', EffA*100), file=pfile)

    HPBW = hpbw_value()
    Omega_mb = np.pi/4/np.log(2)*np.deg2rad(HPBW[0]/3600.0)**2
    Omega_A = (300/F)**2/Ae
    EffB = Omega_mb/Omega_A
    print('%-6s = %8.1f +-%5.1f [arcsec]' %('HPBW', HPBW[0], HPBW[1]))
    print('%-6s = %8.3E [rad^2]' %('Omega_A', Omega_A))
    print('%-6s = %8.3E [rad^2]' %('Omega_mb', Omega_mb))
    print('%-6s = %8.3f [%%]' %('Eff.mb', EffB*100))
    print('')
    print('%-6s = %8.1f +-%5.1f [arcsec]' %('HPBW', HPBW[0], HPBW[1]), file=pfile)
    print('%-6s = %8.3E [rad^2]' %('Omega_A', Omega_A), file=pfile)
    print('%-6s = %8.3E [rad^2]' %('Omega_mb', Omega_mb), file=pfile)
    print('%-6s = %8.3f [%%]' %('Eff.mb', EffB*100), file=pfile)
    print('', file=pfile)

    Tmb_per_S = DPFU/EffB
    print('%-6s = %8.3f [K/Jy]' %('Tmb/s', Tmb_per_S))
    print('%-6s = %8.3f [K/Jy]' %('Tmb/s', Tmb_per_S), file=pfile)

    gain_file = os.listdir('./')
    gain_file = [e for e in gain_file if e.endswith('gain.curve')][0]
    with open(gain_file) as gf:
        gain = gf.readlines()[3].split()[2:]
    gain = [float(e) for e in gain]
    print('Gain = A*x**2 + B*x + C')
    print('%8s %12s %12s' %('A', 'B', 'C'))
    print('%12.5E %12.5E %12.5E' %(gain[0], gain[1], gain[2]))
    print('peak@ %0.1f deg.' %(-gain[1]/2/gain[0]))
    print('Gain = A*x**2 + B*x + C', file=pfile)
    print('%8s %12s %12s' %('A', 'B', 'C'), file=pfile)
    print('%12.5E %12.5E %12.5E' %(gain[0], gain[1], gain[2]), file=pfile)
    print('peak@ %0.1f deg.' %(-gain[1]/2/gain[0]), file=pfile)



if __name__ == '__main__':
    sys_pfm(26.0, 23500.0, 5.0, 5.0)
