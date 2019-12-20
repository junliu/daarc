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
import re
import copy
import numpy as np
from scipy import interpolate, stats
from Misc import readcol, writecol, array_flat
import AstroLib as al
import MathLib as ml
import sample


class FitDat(object):

    def __init__(self, cfg, lnk, flag):

        self.cfg = cfg
        self.lfit = 'LIST.fit'
        self.lflux = 'LIST.flux'
        self.lnk = lnk
        self.flag = flag
        self.aver_tag = 0
        self.proc = './0.point'

    def load_fit(self, scannum):

        self.ffit = os.path.join(self.cfg.fit_path, '%04d.fit' %scannum)

        try:
            data = \
            readcol(self.ffit, fmt='sdsdffffffffffffffff', start=4, flag=True)
        except IndexError:
            print('No data contents in %04d.FITS' %scannum)
            return 0
        self.src = data[0][0]
        self.scannum = data[1][0]
        self.scandir = np.array(data[2])
        self.subnum = np.array(data[3])
        self.mjd = np.array(data[4])
        self.azi = np.array(data[5])
        self.elv = np.array(data[6])
        self.amp = np.array(data[7])
        self.aerr = np.array(data[8])
        self.off = np.array(data[9])
        self.oerr = np.array(data[10])
        self.hpbw = np.array(data[11])
        self.herr = np.array(data[12])
        self.pa = np.array(data[13])
        self.tsys = np.array(data[14])
        self.cal = np.array(data[15])
        self.sunang = np.array(data[16])
        del data
        return 1

    def load_dat(self, scannum):

        with open(os.path.join(self.cfg.fit_path, '%04d.dat' %scannum)) as fdat:
            trdat = fdat.readlines()

        trdat = [e[:-1] for e in trdat if (e[0]!='#' or e[:2]=='#!' or \
                (e[0]=='#' and len(e.split())==6))]
        trdat = trdat[1:]
        tr = [e for e in trdat if re.match('#! SUBSCAN', e)]
        keys = [trdat.index(e) for e in tr]
        keys.append(len(trdat))

        self.dat_offs = []
        self.dat_amps_raw = []
        self.dat_amps_subs = []
        self.dat_residual = []
        self.dat_u = []
        self.dat_q = []

        for i in range(len(tr)):
            trd = trdat[keys[i]+2:keys[i+1]]
            trd = [e.split()[::1] for e in trd]
            trd = np.array(trd).transpose()

            self.dat_offs.append(trd[0].astype(np.float))
            self.dat_amps_raw.append(trd[1].astype(np.float))
            self.dat_amps_subs.append(trd[2].astype(np.float))
            self.dat_residual.append(trd[3].astype(np.float))
            self.dat_u.append(trd[4].astype(np.float))
            self.dat_q.append(trd[5].astype(np.float))

        self.dat_offs = np.array(self.dat_offs)
        self.dat_amps_raw = np.array(self.dat_amps_raw)
        self.dat_amps_subs = np.array(self.dat_amps_subs)
        self.dat_residual = np.array(self.dat_residual)
        self.dat_u = np.array(self.dat_u)
        self.dat_q = np.array(self.dat_q)


    def load_qlt_subscan(self, scannum):

        command = 'grep ' +  "' \<%04d\> '" %scannum  + ' LIST.fit'
        tr = os.popen(command).read()[:-1]
        trfit = tr.split('\n')

        n = len(trfit)
        qlt_lon, qlt_lat = [], []

        for _trfit in trfit:
            if 'ALON' in _trfit:
                if _trfit[0] == '#': qlt_lon.append(False)
                else: qlt_lon.append(True)

            if 'ALAT' in _trfit:
                if _trfit[0] == '#': qlt_lat.append(False)
                else: qlt_lat.append(True)

        return np.array(qlt_lon), np.array(qlt_lat)

    def corr_point(self):

        scannums = os.popen("cat %s/*.fit | grep -v '#' | sort -gk5 | awk '{print $2}' | uniq" %self.cfg.fit_path).read()
        scannums = scannums.split('\n')[:-1]
        scannums = [int(s) for s in scannums]

        head = '#%9s %11s %8s %12s %13s %11s ' \
                '%7s %7s %6s %10s %7s %7s %11s\n#%s\n#' \
                %('MJD','Scan','Source','Amp','Err','Azi','Elv','UT',
                'LST','par.Ang','Tsys','Cal','SunAng','-'*131)
        fmt = '%s %0.8f %5d  %-11s %12.5e %13.5e %9.3f %7.3f %6.2f ' \
                '%6.2f %9.3f %8.3f %9.3f %8.3f'


        if not os.path.isdir(self.proc):
            os.mkdir(self.proc)
        all_src = []

        if self.lnk or not os.path.isfile(self.lfit):
            self.update_lfit()

        f_lflux = open(os.path.join(self.proc, self.lflux), 'w')
        print(head, file=f_lflux)

        for scan in scannums:
            if not self.load_fit(scan): continue
            qlt_lon, qlt_lat = self.load_qlt_subscan(scan)
            n1 = len([e for e in qlt_lon if e == True])
            n2 = len([e for e in qlt_lat if e == True])
            if False in qlt_lon or False in qlt_lat:
                print('%04d %9s %3d %3d xxxxxx' % (scan, self.src, n1, n2))
            else:
                print('%04d %9s %3d %3d' % (scan, self.src, n1, n2))

            qlt_sub = np.append(qlt_lon, qlt_lat)
            qlt_sub = np.array(qlt_sub, dtype=bool)

            mjd = np.mean(self.mjd)
            azi = np.mean(self.azi)
            elv = np.mean(self.elv)
            utc = al.mjd2utc(mjd)
            utc = (utc-int(utc))*24
            lst = al.mjd2lst(mjd, self.cfg.obs_lon)
            pa = np.mean(self.pa)
            if True in qlt_sub:
                tsys = np.mean(self.tsys[qlt_sub])
            else:
                tsys = np.mean(self.tsys)
            cal = np.mean(self.cal)
            sa = np.mean(self.sunang)
            amp, aerr, tmp, tmp, tmp, tmp, qlt = self._point_normal()
            del tmp
            if qlt == True: mark = ' '
            else: mark = '#'

            all_src.append(self.src)

            print(fmt \
            %(mark, mjd, scan, self.src, self.tcal*amp, self.tcal*aerr,
            azi, elv, utc, lst, pa, tsys, cal, sa), file=f_lflux)

        f_lflux.close()

        all_src = list(set(all_src))
        all_src.sort()

        self.cfg.all_src = all_src
        self.cfg.pri_cal = [src for src in all_src if src in sample.PRI_CAL]
        self.cfg.sec_cal = [src for src in all_src if src in sample.SEC_CAL]
        self.cfg.target = [src for src in all_src if src not in
                       self.cfg.pri_cal + self.cfg.sec_cal]

        self.cfg.gain_cal = self.cfg.pri_cal + self.cfg.sec_cal
        self.cfg.time_cal = self.cfg.pri_cal + self.cfg.sec_cal
        self.cfg.gain_cal.sort()
        self.cfg.time_cal.sort()

        ### try to update config file ###
        self.cfg.write_sources()


    def update_lfit(self):

        scannums = os.listdir(self.cfg.fit_path)
        scannums = [int(e[0:4]) for e in scannums if e.endswith('.fit')]
        scannums.sort()

        head = '# %7s %8s %7s %6s %7s %12s %9s %9s %11s ' \
                '%12s %5s %7s %7s %11s %6s %8s %11s\n' \
                % ('Source', 'Scan', 'direct', 'Subsc', 'MJD', 'Azi',
                'Elv', 'Amp', 'err','Offset', 'err', 'FWHM', 'err',
                'Parangle', 'Tsys', 'Cal', 'SunAng')
        head = '%s#%s' %(head, '-'*161)

        f_lfit = open(self.lfit, 'w')
        print(head, file=f_lfit)

        if self.flag:
            if not self.cfg.load_quality_rules():
                self.cfg.write_quality_rules()

        for scannum in scannums:
            if not self.load_fit(scannum): continue
            if self.flag == 0:
                mark = [' ']*len(self.amp)
                q1, q2 = self.qlt_ctrl_none()
                q = np.append(q1, q2)
                mark = [e and ' ' or '#' for e in q]
            else:
                q1, q2 = self.qlt_ctrl(0)
                q = np.append(q1, q2)
                mark = [e and ' ' or '#' for e in q]
            print('#', file=f_lfit)
            for i in range(len(self.amp)):
                print('%s %-11s %04d %6s %4d %15.8f ' \
                        '%9.3f %8.3f %12.5e %12.5e %7.2f ' \
                        '%6.2f %8.2f %7.2f %8.3f %8.3f %10.3f %9.3f' \
                        % (mark[i], self.src, self.scannum, self.scandir[i],
                        self.subnum[i], self.mjd[i], self.azi[i], self.elv[i],
                        self.amp[i], self.aerr[i], self.off[i], self.oerr[i],
                        self.hpbw[i], self.herr[i], self.pa[i], self.tsys[i],
                        self.cal[i], self.sunang[i]), file=f_lfit)

        f_lfit.close()


    def qlt_ctrl_none(self):
        """
        Just output subscan information, without quality control.
        """

        q = np.ones(len(self.amp), dtype=bool)
        n1 = len(self.amp[self.scandir=='ALON'])
        n2 = len(self.amp[self.scandir=='ALAT'])
        print('%04d %9s %3d %3d' %(self.scannum, self.src, n1, n2))

        return np.array([1]*n1, dtype=bool), np.array([1]*n1, dtype=bool)
#        out_fmt = '  %6s %4d %12.5e %7.2f %8.2f %8.3f'
#        for i in range(len(self.amp)):
#            print out_fmt \
#                    % (self.scandir[i], self.subnum[i],
#                    self.amp[i], self.off[i],
#                    self.hpbw[i], self.tsys[i])


    def qlt_ctrl(self, mode):
        """
        Do quality control, either with simple or full mode.
        The simple mode calculate every quality flag except the asymmetry of
        the source.
        """

        if mode == 0:
            qlt_lon, qlt_lat = \
            self.check_scan_simple(self.amp, self.aerr,
                    self.off, self.hpbw, self.scandir)

        else:
            qlt_lon, qlt_lat, idc_sub, goodness, qlt = \
            self.check_scan_full(self.amp, self.aerr,
                    self.off, self.hpbw, self.scandir)

            goodness = goodness + np.arange(len(self.amp))*1j
            goodness.sort()
            goodness = [int(e)+1 for e in goodness.imag]

        q = np.append(qlt_lon ,qlt_lat)
        n1 = len(self.amp[self.scandir=='ALON'][qlt_lon])
        n2 = len(self.amp[self.scandir=='ALAT'][qlt_lat])
        if False in qlt_lon or False in qlt_lat:
            print('%04d %9s %3d %3d xxxxxx' % (self.scannum, self.src, n1, n2))
        else:
            print('%04d %9s %3d %3d' % (self.scannum, self.src, n1, n2))

#        out_fmt0 = '  %6s %4d %12.5e %7.2f %8.2f %8.3f'
#        out_fmt1 = '# %6s %4d %12.5e %7.2f %8.2f %8.3f'
#        for i in range(len(self.amp)):
#            if q[i] == True: out_fmt = out_fmt0
#            else: out_fmt = out_fmt1
#            print out_fmt \
#                        % (self.scandir[i], self.subnum[i],
#                           self.amp[i], self.off[i],
#                           self.hpbw[i], self.tsys[i])

        if mode == 0:
            return qlt_lon, qlt_lat

        else:
            return qlt_lon, qlt_lat, idc_sub, goodness, qlt


    def check_subscan(self, a, da, o, h):

        quality = True ### good subscan
        # 0 for good, 1 for OK and 2 for bad.
        indicator = [0, 0, 0]

        if (a <= 0) or (da/a > self.cfg.rerr):
            quality = False
            indicator[0] = 2
        elif da/a > self.cfg.rerr*0.5:
            indicator[0] = 1

        if np.fabs(o) > self.cfg.doff:
            quality = False
            indicator[1] = 2
        elif np.fabs(o) > self.cfg.doff*0.5:
            indicator[1] = 1

        if np.fabs(h - self.cfg.beam) > self.cfg.beam*self.cfg.dwidth:
            quality = False
            indicator[2] = 2
        elif np.fabs(h - self.cfg.beam) > self.cfg.beam*self.cfg.dwidth*0.5:
            indicator[2] = 1

        goodness = np.fabs(da/a)*(np.fabs(o)/self.cfg.beam)*\
                np.fabs(h/self.cfg.beam-1)

        return quality, indicator, goodness


    def check_outlier(self, a):
        """
        function that remove outlier subscans in one direction,
        caution since it may remove extremely fast variability
        """
        N = len(a)
        quality = np.ones(N, dtype=bool)
        b = a*1 + np.arange(N)*1j
        bb = copy.deepcopy(b)
        bb.sort()

        while (N > 2):
            c = bb.real
            n = len(bb)
            d = np.array(bb.imag, dtype=int).tolist()
            d_cmin = np.fabs(c[0]*n - np.sum(c))/(np.sum(c) - c[0])
            d_cmax = np.fabs(c[-1]*n - np.sum(c))/(np.sum(c) - c[-1])

            if (d_cmin >= d_cmax) and (d_cmin > self.cfg.davg):
                quality[d[0]] = False
            elif (d_cmin <= d_cmax) and (d_cmax > self.cfg.davg):
                quality[d[-1]] = False
            else: break

            bb = b[quality]
            N = len(bb)

        return quality


    def check_scan_simple(self, an, dan, on, hn, dirs):

        flt_lon = dirs == 'ALON'
        flt_lat = dirs == 'ALAT'
        nlon = np.sum(flt_lon*1)
        nlat = np.sum(flt_lat*1)
        n = nlon + nlat
        qlt_sub = np.ones(n, dtype=bool)
        for i in range(n):
            qlt_sub[i], tmp0, tmp1 = \
                self.check_subscan(an[i], dan[i], on[i], hn[i])

        del tmp0, tmp1

        qlt_lon = qlt_sub[flt_lon]
        qlt_lat = qlt_sub[flt_lat]
        idx_lon = np.arange(nlon)[qlt_lon==True]
        idx_lat = np.arange(nlat)[qlt_lat==True]

        ##### caution #####
        if True in qlt_lon:
            tmp_qlt_lon = self.check_outlier(an[flt_lon][qlt_lon==True])
            idx_lon = idx_lon[tmp_qlt_lon]
        if True in qlt_lat:
            tmp_qlt_lat = self.check_outlier(an[flt_lat][qlt_lat==True])
            idx_lat = idx_lat[tmp_qlt_lat]
        #####

        qlt_lon = np.zeros(nlon, dtype=bool)
        qlt_lat = np.zeros(nlat, dtype=bool)
        for i in idx_lon: qlt_lon[i] = True
        for i in idx_lat: qlt_lat[i] = True

        return qlt_lon, qlt_lat


    def check_scan_full(self, an, dan, on, hn, dirs):

        qlt = True
        flt_lon = dirs == 'ALON'
        flt_lat = dirs == 'ALAT'
        nlon = np.sum(flt_lon*1)
        nlat = np.sum(flt_lat*1)
        n = nlon + nlat
        qlt_sub = np.ones(n, dtype=bool)
        idc_sub = [[0, 0, 0]]*n
        goodness = [0]*n
        for i in range(n):
            qlt_sub[i], idc_sub[i], goodness[i] = \
                self.check_subscan(an[i], dan[i], on[i], hn[i])

        qlt_lon = qlt_sub[flt_lon]
        qlt_lat = qlt_sub[flt_lat]
        idx_lon = np.arange(nlon)[qlt_lon==True]
        idx_lat = np.arange(nlat)[qlt_lat==True]

        if (not True in qlt_lon) or (not True in qlt_lat):
            qlt = False

        ##### caution #####
        else:
            tmp_qlt_lon = self.check_outlier(an[flt_lon][qlt_lon==True])
            idx_lon = idx_lon[tmp_qlt_lon]
            tmp_qlt_lat = self.check_outlier(an[flt_lat][qlt_lat==True])
            idx_lat = idx_lat[tmp_qlt_lat]
            #####

            qlt_lon = np.zeros(nlon, dtype=bool)
            qlt_lat = np.zeros(nlat, dtype=bool)
            for i in idx_lon: qlt_lon[i] = True
            for i in idx_lat: qlt_lat[i] = True

            an_lon = an[flt_lon][qlt_lon]
            an_lat = an[flt_lat][qlt_lat]
            dan_lon = dan[flt_lon][qlt_lon]
            dan_lat = dan[flt_lat][qlt_lat]

            if (not True in qlt_lon) or (not True in qlt_lat):
                qlt = False
            else:
                a_lon = np.average(an_lon, weights = dan_lon**-2.0)
                a_lat = np.average(an_lat, weights = dan_lat**-2.0)

                if np.fabs(a_lon - a_lat)/(a_lon + a_lat)*2 > \
                        self.cfg.dsymm: qlt = False

        return qlt_lon, qlt_lat, idc_sub, goodness, qlt


    def _point_all_good(self, amp, aerr, off, oerr, hpbw, herr,
                        nlon, nlat):
        """
        for remaining good subscans (bad are already flagged)
        """

        amp_lon = amp[0:nlon]
        off_lon = off[0:nlon]
        hpbw_lon = hpbw[0:nlon]
        aerr_lon = aerr[0:nlon]
        oerr_lon = oerr[0:nlon]
        herr_lon = herr[0:nlon]

        amp_lat = amp[nlon:nlon+nlat]
        off_lat = off[nlon:nlon+nlat]
        hpbw_lat = hpbw[nlon:nlon+nlat]
        aerr_lat = aerr[nlon:nlon+nlat]
        oerr_lat = oerr[nlon:nlon+nlat]
        herr_lat = herr[nlon:nlon+nlat]

        wt_amp_lon = aerr_lon**(-2.0)
        wt_hpbw_lon = herr_lon**(-2.0)
        wt_amp_lat = aerr_lat**(-2.0)
        wt_hpbw_lat = herr_lat**(-2.0)

        aver_amp_lon = np.average(amp_lon, weights = wt_amp_lon)
        aver_off_lon = np.average(off_lon, weights = wt_amp_lon)
        aver_hpbw_lon = np.average(hpbw_lon, weights = wt_amp_lon)

        aver_amp_lat = np.average(amp_lat, weights = wt_amp_lat)
        aver_off_lat = np.average(off_lat, weights = wt_amp_lat)
        aver_hpbw_lat = np.average(hpbw_lat, weights = wt_amp_lat)

        aver_aerr_lon = (np.sum(wt_amp_lon))**(-0.5)
        aver_aerr_lat = (np.sum(wt_amp_lat))**(-0.5)
        aver_herr_lon = (np.sum(wt_hpbw_lon))**(-0.5)
        aver_herr_lat = (np.sum(wt_hpbw_lat))**(-0.5)

        # in case of just one good subscan remaining in one direction
        if nlon == 1: aver_oerr_lon = oerr_lon[0]
        else: aver_oerr_lon = np.std(oerr_lon)
        if nlat == 1: aver_oerr_lat = oerr_lat[0]
        else: aver_oerr_lat = np.std(oerr_lat)

        if not self.aver_tag:
            f_lon = np.exp(4.0*np.log(2.0)*(aver_off_lat/aver_hpbw_lat)**2.0)
            f_lat = np.exp(4.0*np.log(2.0)*(aver_off_lon/aver_hpbw_lon)**2.0)
            aver_amp_lon *= f_lon
            aver_amp_lat *= f_lat

            aver_aerr_lon = ((aver_aerr_lon*f_lon)**2.0 +
                            (aver_amp_lon*8.0*np.log(2.0)*aver_off_lat
                            /aver_hpbw_lat**2.0*aver_oerr_lat)**2.0 +
                            (-aver_amp_lon*8.0*np.log(2.0)*aver_off_lat**2.0
                            /aver_hpbw_lat**3.0*aver_herr_lat)**2.0)**0.5

            aver_aerr_lat = ((aver_aerr_lat*f_lat)**2.0 +
                            (aver_amp_lat*8.0*np.log(2.0)*aver_off_lon
                            /aver_hpbw_lon**2.0*aver_oerr_lon)**2.0 +
                            (-aver_amp_lat*8.0*np.log(2.0)*aver_off_lon**2.0
                            /aver_hpbw_lon**3.0*aver_herr_lon)**2.0)**0.5

        # check the source symmetry
        quality = True
        if self.flag == 1:
            dsym = np.fabs(aver_amp_lon - aver_amp_lat) \
                    /(aver_amp_lon + aver_amp_lat)*2
            if dsym > self.cfg.dsymm: quality = False

        # natural weighting
        aver_amp = np.array([aver_amp_lon, aver_amp_lat])
        aver_off = np.array([aver_off_lon, aver_off_lat])
        aver_hpbw = np.array([aver_hpbw_lon, aver_hpbw_lat])

        wt_aver_amp = np.array([aver_aerr_lon, aver_aerr_lat])**(-2.0)
        wt_aver_off = np.array([aver_oerr_lon, aver_oerr_lat])**(-2.0)
        wt_aver_hpbw = np.array([aver_herr_lon, aver_herr_lat])**(-2.0)

        final_amp = np.average(aver_amp, weights=wt_aver_amp)
        final_off = np.average(aver_off, weights=wt_aver_off)
        final_hpbw = np.average(aver_hpbw, weights=wt_aver_hpbw)

        final_aerr = np.sum(wt_aver_amp)**(-0.5)
        final_oerr = np.sum(wt_aver_off)**(-0.5)
        final_herr = np.sum(wt_aver_hpbw)**(-0.5)

        return (final_amp, final_aerr, final_off, final_oerr,
                final_hpbw, final_herr, quality)

    def _point_normal(self):

        qutlity = True

        flt_lon = self.scandir == 'ALON'
        flt_lat = self.scandir == 'ALAT'

        if (not True in flt_lon) or (not True in flt_lat):
            if self.flag == 0:
                qlt = np.ones(len(self.scandir), dtype=bool)
            else:
                qlt_lon, qlt_lat = self.load_qlt_subscan(self.scannum)
                qlt = np.append(qlt_lon, qlt_lat)
                qlt = qlt.astype(bool)

            if not True in qlt:
                qlt = np.ones(len(qlt), dtype=bool)

            amp = np.array(list(self.amp[qlt])*2)
            aerr = np.array(list(self.aerr[qlt])*2)
            off = np.array(list(self.off[qlt])*2)
            oerr = np.array(list(self.oerr[qlt])*2)
            hpbw = np.array(list(self.hpbw[qlt])*2)
            herr = np.array(list(self.herr[qlt])*2)

            nlon = nlat = len(self.scandir[qlt])

            (final_amp, final_aerr, final_off, final_oerr,
            final_hpbw, final_herr, quality) =  \
            self._point_all_good(amp, aerr, off, oerr, hpbw, herr,
                    nlon, nlat)
            quality = False

            return (final_amp, final_aerr, final_off, final_oerr,
                    final_hpbw, final_herr, quality)

        if self.flag == 0:
            amp_lon = self.amp[flt_lon]
            amp_lat = self.amp[flt_lat]
            qlt_lon = np.ones(len(amp_lon), dtype=bool)
            qlt_lat = np.ones(len(amp_lat), dtype=bool)

        else:
            qlt_lon, qlt_lat = self.load_qlt_subscan(self.scannum)
            amp_lon = self.amp[flt_lon][qlt_lon]
            amp_lat = self.amp[flt_lat][qlt_lat]

        nlon, nlat = len(amp_lon), len(amp_lat)
        qlt = np.append(qlt_lon, qlt_lat)

        if max(nlon, nlat) == 0:
            return -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, False

        elif max(nlon, nlat) !=0 and min(nlon, nlat)==0:
            nlon = max(nlon, nlat)
            nlat = nlon*1

            amp = np.array(list(self.amp[qlt])*2)
            aerr = np.array(list(self.aerr[qlt])*2)
            off = np.array(list(self.off[qlt])*2)
            oerr = np.array(list(self.oerr[qlt])*2)
            hpbw = np.array(list(self.hpbw[qlt])*2)
            herr = np.array(list(self.herr[qlt])*2)

            (final_amp, final_aerr, final_off, final_oerr,
            final_hpbw, final_herr, quality) =  \
            self._point_all_good(amp, aerr, off, oerr, hpbw, herr,
                    nlon, nlat)
            quality = False

            return (final_amp, final_aerr, final_off, final_oerr,
                    final_hpbw, final_herr, quality)

        else:
            amp = self.amp[qlt]
            aerr = self.aerr[qlt]
            off = self.off[qlt]
            oerr = self.oerr[qlt]
            hpbw = self.hpbw[qlt]
            herr = self.herr[qlt]

            return self._point_all_good(amp, aerr, off, oerr, hpbw, herr,
                    nlon, nlat)


class ListFlux(object):

    def __init__(self, cfg, proc0, proc1):

        self.cfg = cfg
        self.proc0 = proc0
        self.proc1 = proc1
        self.lflux = 'LIST.flux'
        self.head = '#%9s %11s %8s %12s %13s %11s ' \
                    '%7s %7s %6s %10s %7s %7s %11s\n#%s\n#' \
                    %('MJD','Scan','Source','Amp','Err','Azi','Elv','UT',
                    'LST','par.Ang','Tsys','Cal','SunAng','-'*131)
        self.fmt = '%s %0.8f %5d  %-11s %12.5e %13.5e %9.3f %7.3f %6.2f ' \
                    '%6.2f %9.3f %8.3f %9.3f %8.3f'


    def split_listflux(self):

        self.proc0 = self.proc1 # working on calibrated flux
        self.load_listflux()

        if 'flux' in self.proc1:
            for src in self.cfg.all_src:
                flt = self.src == src
                data = [e[flt] for e in self.data]
                mean_f = np.mean(data[4])
                mean_err = np.std(data[4])
                wt = data[5]**-2.0
                avg_f = np.average(data[4], weights=wt)
                avg_err = np.sum(wt)**-0.5

                addhdr = '# %-15s %4s %6.3f %2s %6.3f %s' \
                        %(src, '<S>=', mean_f, '+-', mean_err, 'Jy')
                addhdr = '%s\n# %20s %6.3f %2s %6.3f %s' \
                        %(addhdr, '<S>(gew)=', avg_f, '+-', avg_err, 'Jy')
                hdr = '%s\n%s' %(addhdr, self.head)
                writecol(os.path.join(self.proc1, 'FLUX.%s' %src),
                        hdr, None, data, self.fmt)

        else:
            for src in self.cfg.all_src:
                flt = self.src == src
                data = [e[flt] for e in self.data]
                writecol(os.path.join(self.proc1, 'FLUX.%s' %src),
                        self.head, None, data, self.fmt)


    def load_listflux(self):

        data = readcol(os.path.join(self.proc0, self.lflux),
                start=3, comment='#', flag=False, fmt='fdsffffffffff')
        self.data = [np.array(e) for e in data]
        del data

        self.qlt = np.array([e==' ' and True or False for e in self.data[0]])
#        self.mjd = self.data[1]
#        self.scan = self.data[2]
        self.src = self.data[3]
#        self.amp = self.data[4]
#        self.aerr = self.data[5]
#        self.azi = self.data[6]
#        self.elv = self.data[7]
#        self.utc = self.data[8]
#        self.lst = self.data[9]
#        self.pa = self.data[10]
#        self.tsys = self.data[11]
#        self.cal = self.data[12]
#        self.sunang = self.data[13]


    def write_listflux(self):

        if not os.path.isdir(self.proc1):
            os.mkdir(self.proc1)

        writecol(os.path.join(self.proc1, self.lflux),
                self.head, None, self.data, self.fmt)


    def load_qlt_scan(self, scannum):

        flt = self.data[2] == scannum
        return self.qlt[flt][0]


    def compute_tau(self, tgr, tcal):

        self.cfg.load_telescope()
#        self.load_listflux()
#        elv, tsys = np.array(self.data[7]), np.array(self.data[11])*tcal
#        elv, tsys = elv[self.qlt], tsys[self.qlt]
#        am = 1.0/np.sin(np.deg2rad(elv)) # airmass
        if tgr != 99:
            elv, tsys, temp = readcol('temperature.dat', cols=[3,4,5], fmt='fff')
            elv, tsys, temp = np.array(elv), np.array(tsys), np.array(temp)
            tsys = tsys*tcal
            tatm = (temp+273.15)*1.12 - 50
        else:
            self.load_listflux()
            elv, tsys = np.array(self.data[7]), np.array(self.data[11])*tcal
            tatm = (tgr + 273.15)*1.12 - 50

        am = 1.0/np.sin(np.deg2rad(elv)) # airmass
        AM = am*tatm

#        AM_low, tsys_low = ml.lower_envelope(AM, tsys)
#        func = np.polyfit(AM_low, tsys_low, 1)
        func = ml.lower_envelope(AM[self.qlt], tsys[self.qlt])

        self.cfg.Tgr = np.mean(temp)
        self.cfg.Trec = func[1]
        self.cfg.Tau0 = func[0]

        xmin = min(elv) - 0.1*(max(elv) - min(elv))
        xmax = max(elv) + 0.1*(max(elv) - min(elv))
        ymin = min(tsys) - 0.1*(max(tsys) - min(tsys))
        ymax = max(tsys) + 0.1*(max(tsys) - min(tsys))

        return func, AM



    def write_tau(self):

        self.cfg.write_opacity()


    def load_tau(self):

        return self.cfg.load_opacity()


    def corr_tau(self, tele='eff', off=0):

        ldf = self.load_tau()
        if not ldf:
            self.compute_tau(0, 1)
            self.write_tau()
            print('Loading Tau failed, I will calculate it automatically.')
        else:
            self.load_listflux()

        # read wvr data (effelsberg)
        if tele == 'eff':
            corr = readcol('weather.dat', fmt='fdsfffffff', flag=True)
            corr_scan = corr[1]
            corr_val = []
#            for s in self.data[2]:
            for i in range(len(self.data[2])):
                scan = self.data[2][i]
                try:
                    idx = corr_scan.index(scan)
                    val = np.exp((corr[6][idx]+off)/np.sin(np.deg2rad(self.data[7][i])))
#                    corr_val.append(corr[7][idx])
                    corr_val.append(val)
                except:
                    corr_val.append(1)
                    print('bad scan %d' %scan)
            self.data[4] *= corr_val
            self.data[5] *= corr_val
        else:
            temp = readcol('temperature.dat', cols=[5], fmt='f')[0]
            temp = np.array(temp)
            elv = self.data[7]
            am = 1.0/np.sin(np.deg2rad(elv))
            tatm = (temp + 273.15)*1.12 - 50
            tsys = np.array(self.data[11])
            print(len(tsys), len(tatm))
#            corr_val0 = 1.0 - (tsys - self.cfg.Trec)/tatm
            corr_val = np.exp(1.0/am*np.log(1-(tsys-self.cfg.Trec)/tatm))
#            print(corr_val/corr_val0-1)
#            print(corr_val)
            self.data[4] /= corr_val
            self.data[5] /= corr_val

        self.write_listflux()
        self.split_listflux()


    def compute_gain(self, weight=0, order=2, limit=0, peak0=45):

        self.cfg.load_sources()
        self.load_listflux()
        has_outlier = False

        while 1:

            gain_cal = [src for src in self.cfg.gain_cal \
                        if src in self.src[self.qlt]]

            _x, _y, _yerr = np.array([]), np.array([]), np.array([])
            _src, _scan = np.array([]), np.array([])
            peaks = []

            for gcal in gain_cal:
                flt = self.src == gcal
                flt *= self.qlt

                gelv = self.data[7][flt]
                gamp = self.data[4][flt]
                gaerr = self.data[5][flt]

                if not weight:
                    gpar, gperr, sigma = ml.fit_poly(gelv, gamp, 0, order)
                else:
                    gpar, gperr, sigma = ml.fit_poly(gelv, gamp, gaerr, order)

                pmax = np.roots(np.polyder(gpar))
                pmax = pmax[(pmax>0)*(pmax<90)*(pmax>peak0-15)*(pmax<peak0+15)]
                if not len(pmax):
                    print('Warning: Cannot find an optimized elv. for %s.' \
                            %gcal)
                    pmax = np.median(gamp) + max(gamp)
                    pmax *= 0.5
                    pmax = [pmax]
                    if len(gain_cal) > 1: break
                ymax = max(np.polyval(gpar, pmax))
                gamp /= ymax
                gaerr /= ymax

                _x = np.append(_x, gelv)
                _y = np.append(_y, gamp)
                _yerr = np.append(_yerr, gaerr)
                _src = np.append(_src, self.src[flt])
                _scan = np.append(_scan, self.data[2][flt])
                peaks.append(ymax)

            if not weight:
                gpar, gperr, sigma = ml.fit_poly(_x, _y, 0, order)
            else:
                gpar, gperr, sigma = ml.fit_poly(_x, _y, _yerr, order)

            if limit:
                #npt = len(_x)
                dy = np.fabs(_y-np.polyval(gpar, _x))
                idx = dy.argmax()
                dy = dy[idx]
                if dy > limit*sigma:
                    bad_src = _src[idx]
                    bad_scan = _scan[idx]
                    self.qlt *= self.data[2] != bad_scan
                    print('|==> Scan Removed: SCAN= %04d,   SOURCE= %9s,    '\
                            'ELV= %0.2f,   AMP= %0.4f' \
                            %(bad_scan, bad_src, _x[idx], _y[idx]))
                    has_outlier = True

                else:
                    has_outlier = False
                    break

            else:
                break

        plt_xmin = min(_x) - 0.1*(max(_x) - min(_x))
        plt_xmax = max(_x) + 0.1*(max(_x) - min(_x))
        plt_ymin = min(_y-_yerr) - 0.1*(max(_y+_yerr) - min(_y-_yerr))
        plt_ymax = max(_y+_yerr) + 0.1*(max(_y+_yerr) - min(_y-_yerr))

        pmax = np.roots(np.polyder(gpar))
        pmax = pmax[(pmax>0)*(pmax<90)]
        if not len(pmax):
            print('Warning: Cannot find an optimized elv.')
            pmax = np.median(gamp) + max(gamp)
            pmax *= 0.5
            pmax = [pmax]
        ymax = max(np.polyval(gpar, pmax))
        gamp /= ymax
        gaerr /= ymax

        self.gain = gpar/ymax
        self.gain_err = gperr/ymax

        peaks = np.array(peaks)*ymax

        return gain_cal, peaks, sigma, plt_xmin, plt_xmax, plt_ymin, plt_ymax


    def write_gain(self):

        text_c = ''
        text_c_err = ''

        for i in range(len(self.gain)):
            text_c += '%-+0.5e  ' %(self.gain[i])
            text_c_err += '%-+0.5e  ' %(self.gain_err[i])

        fname = '%d.gain.curve' %(int(self.proc1.split('.')[0]))
        with open(fname, 'w') as f:
            print('#  GAIN-ELEVATION CURVE', file=f)
            print('#' + '-'*56, file=f)
            print('#', file=f)
            print('   COEFFICIENT = ' + text_c, file=f)
            print('   ERROR       = ' + text_c_err, file=f)

        print('#  GAIN-ELEVATION CURVE')
        print('#' + '-'*56)
        print('#')
        print('   COEFFICIENT = ' + text_c)
        print('   ERROR       = ' + text_c_err)


    def load_gain(self):

        fname = '%d.gain.curve' %(int(self.proc1.split('.')[0]))
        with open(fname) as f:
            trs = f.readlines()
            trs = [e[:-1] for e in trs if e[0]!= "#"]
            self.gain = np.array([float(e) for e in trs[0].split()[2:]])
            self.gain_err = np.array([float(e) for e in trs[1].split()[2:]])


    def corr_gain(self):

#        fname = '%d.gain.curve' %(int(self.proc1.split('.')[0]))
#        if not os.path.isfile(fname):
#            self.compute_gain()
#            self.write_gain()
#        else:
        self.load_gain()
        self.load_listflux()

        corr_val = np.polyval(self.gain, self.data[7])
        self.data[4] /= corr_val
        self.data[5] /= corr_val

        self.write_listflux()
        self.split_listflux()


    def compute_time(self, binsize):

        self.cfg.load_sources()
        self.load_listflux()

        time_cal = [src for src in self.cfg.time_cal \
                    if src in self.src[self.qlt]]

        _x, _y, _yerr = np.array([]), np.array([]), np.array([])

        for tcal in time_cal:
            flt = self.src == tcal
            flt *= self.qlt

            ttime = self.data[1][flt]
            tamp = self.data[4][flt]
            taerr = self.data[5][flt]

            pmean = np.mean(tamp)

            tamp /= pmean
            taerr /= pmean

            _x = np.append(_x, ttime)
            _y = np.append(_y, tamp)
            _yerr = np.append(_yerr, taerr)

        idxs = _x.argsort()

        _x = _x[idxs]
        _y = _y[idxs]
        _yerr = _yerr[idxs]
        ##################### Gaussian Process  ##############################
#        x = np.atleast_2d(_x).T
#        kernel = ConstantKernel(1, (1e-2, 1e2))*RBF(1, (1e-2, 1e2)) \
#                 + RBF(1, (1e-2, 1e2))
#        gp = GaussianProcessRegressor(kernel=kernel, alpha=binsize*1700*_yerr**2)
#        gp.fit(x, _y)
#        self.tcurve, self.terror = gp.predict(np.atleast_2d(self.data[1]).T,
#                                              return_std=True)
        # it turns out that Gaussian process regression is very sensitive
        # on noise (alpha), and the result is not as good as cubic spline
        ######################################################################

        x, y, yerr = np.array([]), np.array([]), np.array([])
        i = 1
        du = max(_x) - min(_x)
        x0 = min(_x)
        while i < du/binsize + 1:
            flt = (_x < x0+i*binsize)*(_x>=x0+(i-1)*binsize)
            if True in flt:
                x = np.append(x, np.mean(_x[flt]))
                y = np.append(y, np.mean(_y[flt]))
                yerr = np.append(yerr, 1 and np.std(_y[flt]) or _yerr[flt])
            i += 1

        self.tck = interpolate.splrep(x, y)
        self.tcurve = interpolate.splev(self.data[1], self.tck)
        self.terror = np.interp(self.data[1], x, yerr)

        xmin = min(_x) - 0.1*(max(_x) - min(_x))
        xmax = max(_x) + 0.1*(max(_x) - min(_x))
        ymin = min(_y-_yerr) - 0.1*(max(_y+_yerr) - min(_y-_yerr))
        ymax = max(_y+_yerr) + 0.1*(max(_y+_yerr) - min(_y-_yerr))

        return time_cal, xmin, xmax, ymin, ymax


    def write_time(self):

        # write result into time.curve file
        fname = '%d.time.curve' %(int(self.proc1.split('.')[0]))
        header = '%s\n%s\n%s' %('#  GAIN-TIME CURVE', '#'+'-'*51, '#')
        header = '%s\n%6s %24s %17s' %(header, '#  MJD', 'FACTOR', 'ERR')

        writecol(fname, header, None, [self.data[1], self.tcurve, self.terror],
                '   %0.8f %17.8f %16.4f')


    def load_time(self):

        fname = '%d.time.curve' %(int(self.proc1.split('.')[0]))

        data = readcol(fname)
        self.tcurve = np.array(data[1])
        self.terror = np.array(data[2])


    def corr_time(self, scale):

        self.load_time()
        self.load_listflux()

        self.data[5] = self.data[4]/self.tcurve* \
                        ((self.data[5]/self.data[4])**2 +
                        (scale*self.terror/self.tcurve)**2)**0.5
        self.data[4] /= self.tcurve

        self.write_listflux()
        self.split_listflux()


    def compute_dpfu(self):

        self.cfg.load_sources()
        self.load_listflux()
        freq = self.cfg.freq

        pri_cal = [src for src in self.cfg.pri_cal \
                    if src in self.src[self.qlt]]

        src = np.array([])
        mjd, amp, aerr, elv, scannum = np.array([]), np.array([]), \
                np.array([]), np.array([]), np.array([])

        factor = np.array([])

        for pcal in pri_cal:
            flt = self.src == pcal
            flt *= self.qlt

            absf = al.cal_flux(pcal, freq)
            mjd = np.append(mjd, self.data[1][flt])
            scannum = np.append(scannum, self.data[2])
            src = np.append(src, self.src[flt])
            amp = np.append(amp, self.data[4][flt])
            aerr = np.append(aerr, self.data[5][flt])
            elv = np.append(elv, self.data[7])
            factor = np.append(factor, absf/self.data[4][flt])

            print('%8s: %6.3f %8.4f' %(pcal, absf, np.mean(absf/self.data[4][flt])))

        self.factor = np.mean(factor)
        self.factor_err = np.std(factor)
        self.sensit = np.mean(factor**-1.0)
        self.sensit_err = np.std(factor**-1.0)

        idxs = mjd.argsort()
        mjd = mjd[idxs]
        scannum = scannum[idxs]
        src = src[idxs]
        amp = amp[idxs]
        aerr = aerr[idxs]
        elv = elv[idxs]
        factor = factor[idxs]

        hdr = '# %3s %8s %10s %10s %6s %5s %13s' \
                %('MJD', 'Scan', 'Calibrator', 'Flux',
                'Err', 'Elv', 'Cal-Factor')
        hdr = '%s\n%s' %(hdr, '#'+'-'*62)

        tail = '%s %23s %8.5E %s %7.5E %5s' \
                %('#', 'average Cal.-factor =',
                self.factor, '+-', self.factor_err, 'Jy/K')
        tail = '%s\n%s %23s %8.5E %s %7.5E %5s' \
                %(tail, '#', 'average Sensitivity =',
                self.sensit, '+-', self.sensit_err, 'K/Jy')
        fmt = '%0.3f %4d  %-13s %7.3f %6.3f %5.1f %11.5f'
        writecol('cal.factor', hdr, tail,
                [mjd, scannum, src, amp, aerr, elv, factor], fmt)


    def load_dpfu(self):

        with open('cal.factor') as f:
            trs = f.readlines()
            factor = [e[:-1] for e in trs if 'Jy/K' in e][0]
            sensit = [e[:-1] for e in trs if 'K/Jy' in e][0]

        self.factor = float(factor.split('+-')[0].split()[-1])
        self.factor_err = float(factor.split('+-')[1].split()[0])
        self.sensit = float(sensit.split('+-')[0].split()[-1])
        self.sensit_err = float(sensit.split('+-')[1].split()[0])


    def corr_dpfu(self, scale):

        self.load_dpfu()
        self.load_listflux()

        self.data[5] = self.data[4]*self.factor* \
                        ((self.data[5]/self.data[4])**2 + \
                        (scale*self.factor_err/self.factor)**2)**0.5
        self.data[4] = self.data[4]*self.factor

        self.write_listflux()
        self.split_listflux()


    def get_m0(self, m0=0):

        stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)

        self.load_listflux()
        srcs = list(set(self.src))
        srcs.sort()

        if not m0:
            cals = list(set(self.cfg.pri_cal+\
                    self.cfg.gain_cal+\
                    self.cfg.time_cal))
            cal_amps = np.array([])
            for cal in cals:
                flt = (self.src==cal)*self.qlt
                dd = self.data[4][flt]
                dd /= np.mean(dd)
                cal_amps = np.append(cal_amps, dd)

            m0 = np.std(cal_amps)

        ut0 = al.mjd2utc(self.data[1][0])
        ut1 = al.mjd2utc(self.data[1][-1])
        dhr = int((ut1 - ut0)*24)
        ut0 = str(ut0)[:-3]
        ut1 = str(ut1)[:-3]
        eff = np.sum(self.qlt)/float(len(self.qlt))*100

        fm0 = open('mod.index', 'w')
        print('# %s - %s  %13s%3.0f%s %6s %d hrs %7s %0.2f' \
            % (ut0, ut1, 'Efficiency:', eff, '%', 'T =', dhr, 'm0 =',100*m0))
        print('#' + '-'*92)
        print('%-13s %5s %7s %6s %7s %7s %8s %10s %12s %9s' \
            %('# Source', 'Scans', '<S>[Jy]', 'sigma', 'm[%]', 'Y[%]', \
            'Chi^2', 'red.Chi^2', 'Probability', 'C.I.'))
        print('#' + '-'*92)

        print('# %s - %s  %13s%3.0f%s %6s %d hrs %7s %0.2f' \
            % (ut0, ut1, 'Efficiency:', eff, '%', 'T =', dhr, 'm0 =',100*m0), file=fm0)
        print('#' + '-'*92, file=fm0)
        print('%-13s %5s %7s %6s %7s %7s %8s %10s %12s %9s' \
            %('# Source', 'Scans', '<S>[Jy]', 'sigma', 'm[%]', 'Y[%]', \
            'Chi^2', 'red.Chi^2', 'Probability', 'C.I.'), file=fm0)
        print('#' + '-'*92, file=fm0)

        for src in srcs:

            flt0 = self.src == src
            flt1 = flt0*self.qlt

            tmp = self.data[4][flt1]
            if not len(tmp):
                count = 0
                flt = flt0*(self.data[4]>0)
                if not True in flt:
                    flt = flt0
            else:
                flt = flt1
                count  = len(tmp)

            del tmp
            amps = self.data[4][flt]
            errs = self.data[5][flt]
            flux = np.mean(amps)
            rms = np.std(amps)
            m = rms / flux

            if m > m0: Y = 3*(m*m - m0*m0)**0.5
            else: Y = 0

            dof = len(amps) - 1
            x2 = np.sum(((amps-flux)/errs)**2)
            xr2 = x2 / dof
            prob = stats.chisqprob(x2, dof)
            ci = stats.t.ppf((2-prob)/2, dof)

            print('%-13s %3d %8.3f %7.3f %7.2f %7.2f %9.3f %7.3f %13.4E %10.2E'\
                % (src, count, flux, rms, 100*m, 100*Y, x2, xr2, prob, ci))
            print('%-13s %3d %8.3f %7.3f %7.2f %7.2f'\
                    ' %9.3f %7.3f %13.4E %10.2E'\
                % (src, count, flux, rms, 100*m, 100*Y, x2, xr2, prob, ci), file=fm0)

        fm0.close()



    def performance(self):
        pass


if __name__ == '__main__':
    pass
