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
from cmd import *
from ppgplot import *
import numpy as np
from scipy import interpolate, stats
import ObsData
import Calibration
import MathLib as ml
import AstroLib as al
import perform



class DRP(Cmd):

    def __init__(self):

        Cmd.__init__(self)
        self.prompt = '|-=DrCont@NSRT=->> '
        self.all_pid = []
        self.gs_pid = 0
        self.window = 0
        self.load_cfg()

    def parm_exact(self, line, key):

            return line.split('-'+key)[1].split()[0]

    def load_cfg(self):

        self.cfg = ObsData.Config()
        if not self.cfg.load_path():
            self.cfg.write_path()

        if not self.cfg.load_telescope():
            self.cfg.write_telescope()

        if not self.cfg.load_quality_rules():
            self.cfg.write_quality_rules()


    def all_procs(self):

        procs = os.listdir('./')
        procs = [e for e in procs if e[0].isdigit() and len(e)<=7]
        procs.sort()
        return procs

    def flag_all_listflux(self, idx):

        all_procs = self.all_procs()
        print('#----------')
        for proc in all_procs:
            clf = Calibration.ListFlux(self.cfg, proc, proc)
            clf.load_listflux()
            if clf.data[0][idx] == ' ':
                clf.data[0][idx] = '#'
                print("  Disable Point: %s: %04d in %s" %(clf.src[idx],
                        clf.data[2][idx], proc))
            else:
                clf.data[0][idx] = ' '
                print("  Enable Point: %s: %04d in %s" %(clf.src[idx],
                        clf.data[2][idx], proc))
            clf.write_listflux()
            del clf


    def do_exit(self, arg):
        return 1


    def do_EOF(self, arg):
        print('')
        return 1

    def emptyline(self):
        pass

    def do_cd(self, arg):
        curpth = os.getcwd()
        if not arg:
            env = os.environ
            pth = env["HOME"]
        elif arg[0] == "/": pth = arg
        elif arg == "..": pth = os.path.split(curpth)[0]
        else: pth = curpth + "/" + arg
        if not os.path.isdir(pth):
            print("Error: NO such directory.")
        else: os.chdir(pth)


    def do_pwd(self, arg):
        print(os.getcwd())


    def do_ls(self, arg):
        print(os.popen("ls "+arg).read()[:-1])


    def do_set_window(self, arg):

        window = arg.split()
        self.window = [float(w) for w in window]


    def do_openX(self, arg):
        """
        Open a pgplot device.

        Examples:
            openX -- open a xw window, the same as "openX xw".
            openX ps -- open a ps device.
            openX cps -- open a color ps device.
            openX ps  fname 7 4 -- open a ps device with width=7 and height=4.
            openX cps fname 7 4 -- open a color ps device with width=7
                                   and height=4.
        """
        arg = arg.split()
        if not arg: arg.append("xw")
        device = "/" + arg[0]

        if device.upper() =="/PS" or device.upper() == "/CPS":
            if len(arg) == 1: pid = pgopen(device)
            if len(arg) >= 2:
                filename = arg[1]
                pid = pgopen(filename+device)
                if len(arg) >= 4:
                    width = float(arg[2])
                    height = float(arg[3])
                    height/=width
                    pgpap(width, height)

        else: pid = pgopen(device)
        self.all_pid.append(pid)


    def do_slctX(self, arg):
        """
        Select a specified pgplot device. If the device is not open,
        skip to current one and raise a warning.
        """
        if not arg: pgslct(1)
        else:
            pid = int(arg)
            pgslct(pid)

    def do_which(self, arg):
        """
        Find out which pgplot device is currently in use.
        """
        print(pgqid())

    def do_closeX(self, arg):
        """
        closeX Agrgument

        Close pgplot device(s).

        Arguments:
            It accept either an integer, or a string "all", or None.
            <1> if integer, close the specified pgplot device.
            <2> if "all", close all pgplot device.
            <3> if None, close currently selected pgplot device.

        Examples:
            closeX 3
        """
        if not arg:
            pid =  pgqid()
            pgask(0)
            pgclos(pid)
            if pid == self.gs_pid: self.gs_pid = 0
            self.all_pid.remove(pid)
        elif arg == "all":
            pgask(0)
            pgend()
        else:
            pid = int(arg)
            pgask(0)
            pgclos(pid)
            if pid == self.gs_pid: self.gs_pid = 0
            self.all_pid.remove(pid)


    def do_chk_fit(self, line):
        """
        Check the Gaussian fitting conditions in a pgplot device.
        The argument should be an integer.
        Example: ckfit 123
        """
        scannum = int(line)
        cf = Calibration.FitDat(self.cfg, 0, 1)
        cf.load_fit(scannum)
        cf.load_dat(scannum)
        if not self.window:
            print('using default window 1.35 1.35')
            window = [1.35, 1.35]
        else:
            window = self.window

        subs = len(cf.mjd)
        n1 = int((subs+1)**0.5)
        n2 = n1*1
        while n1*n2 < subs+1:
            n2 = n2 +1
        flag1, flag2 = cf.load_qlt_subscan(scannum)
        sub_flag = np.append(flag1, flag2)
        del flag1, flag2

        pid = pgqid()
        pgslct(pid)
        pgask(0)
        pgpage()
        pgsubp(n1, n2)

        if pgqinf('device') == '':
            color0 = 3
            color1 = 7
            color2 = 2
        else:
            color0 = 4
            color1 = 8
            color2 = 2

        for i in range(subs):
            dataX = cf.dat_offs[i]
            dataY = cf.dat_amps_subs[i]
            if cf.amp[i] > 0:
                fitY = ml.Gauss(1, 1).fGauss(
                        np.array([cf.amp[i], cf.off[i], cf.hpbw[i], 0, 0]),
                        dataX)
            else:
                fitY = dataY*1.0

            ampcolor, offcolor, hpbwcolor = color0, color0, color0

            xmin, xmax = min(dataX), max(dataX)
            dy = max(dataY) - min(dataY)
            if (cf.amp[i] <= 0) or (cf.aerr[i]/cf.amp[i] > self.cfg.rerr):
                ampcolor = color2
            elif cf.aerr[i]/cf.amp[i] > 0.5*self.cfg.rerr:
                ampcolor = color1
            if np.fabs(cf.off[i]) > self.cfg.doff:
                offcolor = color2
            elif np.fabs(cf.off[i]) > 0.5*self.cfg.doff:
                offcolor = color1
            if np.fabs(cf.hpbw[i]/self.cfg.beam-1) > self.cfg.dwidth:
                hpbwcolor = color2
            elif np.fabs(cf.hpbw[i]/self.cfg.beam-1) > 0.5*self.cfg.dwidth:
                hpbwcolor = color1

            if not sub_flag[i]:
                ymin = min(dataY)
                ymax = max(dataY) +0.1*dy

            else:
                ymin = max(min(dataY), -0.1*dy)
                ymax = min(max(dataY)+0.1*dy, 1.5*cf.amp[i])

            if cf.amp[i] > 0:
                base_flt = (dataX >= cf.off[i]-window[0]*cf.hpbw[i])* \
                        (dataX <= cf.off[i]+window[1]*cf.hpbw[i])
                baseX = dataX[base_flt]
                baseY = np.ones(len(baseX))*0.1*(ymax-ymin)
                baseY[0] = ymin
                baseY[-1] = ymin

            pgsch(2)
            pgslw(2)
            pgsci(1)
            pgenv(xmin, xmax, ymin, ymax, 0, 0)
            pglab(cf.scandir[i]+'-OFF ["]', 'COUNTS', \
                    'Scan: %04d  %s  %0.1fMHz' \
                    %(cf.scannum, cf.src, self.cfg.freq))
            pgline(dataX, dataY)

            if cf.amp[i] > 0:
                pgsci(color0)
                pgslw(3)
                pgline(dataX, fitY)
                pgsci(color1)
                pgslw(2)
                pgline(baseX, baseY)

                pgsci(color0)
                pgptxt(0.95*xmax+0.05*xmin, 0.9*ymax+0.1*ymin, 0, 1, \
                        'Az=%0.1f' %cf.azi[i])

                pgptxt(0.95*xmax+0.05*xmin, 0.815*ymax+0.185*ymin, 0, 1, \
                        'El=%0.1f' %cf.elv[i])

                pgptxt(0.95*xmax+0.05*xmin, 0.73*ymax+0.27*ymin, 0, 1, \
                        'Sun=%0.1f' %cf.sunang[i])

                pgsci(ampcolor)
                pgptxt(0.05*xmax+0.95*xmin, 0.9*ymax+0.1*ymin, 0, 0, \
                        'Amp=%0.4f' %cf.amp[i])

                pgsci(offcolor)
                pgptxt(0.05*xmax+0.95*xmin, 0.815*ymax+0.185*ymin, 0, 0, \
                        'Off=%0.2f' %cf.off[i])

                pgsci(hpbwcolor)
                pgptxt(0.05*xmax+0.95*xmin, 0.73*ymax+0.27*ymin, 0, 0, \
                        'HPBW=%0.2f' %cf.hpbw[i])

            if sub_flag[i] == False:
                pgsci(color2)
                pgrect(xmax, xmax-(xmax-xmin)/25.0, ymin, ymin+(ymax-ymin)/25.0)

        xmin = -0.1*(max(cf.mjd) - min(cf.mjd))*1440
        xmax = 1.1*(max(cf.mjd) - min(cf.mjd))*1440
        if True in sub_flag:
            min_y = min(cf.amp[sub_flag] - cf.aerr[sub_flag])
            max_y = max(cf.amp[sub_flag] + cf.aerr[sub_flag])
            ymin = min_y - 0.1*(max_y - min_y)
            ymax = max_y + 0.1*(max_y - min_y)
            mean_sub_flag = np.mean(cf.amp[sub_flag])
        else:
            ymin = -1
            ymax = 1
            mean_sub_flag = 0

        cl = Calibration.ListFlux(self.cfg, '0.point', '0.point')
        cl.load_listflux()
        qlt = cl.load_qlt_scan(scannum)

        pgsci(1)
        pgsch(2)
        pgenv(xmin, xmax, ymin, ymax, 0, 0)
        if not qlt:
            pgsci(color2)
            pgrect(xmax, xmax-(xmax-xmin)/25.0, ymin, ymin+(ymax-ymin)/25.0)
            pgsci(1)
        pgline(np.array([xmin, xmax]),
                np.array([mean_sub_flag, mean_sub_flag]))
        pglab('DURATION', 'COUNTS', 'OBSTIME: %s' %al.mjd2utc(np.mean(cf.mjd)))
        pgsch(3)
        pgsci(color0)
        for i in range(subs):
            if sub_flag[i]:
                pgpt(np.array([(cf.mjd[i]-min(cf.mjd))*1440]),
                        np.array([cf.amp[i]]), 17)
                pgerrb(2, np.array([(cf.mjd[i]-min(cf.mjd))*1440]),
                        np.array([cf.amp[i]]), np.array([cf.aerr[i]]), 2)
                pgerrb(4, np.array([(cf.mjd[i]-min(cf.mjd))*1440]),
                        np.array([cf.amp[i]]), np.array([cf.aerr[i]]), 2)
            else:
                pgsci(color2)
                pgpt(np.array([(cf.mjd[i]-min(cf.mjd))*1440]),
                        np.array([cf.amp[i]]), 10)
                pgsci(color0)
        pgsci(1)
        pgsch(1)
        pgslw(1)


    def do_plot_ps(self, arg):
        """
        Plot Gaussian profile in ps files.
        """
        scans = os.listdir(self.cfg.fit_path)
        scans = [e[:-4] for e in scans if e[-4:] == '.fit']
        scans.sort()
        if arg == '-single':
            self.do_openX('cps fit.ps 16 11')
            for scannum in scans: self.do_chk_fit(scannum)
            self.do_closeX(None)
        else:
            for scannum in scans:
                #os.chdir(self.cfg.fit_path)
                #self.do_openX('cps %04d.ps 16 11' %int(scannum))
                self.do_openX('cps %04d.ps 16 11' %int(scannum))
                #os.chdir("..")
                self.do_chk_fit(scannum)
                os.rename('%04d.ps' %int(scannum), '%s/%04d.ps' \
                        %(self.cfg.fit_path, int(scannum)))
                self.do_closeX(None)


    def do_light_curve(self, arg):

        import Gnuplot

        """
        example: light_curve 5.flux
        Produce a ps file: fluxes.ps in which all light curve are plotted.
        """

#       first find the file names and sources names
        idpth = arg
#        idpth = core.locate_proc("flux")[1]
#        line = arg.split(" -")
#        if core.auto_find("idpth", line)[0]:
#            idpth = core.auto_find("idpth", line)[1][0].split()[1]
#            ver = core.locate_proc(idpth)[0]
#
        with open("mod.index") as modfile:
            mod = modfile.readlines()
        m0 = float(mod[0].split()[-1])
        mod = mod[4:]
        n = len(mod)
        mod = [e[:-1].split() for e in mod]
        srcs = [mod[i][0] for i in range(n)]
        S = [float(mod[i][2]) for i in range(n)]
        sigma = [float(mod[i][3]) for i in range(n)]
        m = [float(mod[i][4]) for i in range(n)]
        Y = [float(mod[i][5]) for i in range(n)]
        chi = [float(mod[i][6]) for i in range(n)]
        prob = [float(mod[i][8]) for i in range(n)]
        CI = [float(mod[i][9]) for i in range(n)]

        infname = [idpth + "/FLUX." + e for e in srcs]
#        infname.sort()
        g = Gnuplot.Gnuplot(debug=1)
        g("set bar 0")
        g("set border lw 0")
        g("set xtics scale 1.2,0.6 off 0,0.2 1")
        g("set ytics nomirror scale 1.6,0.8 off 0.2,0")
        g("set y2tics scale 1.6,0.8")
        g("set format y '%0.2f'")
        g("set format y2 '%0.2f'")
        g("set terminal postscript size 6.5,8.8 'Times-Roman' 10 enhanced")
        g("set output 'fluxes.ps'")
        g("set multiplot")

        j = -1
        for i in range(n):
            try:
                inf  = open(infname[i], "r")
            except IOError:
                continue
            data = inf.readlines()
            inf.close()
            data = [e[:-1] for e in data if len(e)>1]
            srcname = data[0].split()[1]

            data = [e for e in data if e[0]!= "#"]
            if len(data) > 1:
                j+=1
                y1 = np.array([float(e.split()[3]) for e in data])
                err_y1 = np.array([float(e.split()[4]) for e in data])
                min_y1 = np.floor(min(y1-err_y1)*100)/100.
                max_y1 = np.ceil(max(y1+err_y1)*100)/100.
                dy1 = max_y1 - min_y1

                x1 = [float(e.split()[0]) for e in data]
                min_x1 = np.floor(min(x1))
                max_x1 = np.ceil(max(x1))
                if min_x1+0.5<min(x1):min_x1+=0.5
                if max_x1-0.5>max(x1):max_x1-=0.5
                if min(x1)-min_x1<0.2:min_x1-=0.5
                if max_x1-max(x1)<0.2:max_x1+=0.5

                if (max_y1-min_y1)>0.1: g("set mytics 5")
                else: g("set mytics 10")
                if (max_y1-min_y1)/S[i]>0.1: g("set my2tics 5")
                else: g("set my2tics 10")

                g("set mxtics 10")

                title = "{/Helvetica-Bold=12 %s} {/Times-Roman=15: %-15s <I> =  %8.3f  +-  %-7.3f Jy}" \
                %(srcname,"Intensity", S[i], sigma[i])
#               g.title(s=title,offset=(0, -0.4),font=None)
                g.title(s=title,font=None)
                g("set xrange [%0.1f:%0.1f]" %(min_x1,max_x1))
#               g("set yrange [%0.2f:%0.2f]" %(min_y1-0.1*dy1, max_y1+0.1*dy1))
                g("set yrange [%0.2f:%0.2f]" %(min_y1, max_y1))
                g("set y2range[%0.2f:%0.2f]" %(min_y1/S[i], max_y1/S[i]))
#               g.ylabel("I [Jy]",offset=(1, 0))
                g.ylabel("I [Jy]")
                dx = (max_x1-min_x1)*0.1
                dy = (max_y1-min_y1)*0.12
                g("set size 1.0,0.34")
                g("set origin 0,%f" %((2-j%3)*0.33))
                if j!=0 and j%3==0:
                    g("set nomultiplot")
                    g("set multiplot")

                g("set label 'm = %0.1f[%%]' at graph 0.05, graph 0.92" %(m[i]))
                g("set label 'chi-square = %0.2f' at graph 0.05, graph 0.86" %(chi[i]))
                g("set label 'p-value = %0.3E' at graph 0.05, graph 0.80" %(prob[i]))
                g("set label 'CI = %0.2f' at graph 0.05, graph 0.74" %(CI[i]))
                f = Gnuplot.File(infname[i])
                f.set_option(title = None)
                f.set_option(with_="errorbars lw 3 lt rgb 'blue' pt 1 ps 0.8")
                f.set_option(using=(1,4,5))
                g.plot(f)
                g.plot(Gnuplot.Data([min_x1, max_x1], [1, 1], with_ = "l lt 2 axes x1y2"))
                g.plot(Gnuplot.Data([min_x1, max_x1], [1+3*m[i]/100., 1+3*m[i]/100.], with_ = "l lt 2 lc rgb 'blue' axes x1y2"))
                g.plot(Gnuplot.Data([min_x1, max_x1], [1-3*m[i]/100., 1-3*m[i]/100.], with_ = "l lt 2 lc rgb 'blue' axes x1y2"))
                g.plot(Gnuplot.Data([min_x1, max_x1], [1+3*m0/100., 1+3*m0/100.], with_ = "l lt 2 lc rgb 'red' axes x1y2"))
                g.plot(Gnuplot.Data([min_x1, max_x1], [1-3*m0/100., 1-3*m0/100.], with_ = "l lt 2 lc rgb 'red' axes x1y2"))
                g("set nolabel")


    def do_fit(self, line):
        """
fit -scan [0,123] -chan 0.5*(L+R) -nchan 4 -nphase 4 -denoise 0
    -cut [0.1,0.1] -window [1.35,1.35] -npeak 1 -bo 1 -smodel [0,0,0,0,0]

    Gauss fitting to the FITS data.

    Description:

      Since the observations are performed in "cross-scan" mode, for each
      sub-scan, the received power is the convolution of the source brightness
      distribution and the antenna beam pattern. For an unresolved
      point-source, this convolution can be well approximated by a Gaussian
      profile.

    Arguments:

      -scan    -- scan number, or [S,E]. 'S' and 'E' indicate 'start number' and
                  'end number' respectively. None -> all.
      -chan    -- determin which channel will be reduced. Can be any kind of
                  function about L and R, such as 0.5*(L+R), which means average
                  of left and right channel. None -> 0.5*(L+R).
      -nchan   -- number of channel. None -> 4.
      -nphase  -- number of modulation phase. None -> 4.
      -denoise -- options to denoise data with interference. It's value is the
                  hight of spike in multiples of the RMS. None -> 0, which
                  means turn of denoise.
      -cut     -- restricts the range of the Gaussian fit in scanning direction
                  on the left and right. E.g. cut=[0.1,0.1] will reduce the
                  range by 10% on both sides. None -> [0.1,0.1].
      -window  -- defines the window size in multiples of the beam width for the
                  second iteration of the Gaussian fit. None -> [1.35, 1.35].
      -npeak   -- number of Gaussian peaks. None -> 1.
      -bo      -- order of baseline. None -> 1.
      -smodel  -- source model. None -> [0,0,0,0,0].

    Outputs: Generate xxxx.fit and xxxx.dat files.

        """

        self.cfg.load_path()
        self.cfg.load_telescope()

        args = {}
        scan = [0, 9999]
        if not os.path.isdir(self.cfg.fit_path):
            os.mkdir(self.cfg.fit_path)

        if 'scan' in line:
            scan = eval(self.parm_exact(line, 'scan'))
        if 'chan' in line:
            args['chan'] = self.parm_exact(line, 'chan')
        if 'nchan' in line:
            args['nchan'] = int(self.parm_exact(line, 'nchan'))
        if 'nphase' in line:
            args['nphase'] = int(self.parm_exact(line, 'nphase'))
        if 'denoise' in line:
            args['despike'] = float(self.parm_exact(line, 'denoise'))
        if 'cut' in line:
            args['cut'] = eval(self.parm_exact(line, 'cut'))
        if 'window' in line:
            args['window'] = eval(self.parm_exact(line, 'window'))
            self.window = args['window']
        if 'npeak' in line:
            args['npeak'] = int(self.parm_exact(line, 'npeak'))
        if 'bo' in line:
            args['bo'] = int(self.parm_exact(line, 'bo'))
        if 'smodel' in line:
            args['smodel'] = eval(self.parm_exact(line, 'smodel'))

        fname = os.listdir(self.cfg.fits_path)
        fname.sort()
        fname = [int(e[:4]) for e in fname if e.endswith('.FITS')]
        fname = [e for e in fname if (e>=scan[0] and e<=scan[1])]
        fname = ['%04d.FITS' %e for e in fname]

        rd = ObsData.RawData
        for fn in fname:
            rd(self.cfg, fn, **args).fit_data()

        print('Fitting complete!')


    def do_pointing(self, line):
        """
pointing   -link 0 -flag 1

    Correct the lost flux which is due to pointing offset.

    Description:

      As the antenna is not always pointing to the source position acurately,
      the measured source fluxes are lower than the real. For an estimation of
      Gaussian profile, the observed fluxe was corrected by computing the
      difference between the measured and real (model predicted) flux.

    Arguments:

      -link -- Determine whether to generate the LIST.fit file or not. Shoule be
               0 or 1. The LIST.fit files, which are used in pointing, locate in
               current directory. If this file already exist, to save time, you
               can set link to 0. Default is 0.
      -flag -- Determine whether to flag the bad subscans or not. Should be 0 or
               1. If 0, all subscans in one scan are considered good and will be
               taken into account. While if flag == 1, it judge the fitting
               conditions (by offset, HPBW, relative error, and source
               symmetry). The bad subscans are flagged and won't be taken into
               account in pointing correction.

    Outputs: Generate LIST.raw file and write data into Raw.* file.
             Write source names to config file.
        """

        tcal = 1.0
        cf = Calibration.FitDat(self.cfg, 0, 1)

        if 'link' in line:
            cf.lnk = int(self.parm_exact(line, 'link'))

        if 'flag' in line:
            cf.flag = int(self.parm_exact(line, 'flag'))

        if 'tcal' in line:
            cf.tcal = float(self.parm_exact(line, 'tcal'))
        else:
            cf.tcal = tcal

        if cf.flag:
            self.cfg.load_quality_rules()

        cf.corr_point()

        print('Pointing correction complete.')

        cl = Calibration.ListFlux(self.cfg, '0.point', '0.point')
        cl.load_listflux()
        cl.split_listflux()


    def do_split_point_tau(self, line):

        clf = Calibration.ListFlux(self.cfg, '0.point', '0.point')
        clf.load_listflux()
        clf.split_listflux()
        del clf

        clf = Calibration.ListFlux(self.cfg, '0.point', '1.tau')
        clf.load_listflux()
        clf.split_listflux()
        del clf


    def do_split_flux(self, line):

#        clf = Calibration.ListFlux(self.cfg, '5.flux', '5.flux')
        clf = Calibration.ListFlux(self.cfg, line, line)
        clf.load_listflux()
        clf.split_listflux()
        del clf


    def do_get_tau(self, line):

        """
get_tau -tgr 0.0 -tcal 1.0

    Calculate the best air opacity by fitting the lower envelop of
    tsys vs. airmass.

    Arguments:

      -tgr -- Ground temperature in unit of degree centigrade.
              Default is 0 deg.cent.
      -tcal -- Default 1.0

    Output: Write information of tau to config file.
        """

        tgr = 0.0
        tcal = 1.0
        if 'tgr' in line:
           tgr = float(self.parm_exact(line, 'tgr'))
        if 'tcal' in line:
           tcal = float(self.parm_exact(line, 'tcal'))

        pid = pgqid()
        if not pid:
            print('No PGPLOT device is selected!')
            return
        pgslct(pid)
        pgask(0)

        def plot_tau(clf):

            if pgqinf('device') == '':
                ptcolor = 3
            else:
                ptcolor = 4

#            func = np.array([clf.cfg.Tau0*tatm, clf.cfg.Trec])
#
#            x = np.array(clf.data[7])[clf.qlt]
#            x = 1.0/np.sin(np.deg2rad(x))
#            clf = Calibration.ListFlux(self.cfg, '0.point', '1.tau')
            func, x = clf.compute_tau(tgr, tcal)
            x = x[clf.qlt]
            y = np.array(clf.data[11])[clf.qlt]*tcal

            xmin = min(x) - 0.1*(max(x) - min(x))
            xmax = max(x) + 0.1*(max(x) - min(x))
            ymin = min(y) - 0.1*(max(y) - min(y))
            ymax = max(y) + 0.1*(max(y) - min(y))
            XYscale = (xmax - xmin)/(ymax - ymin)

            pgpage()
            pgsch(1.5)
            pgsls(1)
            pgslw(2)
            pgsci(1)

            pgsvp(6.0/80, 1-1.0/80, 8.0/80, 1-5.0/80)
            pgswin(xmin, xmax, ymin, ymax)
            pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
            pgmtxt('L', 2, 0.5, 0.5, 'Tsys')
            pgmtxt('B', 2.4, 0.5, 0.5, 'sin(elv.\\u-1\d)*T\datm\\u')
            pgmtxt('T', 1, 0.5, 0.5, 'T\dsys\\u=%0.2f + %0.4f*Airmass*T\datm\\u' \
                    %(func[1], func[0]))
            pgsch(2)
            pgsci(ptcolor)
            pgpt(x, y, 17)
            pgsch(1.5)
            pgsci(1)

            tsys0 = np.polyval(func, 1)
            pgptxt(0.95*xmin+0.05*xmax, 0.9*ymax+0.1*ymin, 0, 0, \
                    'T\damb\\u=%0.2f' %clf.cfg.Tgr)

            pgptxt(0.95*xmin+0.05*xmax, 0.815*ymax+0.185*ymin, 0, 0, \
                    'T\dsys\\u0=%0.2f' %tsys0)

            pgptxt(0.95*xmin+0.05*xmax, 0.73*ymax+0.27*ymin, 0, 0, \
                    '\gt\d0\\u=%0.4f' %clf.cfg.Tau0)

            pltx = np.linspace(xmin, xmax, 3)
            plty = np.polyval(func, pltx)
            pgsci(2)
            pgslw(5)
            pgsls(2)
            pgline(pltx, plty)

            return x, y, XYscale

        clf = Calibration.ListFlux(self.cfg, '0.point', '1.tau')
        clf.load_listflux()
        clf.compute_tau(tgr, tcal)
        x, y, XYscale = plot_tau(clf)
        clf.write_tau()
        while 1:
            if pgqinf('cursor') == 'NO': # none-interactive device
                break

            gesture = pgband(0)
            if gesture[2].lower() == 'a': # left click
                x, y, XYscale = plot_tau(clf)
                gx, gy = gesture[0], gesture[1]
                dist = (x-gx)**2 + XYscale**2*(y-gy)**2
                dist = dist**0.5
                idx = dist.argmin()
                pick_scannum = clf.data[2][clf.qlt][idx]
                idx_scan = list(clf.data[2]).index(pick_scannum)
                pgsci(1)
                pgpt(np.array([x[idx]]), np.array([y[idx]]), 23)

            elif gesture[2].lower() == 'x': # right click
                try: pick_scannum
                except: return

                old_id = pgqid()
                if not self.gs_pid:
                    self.gs_pid = pgopen('/xw')
                    self.all_pid.append(self.gs_pid)
                pgslct(self.gs_pid)
                self.do_chk_fit(str(pick_scannum))
                pgslct(old_id)

            elif gesture[2].lower() == 'd': # middle click, delet
                try: idx_scan
                except: return
                self.flag_all_listflux(idx_scan)
                clf.compute_tau(tgr, tcal)
                x, y, XYscale = plot_tau(clf)
                clf.write_tau()
                pgsci(1)
                pgpt(np.array([1.0/np.sin(np.deg2rad(clf.data[7][idx_scan]))]),
                        np.array([clf.data[11][idx_scan]]), 23)

            elif gesture[2].lower() == "q" or gesture[2] == "": break


    def do_tau_corr(self, line):
        """
tau_cor

    Perform air opacity correction for observing data.
    tau_corr 0.2
        """

        off = 0
        tele = 'NSRT'

        clf = Calibration.ListFlux(self.cfg, '0.point', '1.tau')
        if 'off' in line:
            off = float(self.parm_exact(line, 'off'))
        if 'tele' in line:
            tele = float(self.parm_exact(line, 'tele'))
        clf.corr_tau(tele=tele, off=off)
        print('tau correction done.')


    def do_get_gain_curve(self, line):
        """
get_gain_curve -inpth 1.tau -order 2 -weight 0 -limit 3

    Find the Gain-Elevation curve for observed calibrators by a polynomial
    fitting to their normalized flux density versus elevation.

    Description:

      Fit a least squares polynomial "p(x) = p[0] * x**deg + ... + p[deg]"
      Tolerance the given set of data-points (Elv[i], Amp[i], Err[i]).

    Arguments:

      -inpth  -- Input data path. The location of data files which will be used
                 for the original data input. Default value is 1.tau, that means
                 we do this based on the tau corrected data.
      -order  -- Order of degree for the fitting polynomial. Default is 2.
      -limit  -- Tolerance of auto judgement of bad data-points in determining
                 the Gain-Elevation curve. The function calculate the parameters
                 of the curve (value and its error, which contributes the sigma
                 of the curve). Data-points that lie beyond -/+limit*sigma will
                 not be taken into account.
      -weight -- Strictly positive number.
                 The weights are used in computing the weighted least-square
                 polynomial fit. If the errors in the Amp values given by the
                 array Err, then weight should be (1/Err)^weight. Default is 0.
      -peak0  -- The expected peak position in elevation for the gain curve.
                 Default value: 45.

    Outputs: Write Gain-Elevation curve parameters into a text file named
             gain.curve.

        """

        proc0 = '1.tau'
        weight = 0
        order = 2
        limit = 3
        peak0 = 45.0

        if 'inpth' in line:
            proc0 = self.parm_exact(line, 'inpth')

        if 'weight' in line:
            weight = int(self.parm_exact(line, 'weight'))

        if 'order' in line:
            order = int(self.parm_exact(line, 'order'))

        if 'limit' in line:
            limit = float(self.parm_exact(line, 'limit'))

        if 'peak0' in line:
            limit = float(self.parm_exact(line, 'peak0'))

        proc1 = '%d.gain' %(int(proc0[0])+1)
        pid = pgqid()
        if not pid:
            print('No PGPLOT device is selected!')
            return
        pgslct(pid)
        pgask(0)

        def plot_gain(clf, gain_cal, gain_peaks, sigma, xmin, xmax, ymin, ymax):

            XYscale = (xmax - xmin)/(ymax - ymin)

            pgpage()
            pgsch(1.5)
            pgsls(1)
            pgslw(1)
            pgsci(1)

            pgsvp(6.0/80, 1-1.0/80, 8.0/80, 1-5.0/80)
            pgswin(xmin, xmax, ymin, ymax)
            pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
            pgmtxt('L', 2, 0.5, 0.5, 'Norm. Flux')
            pgmtxt('B', 2.4, 0.5, 0.5, 'Elevation')
            pgmtxt('T', 1, 0.5, 0.5, 'Gain-Elevation Curve')

            ci = 0
            xx = np.array([])
            yy = np.array([])
            scans = np.array([])
            for gcal in gain_cal:
                ci += 1
                pgsci(ci+1)
                flt = (clf.src == gcal)*clf.qlt
                x = clf.data[7][flt]
                y = clf.data[4][flt]/gain_peaks[ci-1]
                yerr = clf.data[5][flt]/gain_peaks[ci-1]
                pgpt(x, y, 17)
                pgerrb(2, x, y, yerr, 1)
                pgerrb(4, x, y, yerr, 1)

                py0 = 1.025-0.085*ci
                py1 = 1 - py0
                pgsch(1)
                pgptxt(0.98*xmin+0.02*xmax, py0*ymax+py1*ymin, 0, 0, \
                    '%s' %gcal)
                pgsch(1.5)

                xx = np.append(xx, x)
                yy = np.append(yy, y)
                scans = np.append(scans, clf.data[2][flt])

            pltx = np.linspace(xmin, xmax, 100)
            plty = np.polyval(clf.gain, pltx)
            pgsci(1)
            pgsls(2)
            pgline(np.array([xmin, xmax]), np.ones(2))
            if limit:
                pgline(pltx, plty-limit*sigma)
                pgline(pltx, plty+limit*sigma)
            pgsls(1)
            pgslw(2)
            pgline(pltx, plty)

            scans = [int(e) for e in scans]

            return scans, xx, yy, XYscale

        # compute and then plot
        clf = Calibration.ListFlux(self.cfg, proc0, proc1)
        gain_cal, gain_peaks, sigma, xmin, xmax, ymin, ymax \
                = clf.compute_gain(weight, order, limit, peak0)
        scans, x, y, XYscale = \
        plot_gain(clf, gain_cal, gain_peaks, sigma, xmin, xmax, ymin, ymax)
        clf.write_gain()
        while 1:
            if pgqinf('cursor') == 'NO': # none-interactive device
                break

            gesture = pgband(0)
            if gesture[2].lower() == 'a': # left click
                # plot without computing
                scans, x, y, XYscale = \
                plot_gain(clf, gain_cal, gain_peaks,
                        sigma, xmin, xmax, ymin, ymax)
                gx, gy = gesture[0], gesture[1]
                dist = (x-gx)**2 + XYscale**2*(y-gy)**2
                dist = dist**0.5
                idx = dist.argmin()
                pick_scannum = scans[idx]
                idx_scan = list(clf.data[2]).index(pick_scannum)
                pgsci(1)
                pgpt(np.array([x[idx]]), np.array([y[idx]]), 23)

            elif gesture[2].lower() == 'x': # right click
                try: pick_scannum
                except: return

                old_id = pgqid()
                if not self.gs_pid:
                    self.gs_pid = pgopen('/xw')
                    self.all_pid.append(self.gs_pid)
                pgslct(self.gs_pid)
                self.do_chk_fit(str(pick_scannum))
                pgslct(old_id)

            elif gesture[2].lower() == 'd': # middle click, delet
                try: idx_scan
                except: return
                self.flag_all_listflux(idx_scan)
                # compute and plot
                gain_cal, gain_peaks, sigma, xmin, xmax, ymin, ymax \
                        = clf.compute_gain(weight, order, limit, peak0)
                scans, x, y, XYscale = \
                plot_gain(clf, gain_cal, gain_peaks,
                        sigma, xmin, xmax, ymin, ymax)
                clf.write_gain()

                pgsci(1)
                pgpt(np.array([1.0/np.sin(np.deg2rad(clf.data[7][idx_scan]))]),
                        np.array([clf.data[11][idx_scan]]), 23)

            elif gesture[2].lower() == "q" or gesture[2] == "": break


    def do_gain_corr(self, line):
        """
gain_corr -inpth 1.tau -err 0

    Gain-elevation correction. Read gain-elevation curve parameters from
    gain.curve, caculate the antenna gain at source elevation.

    Arguments:
      -inpth  -- Input data path. The location of data files which will be used
                 for the original data input. Default value is 1.tau, that means
                 we do this based on the tau corrected data.
      -err    -- err == 0 neglect the error of gain-elevation curve, just
                 consider the parameter value in gain_corr.
                 err == 1 take both the gain-elevation curve value and its error
                 into account.
        """

        proc0 = '1.tau'
        err = 0

        if 'inpth' in line:
            proc0 = self.parm_exact(line, 'inpth')

        if 'err' in line:
            err = float(self.parm_exact(line, 'err'))

        proc1 = '%d.gain' %(int(proc0[0])+1)
        clf = Calibration.ListFlux(self.cfg, proc0, proc1)
        clf.corr_gain()
        print('gain correction done.')


    def do_get_time_curve(self, line):
        """
get_time_curve -inpth 2.gain -bs 0.08

    Find the Gain-Time curve of observed calibrators.

    Description:

      Given the set of data-points (MJD[i], AMP[i], ERR[i]) and knots, evaluate
      the value of the peicewise cubic B-spline.

    Arguments:

      -inpth  -- Indata path. The location of data files which will be used for
                 the original data input. Default value is 2.gain, that means we
                 do this based on the gain corrected data.
      -bs     -- Binsize for data binning, in unit of day.
                 It works as smooth factor. The user can use bin to control the
                 tradeoff between closeness and smoothness of fit.  Larger bin
                 means more smoothing while smaller values of bin indicate less
                 smoothing. Default is 0.08, which is ~ 2 hrs.

    Outputs: Write Gain-Time curve parameters into a text file named time.curve.

        """

        proc0 = '2.gain'
        binsize = 0.08

        if 'inpth' in line:
            proc0 = self.parm_exact(line, 'inpth')

        if 'bs' in line:
            binsize = float(self.parm_exact(line, 'bs'))

        proc1 = '%d.time' %(int(proc0[0])+1)
        pid = pgqid()
        if not pid:
            print('No PGPLOT device is selected!')
            return
        pgslct(pid)
        pgask(0)

        def plot_time(clf, time_cal, xmin, xmax, ymin, ymax):

            XYscale = (xmax - xmin)/(ymax - ymin)

            pgpage()
            pgsch(1.5)
            pgsls(1)
            pgslw(1)
            pgsci(1)

            x0 = clf.data[1][0]
            xmin -= x0
            xmax -= x0

            pgsvp(6.0/80, 1-1.0/80, 8.0/80, 1-5.0/80)
            pgswin(xmin, xmax, ymin, ymax)
            pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
            pgmtxt('L', 2, 0.5, 0.5, 'Norm. Flux')
            pgmtxt('B', 2.4, 0.5, 0.5, 'MJD-%f' %x0)
            pgmtxt('T', 1, 0.5, 0.5, 'Gain-Time Curve')

            ci = 0
            xx = np.array([])
            yy = np.array([])
            scans = np.array([])
            for tcal in time_cal:
                ci += 1
                pgsci(ci+1)
                flt = (clf.src == tcal)*clf.qlt
                x = clf.data[1][flt]
                ymean = np.mean(clf.data[4][flt])
                y = clf.data[4][flt]/ymean
                yerr = clf.data[5][flt]/ymean
                pgpt(x-x0, y, 17)
                pgerrb(2, x-x0, y, yerr, 1)
                pgerrb(4, x-x0, y, yerr, 1)

                py0 = 1.025-0.085*ci
                py1 = 1 - py0
                pgsch(1)
                pgptxt(0.98*xmin+0.02*xmax, py0*ymax+py1*ymin, 0, 0, '%s' %tcal)
                pgsch(1.5)

                xx = np.append(xx, x)
                yy = np.append(yy, y)
                scans = np.append(scans, clf.data[2][flt])

            pgsci(1)
            pgsls(2)
            pgline(np.array([xmin, xmax]), np.ones(2))
            pgsls(1)
            pgslw(2)
            #pltx = np.linspace(clf.data[1][0], clf.data[1][-1], \
            #        3*len(clf.data[1]))
            #plty = interpolate.splev(pltx, clf.tck)
            #pgpt(pltx, plty, 17)
            pgline(clf.data[1]-x0, clf.tcurve)
            #pgpt(clf.data[1]-x0, clf.tcurve, 23)

            scans = [int(e) for e in scans]

            return scans, xx, yy, XYscale

        # compute and then plot
        clf = Calibration.ListFlux(self.cfg, proc0, proc1)
        time_cal, xmin, xmax, ymin, ymax  = clf.compute_time(binsize)
        scans, x, y, XYscale = \
        plot_time(clf, time_cal, xmin, xmax, ymin, ymax)
        clf.write_time()
        while 1:
            if pgqinf('cursor') == 'NO': # none-interactive device
                break

            gesture = pgband(0)
            if gesture[2].lower() == 'a': # left click
                # plot without computing
                scans, x, y, XYscale = \
                plot_time(clf, time_cal, xmin, xmax, ymin, ymax)
                x -= clf.data[1][0]
                gx, gy = gesture[0], gesture[1]
                dist = (x-gx)**2 + XYscale**2*(y-gy)**2
                dist = dist**0.5
                idx = dist.argmin()
                pick_scannum = scans[idx]
                idx_scan = list(clf.data[2]).index(pick_scannum)
                pgsci(1)
                pgpt(np.array([x[idx]]), np.array([y[idx]]), 23)

            elif gesture[2].lower() == 'x': # right click
                try: pick_scannum
                except: return

                old_id = pgqid()
                if not self.gs_pid:
                    self.gs_pid = pgopen('/xw')
                    self.all_pid.append(self.gs_pid)
                pgslct(self.gs_pid)
                self.do_chk_fit(str(pick_scannum))
                pgslct(old_id)

            elif gesture[2].lower() == 'd': # middle click, delet
                try: idx_scan
                except: return
                self.flag_all_listflux(idx_scan)
                # compute and plot
                time_cal, xmin, xmax, ymin, ymax = clf.compute_time(binsize)
                scans, x, y, XYscale = plot_time(clf, time_cal,
                        xmin, xmax, ymin, ymax)
                clf.write_time()

                src = clf.src[idx_scan]
                flt = (clf.src == src)*clf.qlt
                tmpx, tmpy = clf.data[1][idx_scan], \
                        clf.data[4][idx_scan]/np.mean(clf.data[4][flt])
                tmpx -= clf.data[1][0]

                pgpt(np.array([tmpx]), np.array([tmpy]), 23)

            elif gesture[2].lower() == "q" or gesture[2] == "": break



    def do_time_corr(self, line):
        """
time_corr -inpth 2.gain -err 1.0

    gain-time correction. Read gain-time curve parameters from time.curve,
    caculate the antenna gain at source observing time.

    Arguments:
      -inpth  -- Indata path. The location of data files which will be used for
                 the original data input. Default value is 2.gain, that means we
                 do this based on the gain corrected data.
      -outpth -- Output data path. The location where time corrected result will
                 be saved. Default is 3.time.
      -err    -- err == 0 neglect the error of gain-time curve, just
                 consider the parameter value for time_corr.
                 err == 1 take both the gain-time curve value and its error
                 into account.
        """

        proc0 = '2.gain'
        err = 1.0

        if 'inpth' in line:
            proc0 = self.parm_exact(line, 'inpth')

        if 'err' in line:
            err = float(self.parm_exact(line, 'err'))

        proc1 = '%d.time' %(int(proc0[0])+1)
        clf = Calibration.ListFlux(self.cfg, proc0, proc1)
        clf.corr_time(err)
        print('time correction done.')


    def do_get_flux_factor(self, line):
        """
get_flux_factor -inpth 3.time

    A link between the measured temperature and absolute flux density.
    It's determined by the flux density of primary calibrators 3C286, 3C48 and
    NGC7027 (flux density 7.53, 5.53 and 5.47 at 4800MHz, respectively).
    Return the factor and also its error.

    Arguments:
      -inpth -- Indata path. The location of data files which will be used for
                the original data input. Default value is 3.time, that means we
                do this based on the time corrected data.

    Output: Write factor to cal.factor.
        """

        proc0 = '3.time'

        if 'inpth' in line:
            proc0 = self.parm_exact(line, 'inpth')

        proc1 = '%d.flux' %(int(proc0[0])+1)
        clf = Calibration.ListFlux(self.cfg, proc0, proc1)
        clf.compute_dpfu()


    def do_abs_flux_corr(self, line):
        """
abs_flux_cor -inpth 3.time -err 1

    Convert the measured temperature to absolute flux density.

    Arguments:
      -inpth -- Indata path. The location of data files which will be used for
                the original data input. Default value is 0.raw, that means we
                do this based on the pointing corrected data.
      -err   -- Determine how much of the factor_error will be propagated,
                ranging from 0 to 1.
        """

        proc0 = '3.time'
        err = 1.0

        if 'inpth' in line:
            proc0 = self.parm_exact(line, 'inpth')

        if 'err' in line:
            err = float(self.parm_exact(line, 'err'))

        proc1 = '%d.flux' %(int(proc0[0])+1)
        clf = Calibration.ListFlux(self.cfg, proc0, proc1)
        clf.corr_dpfu(err)
        print('absolute flux density correction done.')


    def do_summary(self, line):
        """
summary -inpth 4.flux -m0 0.5

    Creat a summary text file: mod.index, in which some essential parameters
    of the observation are given.
        """

        proc0 = '0.point'
        proc1 = '1.tau'
        m0 = 0

        if 'inpth' in line:
            proc0 = self.parm_exact(line, 'inpth')

        if 'm0' in line:
            m0 = float(self.parm_exact(line, 'm0'))/100.0

        clf = Calibration.ListFlux(self.cfg, proc0, proc0)
        clf.get_m0(m0)



    def do_view(self, arg):
        """
view -show 0 -line 0 -sets 0.point mjd amp err 3c286,0893+710

    Plot data to a pgplot device.

    Arguments:
      -show -- Determine whether to show the flagged data-points or not.
      -line -- Determine whether to draw a line passing through all data
               points or not.
      -sets -- <1> Data location
               <2>, <3> x and y axis.
               <4> Error of <3> or source name
               <5> Only exist if <4> is error. The source name can be one
                   single source name or several source name divided by comma,
                   or "calibrator".
        """
        stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)

        if not pgqid():
            print('No PGPLOT device is selected!')
            return
        pid = pgqid()
        if not pid:
            print('No PGPLOT device is selected!')
            return
        pgslct(pid)
        pgask(0)

        show, line = 0, 0
        sets = ["1.tau", "mjd", "amp", "err", "calibrator"]

        if 'show' in arg:
            show = int(self.parm_exact(arg, 'show'))

        if 'line' in arg:
            line = int(self.parm_exact(arg, 'line'))

        if '-sets' in arg:
            tmp = arg.split()
            idx0 = tmp.index('-sets') + 1
            tmp = tmp[idx0:]
            tmp0 = [e for e in tmp if e.startswith('-')]
            if tmp0:
                idx1 = tmp.index(tmp0[0])
                tmp = tmp[:idx1]
            sets = tmp

        if ',' in sets[-1]: tmp = sets[-1].split(',')
        else: tmp = [sets[-1]]
        srcs = []
        for s in tmp:
            if 'calib' in s.lower():
                ss = self.cfg.pri_cal + self.cfg.gain_cal + self.cfg.time_cal
                ss = list(set(ss))
                ss.sort()
                srcs += ss
            else:
                srcs += [s.upper()]

        dict_idx = {'mjd':1,
                    'scan': 2,
                    'amp': 4,
                    'err': 5,
                    'azi': 6,
                    'elv': 7,
                    'utc': 8,
                    'lst': 9,
                    'pa': 10,
                    'tsys': 11,
                    'cal': 12,
                    'sunang': 13}

        idxx = dict_idx[sets[1]]
        idxy = dict_idx[sets[2]]
        if len(sets) == 4: idxz = None
        else: idxz = dict_idx[sets[3]]

        def view(clf, show, line, srcs, idxx, idxy, idxz,
                xmin, xmax, ymin, ymax):

            XYscale = (xmax - xmin)/(ymax - ymin)

            pgpage()
            pgsch(1.5)
            pgsls(1)
            pgslw(1)
            pgsci(1)

            pgsvp(6.0/80, 1-1.0/80, 8.0/80, 1-5.0/80)
            if idxx == 1: x0 = clf.data[idxx][0]
            else: x0 = 0
            pgswin(xmin-x0, xmax-x0, ymin, ymax)
            pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
            if 'flux' in clf.proc0:
                if len(srcs) == 1:
                    pgmtxt('L', 2, 0.5, 0.5, 'Flux Density [Jy]')
                else:
                    pgmtxt('L', 2, 0.5, 0.5, 'Norm. Flux Density')
            else:
                if len(srcs) == 1:
                    pgmtxt('L', 2, 0.5, 0.5, 'Counts')
                else:
                    pgmtxt('L', 2, 0.5, 0.5, 'Norm. Counts')
            if not x0:
                pgmtxt('B', 2.4, 0.5, 0.5, '%s' %sets[1].upper())
            else:
                pgmtxt('B', 2.4, 0.5, 0.5, '%s-%f' %(sets[1].upper(), x0))
            pgmtxt('T', 1, 0.5, 0.5, 'Light-curve for %s' %sets[-1])

            ci = 0
            xx = np.array([])
            yy = np.array([])
            scans = np.array([])
            for src in srcs:
                ci += 1
                pgsci(ci+1)
                if show: flt = clf.src == src
                else: flt = (clf.src==src)*clf.qlt
                if not True in flt:
                    ci -= 1
                    continue
                x = clf.data[idxx][flt]
                if len(srcs) == 1: ymean = 1.0
                else: ymean = np.mean(clf.data[idxy][flt])
                y = clf.data[idxy][flt]/ymean

                qlt = clf.qlt[flt]
                qlt_inv = np.array([not e for e in qlt], dtype=bool)

                if line: pgline(x-x0, y)
                if show: pgpt(x[qlt]-x0, y[qlt], 17)
                else: pgpt(x-x0, y, 17)
                if idxz:
                    yerr = clf.data[idxz][flt]/ymean
                    if show:
                        pgerrb(2, x[qlt]-x0, y[qlt], yerr[qlt], 1)
                        pgerrb(4, x[qlt]-x0, y[qlt], yerr[qlt], 1)
                        pgsls(2)
                        pgerrb(2, x[qlt_inv]-x0, y[qlt_inv], yerr[qlt_inv], 1)
                        pgerrb(4, x[qlt_inv]-x0, y[qlt_inv], yerr[qlt_inv], 1)
                        pgsls(1)
                    else:
                        pgerrb(2, x-x0, y, yerr, 1)
                        pgerrb(4, x-x0, y, yerr, 1)

                px = 0.98*xmin+0.02*xmax-x0
                if len(srcs) > 1:
                    py0 = 1.0-0.045*ci
                    py1 = 1 - py0
                    pgsch(1)
                    pgptxt(px, py0*ymax+py1*ymin, 0, 0, '%s' %src)
                else:
                    if not idxz: err = clf.data[5][flt]
                    else: err = yerr
                    meany = np.mean(y[qlt])
                    x2 = np.sum(((y[qlt]-meany)/err[qlt])**2)
                    xr2 = x2/(len(err[qlt])-1)
                    prob = stats.chisqprob(x2, len(err[qlt])-1)

                    py1 = ymax - (ymax-ymin)*0.05
                    py2 = ymax - (ymax-ymin)*0.08
                    py3 = ymax - (ymax-ymin)*0.11
                    py4 = ymax - (ymax-ymin)*0.14
                    pgsch(1)
                    pgsci(3)
                    pgtext(px, py1, '%-12s %0.3f +- %0.3f'
                            %('COUNTS    =', meany, np.std(y[qlt])))
                    pgtext(px, py2, '%-13s %0.2f' %(r'red.Chi.\u2\d   = ', xr2))
                    pgtext(px, py3, '%-13s %0.4E' %('Probability = ', prob))
                    pgtext(px, py4, '%-13s %0.2f %%'
                            %('Mod. Index =', 100*np.std(y[qlt]/np.mean(y[qlt]))))

                pgsch(1.5)
                pgsci(ci)
                xx = np.append(xx, x)
                yy = np.append(yy, y)
                scans = np.append(scans, clf.data[2][flt])

            pgsci(1)
            pgsls(2)
            pgline(np.array([xmin-x0, xmax-x0]), np.ones(2))
            pgsls(1)
            pgslw(2)

            scans = np.array([int(e) for e in scans])
            idxs = scans.argsort()
            scans = list(scans[idxs])
            xx = xx[idxs]
            yy = yy[idxs]

            return scans, xx, yy, XYscale


        def get_xybound(clf, srcs, idxx, idxy, idxz, show):

            xx = np.array([])
            yy = np.array([])
            yyerr = np.array([])
            for src in srcs:
                if show:
                    flt = clf.src == src
                else:
                    flt = (clf.src == src)*clf.qlt
                x = clf.data[idxx][flt]
                if len(srcs) == 1: ymean = 1.0
                else: ymean = np.mean(clf.data[idxy][flt])
                y = clf.data[idxy][flt]/ymean
                xx = np.append(xx, x)
                yy = np.append(yy, y)

                if idxz:
                    yerr = clf.data[idxz][flt]/ymean
                    yyerr = np.append(yyerr, yerr)

            dx = max(xx) - min(xx)
            if idxz:
                dy = max(yy+yyerr) - min(yy-yyerr)
                ymin = min(yy-yyerr) - 0.1*dy
                ymax = max(yy+yyerr) + 0.1*dy
            else:
                dy = max(yy) - min(yy)
                ymin = min(yy) - 0.1*dy
                ymax = max(yy) + 0.1*dy
            xmin = min(xx) - 0.1*dx
            xmax = max(xx) + 0.1*dx
            if xmin == xmax:
                xmin -= 0.1
                xmax += 0.1

            return xmin, xmax, ymin, ymax

        # compute and then plot
        clf = Calibration.ListFlux(self.cfg, sets[0], sets[0])
        clf.load_listflux()
        if sets[1] == 'mjd': x0 = clf.data[1][0]
        else: x0 = 0
        xmin, xmax, ymin, ymax  = get_xybound(clf, srcs, idxx, idxy, idxz, show)
        scans, x, y, XYscale = \
        view(clf, show, line, srcs, idxx, idxy, idxz, xmin, xmax, ymin, ymax)
        while 1:
            if pgqinf('cursor') == 'NO': # none-interactive device
                break

            gesture = pgband(0)
            if gesture[2].lower() == 'a': # left click
                # plot without computing
                scans, x, y, XYscale = view(clf, show, line, srcs, \
                        idxx, idxy, idxz, xmin, xmax, ymin, ymax)
                gx, gy = gesture[0], gesture[1]
                dist = (x-gx-x0)**2 + XYscale**2*(y-gy)**2
                dist = dist**0.5
                idx = dist.argmin()
                pick_scannum = scans[idx]
                pgsci(1)
                pgpt(np.array([x[idx]-x0]), np.array([y[idx]]), 23)

            elif gesture[2].lower() == 'x': # right click
                try: pick_scannum
                except: return

                if os.path.isfile(self.cfg.fit_path+'/%04d.dat' %pick_scannum):
                    old_id = pgqid()
                    if not self.gs_pid:
                        self.gs_pid = pgopen('/xw')
                        self.all_pid.append(self.gs_pid)
                    pgslct(self.gs_pid)
                    self.do_chk_fit(str(pick_scannum))
                    pgslct(old_id)
                else:
                    os.popen('gpicview %s/%04d.png &' %(self.cfg.fit_path,
                        pick_scannum)).read()

            elif gesture[2].lower() in ['f', 'b']:
                try: pick_scannum
                except: return
                scans, x, y, XYscale = view(clf, show, line, srcs, \
                        idxx, idxy, idxz, xmin, xmax, ymin, ymax)
                if gesture[2].lower() == 'f': idx += 1
                else: idx -= 1
                pick_scannum = scans[idx]
                pgsci(1)
                pgpt(np.array([x[idx]-x0]), np.array([y[idx]]), 23)

                if os.path.isfile(self.cfg.fit_path+'/%04d.dat' %pick_scannum):
                    old_id = pgqid()
                    if not self.gs_pid:
                        self.gs_pid = pgopen('/xw')
                        self.all_pid.append(self.gs_pid)
                    pgslct(self.gs_pid)
                    self.do_chk_fit(str(pick_scannum))
                    pgslct(old_id)
                else:
                    os.popen('gpicview %s/%04d.png &' %(self.cfg.fit_path,
                        pick_scannum)).read()

            elif gesture[2].lower() == 'd': # middle click, delet
                try: pick_scannum
                except: return
                idx_scan = list(clf.data[2]).index(pick_scannum)
                self.flag_all_listflux(idx_scan)
                clf.load_listflux()
                # compute and plot
#                xmin, xmax, ymin, ymax  = \
#                        get_xybound(clf, srcs, idxx, idxy, idxz, show)
                scans, x, y, XYscale = view(clf, show, line, srcs, \
                        idxx, idxy, idxz, xmin, xmax, ymin, ymax)

                src = clf.src[idx_scan]
                flt = (clf.src == src)*clf.qlt
                tmpx, tmpy = clf.data[1][idx_scan], \
                        clf.data[4][idx_scan]/np.mean(clf.data[4][flt])

                pgpt(np.array([tmpx-x0]), np.array([tmpy]), 23)

            elif gesture[2].lower() == "z":
                x00, y00, cursor = pgband(7)
                x11, y11, cursor = pgband(7)
                xmin, xmax, ymin, ymax = \
                        min(x00, x11), max(x00, x11), \
                        min(y00, y11), max(y00, y11)
                xmin += x0
                xmax += x0
                view(clf, show, line, srcs, idxx, idxy, idxz, \
                        xmin, xmax, ymin, ymax)

            elif gesture[2].lower() == "r":
                xmin, xmax, ymin, ymax  = \
                    get_xybound(clf, srcs, idxx, idxy, idxz, show)
                scans, x, y, XYscale = \
                    view(clf, show, line, srcs, idxx, idxy, idxz, \
                    xmin, xmax, ymin, ymax)

            elif gesture[2].lower() == "q" or gesture[2] == "": break


    def creat_pipe_script(self, clf):

        clf.load_listflux()
        ut0 = al.mjd2utc(clf.data[1][0])
        ut1 = al.mjd2utc(clf.data[1][-1])
        dhr = int((ut1 - ut0)*24)

        cals = self.cfg.pri_cal + self.cfg.gain_cal + self.cfg.time_cal
        cals = list(set(cals))
        cals.sort()
        flt = clf.src == cals[0]
        for cal in cals:
            flt += clf.src == cal
        count_cal = len(clf.src[flt])

        with open('pipeline', 'w') as fpipe:
            if dhr < 24 or count_cal < 20:
                print('# pipeline script procedure 1:', file=fpipe)
                print('# -----------------------------------', file=fpipe)
                print('  get_tau -tgr 0.0', file=fpipe)
                print('  tau_corr', file=fpipe)
                print('  get_gain_curve -inpth 1.tau -order 2 -weight 0 -limit 3', file=fpipe)
                print('  gain_corr -inpth 1.tau -err 0', file=fpipe)
                print('  get_flux_factor -inpth 2.gain', file=fpipe)
                print('  abs_flux_corr -inpth 2.gain -err 1.0', file=fpipe)
                print('  summary -inpth 3.flux -m0 0.5', file=fpipe)
                print('  view -line 1 -sets 2.gain mjd amp err calib', file=fpipe)

                print('# pipeline script procedure 2:', file=fpipe)
                print('# -----------------------------------', file=fpipe)
                print('# get_tau -tgr 0.0', file=fpipe)
                print('# tau_corr', file=fpipe)
                print('# get_gain_curve -inpth 1.tau -order 2 -weight 0 -limit 3', file=fpipe)
                print('# gain_corr -inpth 1.tau -err 0', file=fpipe)
                print('# get_time_curve -inpth 2.gain -bs 0.08', file=fpipe)
                print('# time_corr -inpth 2.gain -err 1.0', file=fpipe)
                print('# get_flux_factor -inpth 3.time', file=fpipe)
                print('# abs_flux_corr -inpth 3.time -err 1.0', file=fpipe)
                print('# summary -inpth 4.flux -m0 0.5', file=fpipe)
                print('# view -line 1 -sets 3.time mjd amp err calib', file=fpipe)

            else:
                print('# pipeline script procedure 1:', file=fpipe)
                print('# -----------------------------------', file=fpipe)
                print('  get_tau -tgr 0.0', file=fpipe)
                print('  tau_corr', file=fpipe)
                print('  get_gain_curve -inpth 1.tau -order 2 -weight 0 -limit 3', file=fpipe)
                print('  gain_corr -inpth 1.tau -err 0', file=fpipe)
                print('  get_time_curve -inpth 2.gain -bs 0.08', file=fpipe)
                print('  time_corr -inpth 2.gain -err 1.0', file=fpipe)
                print('  get_flux_factor -inpth 3.time', file=fpipe)
                print('  abs_flux_corr -inpth 3.time -err 1.0', file=fpipe)
                print('  summary -inpth 4.flux -m0 0.5', file=fpipe)
                print('  view -line 1 -sets 3.time mjd amp err calib', file=fpipe)

                print('# pipeline script procedure 2:', file=fpipe)
                print('# -----------------------------------', file=fpipe)
                print('# get_tau -tgr 0.0', file=fpipe)
                print('# tau_corr', file=fpipe)
                print('# get_time_curve -inpth 1.tau -bs 0.25 -err 0.0', file=fpipe)
                print('# time_corr -inpth 1.tau -err 0', file=fpipe)
                print('# get_gain_curve -inpth 2.time -order 2 -weight 0 -limit 3', file=fpipe)
                print('# gain_corr -inpth 2.time -err 0', file=fpipe)
                print('# get_time_curve -inpth 3.gain -bs 0.08', file=fpipe)
                print('# time_corr -inpth 3.gain -err 1.0', file=fpipe)
                print('# get_flux_factor -inpth 4.time', file=fpipe)
                print('# abs_flux_corr -inpth 4.time -err 0.8', file=fpipe)
                print('# summary -inpth 5.flux -m0 0.5', file=fpipe)
                print('# view -line 1 -sets 4.time mjd amp err calib', file=fpipe)


    def do_pipe(self, line):
        """
pipe

    Data reduction in pipeline mode. It excute the command listed in pipline
    file line by line (exclude lines which are comment out). The operator can
    manually modified the pipline file for a better quality of the result.
        """

        clf = Calibration.ListFlux(self.cfg, '0.point', '0.point')
        clf.load_listflux()
        if not os.path.isfile('pipeline'): self.creat_pipe_script(clf)
        command_list = [e[3:] for e in dir(self) if e[:3] == "do_"]
        with open('pipeline') as f:
            trs = f.readlines()
        trs = [e[:-1] for e in trs if e[0]!="#"]

        for tr in trs:
            command = tr.split()
            if len(command) > 1:
                arg = command[1:]
                arg = (' ').join(arg)
            else: arg = ''
            command = command[0]
            print('\n' + command + ' ' + arg)
            func = getattr(self, 'do_' + command)
            func(arg)


    def do_fdvar(self, line):

        cal = self.cfg.pri_cal + self.cfg.gain_cal + self.cfg.time_cal
        cal = list(set(cal))
        cal.sort()

        targets = self.cfg.target
        targets.sort()

#   command = '''cat mod.index | awk '$0!~"#" {print}' | sort -g -k9'''
        command = 'cat mod.index | ' + "awk '$0!~" + '"#" {print}' + "'" + " | sort -g -k9"
        trs = os.popen(command).read()[:-1]
        trs = trs.split("\n")

        tr_cal, tr_Nvar, tr_var = [], [], []

        for tr in trs:
            _tr = tr.split()
            if _tr[0] in cal:
                tr_cal.append(tr)
            elif _tr[0] in targets:
                if np.float(_tr[8]) < 0.001: tr_var.append(tr)
                else: tr_Nvar.append(tr)
        n_cal, n_Nvar, n_var = len(tr_cal), len(tr_Nvar), len(tr_var)
        tr_cal = ("\n").join(tr_cal)
        tr_Nvar = ("\n").join(tr_Nvar)
        tr_var = ("\n").join(tr_var)

        print("\n# ---------- %d Calibrators: ----------" %n_cal)
        print(tr_cal)
        print("\n# ---------- %d None-Variables: ----------" %n_Nvar)
        print(tr_Nvar)
        print("\n# ---------- %d Variables: ----------" %n_var)
        print(tr_var)


    def do_performance(self, line):
        """
performance -d 26.0 -f 4800.0 -tcal 1.0 -tg 10.0

    Calculate the system performance of an antenna.

    Arguments:

      -d -- diameter of the antenna in unit of meter. Default is 26.0 meters.
      -f -- Observing frequency in unit of MHz. Default is 4800.0 MHz.
      -tcal -- Strength of calibration signal generated by noise diode, in
               unit of K. Default is 1.0 K.
      -tg -- Ground temperature in unit of degree centigrade.
             Default is 10 deg.cent.

    Outputs: Generate following figures:
             hpbw.png, off_hpbw.png, pointing.png, offset.png tsys.png;
             and display the summary of the performance of the antenna.
        """

        D = 26.0
        F = 4800.0
        Tcal = 1.0
        Tgr = 10.0

        if '-d' in line:
            D = float(self.parm_exact(line, 'd'))

        if '-f' in line:
            F = float(self.parm_exact(line, 'f'))

        if '-tcal' in line:
            Tcal = float(self.parm_exact(line, 'tcal'))

        if '-tg' in line:
            Tgr = float(self.parm_exact(line, 'tg'))

        perform.sys_pfm(D, F, Tcal, Tgr)


    def do_how(self, line):

        how = []
        how.append('\nSteps of Continuum Data Reduction for Antenna System Measurement')
        how.append('*** 1. -> show_telescope')
        how.append('       Check the setup before data reduction.')
        how.append('*** 2. -> set_tele_frq 22400.0')
        how.append('***    -> set_tele_tcal 5.0')
        how.append('       Modify the default vaule of the parameters.')
        how.append('*** 3. -> fit -scan [100,200] -smodel [0,0,120,0,0]')
        how.append('       Gaussfit to the data, run "help fit" for details.')
        how.append('*** 4. -> pointing -link 1 -flag 1')
        how.append('       Perform pointing correction (with flagging) to the data.')
        how.append('       You may want to modify the flag options in config file,')
        how.append('       and then rerun "pointing -link 1" again.')
        how.append('*** 5. -> openX')
        how.append('       Open a pgplot window.')
        how.append('*** 6. -> pipe')
        how.append('       Do the rest steps by a pipeline.')
        how.append('       You may want to modify the script in pipline file,')
        how.append('       and then rerun "pipe" again.')
        how.append('*** 7. -> performance -f 22400.0 -tcal 5.0')
        how.append('       Calculate the performance of the antenna.')
        how.append('Run "help" to see more details of individual function.\n')

        for doc in how:
            print(doc)




if __name__ == "__main__":
    app = DRP()
    app.cmdloop('just a test')
