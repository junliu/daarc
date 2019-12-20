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
import numpy as np
from astropy.io import fits as pyfits
import AstroLib as al
import MathLib as ml
from Misc import readcol, writecol


class Config(object):

    def __init__(self, **kwargs):

        args = {'fits_path': './fits',
                'fit_path': './fit',
                'tele_name': 'NSRT',
                'obs_lon': 87.178108,
                'obs_lat': 43.4709389,
                'diam': 26.0,
                'freq': 4800.0,
                'tcal': 1.7}

        for key in list(args.keys()):
            if key in kwargs:
                setattr(self, key, kwargs[key])
            else:
                setattr(self, key, args[key])

        self.cfn = 'config'
        self.load_path()
        self.load_telescope()
        self.load_quality_rules()
        try: self.load_sources()
        except: pass
        try: self.load_opacity()
        except: pass

    def load_path(self):

        ldf = 1 # load from file
        config_fname = self.cfn

        if not os.path.isfile(config_fname):
            ldf = 0
        else:
            with open(config_fname) as config_fobj:
                line = config_fobj.read()
            n = len(line)
            line = line[0:n-1]
            line = re.sub('\n', '^Z', line)

            line = self.key_matcher('DATA PATH', line)
            if line:
                line = [e.split()[0:2] for e in line]

                tp0 = dict(line)
                tp1 = {'Raw': 'fits_path',
                        'Fit': 'fit_path'}

                for key in list(tp0.keys()):
                    setattr(self, tp1[key], tp0[key])
                ldf = 1

            else: ldf = 0

        if not ldf:
            self.fits_path = './fits'
            self.fit_path = './fit'

        return ldf


    def write_path(self):

        config_fname = self.cfn

        S = '\n  Raw    %s' %self.fits_path
        S = '%s\n  Fit    %s' %(S, self.fit_path)
        S = '# BEGIN DATA PATH%s' %S
        S = '%s\n# END DATA PATH' %S

        if os.path.isfile(config_fname):
            with open(config_fname) as config_fobj:
                line = config_fobj.read()
            n = len(line)
            line = line[0:n-1]
            line = re.sub('\n', '^Z', line)

            p = re.compile('.*DATA PATH['', ' ']*\^Z')
            line0 = p.search(line)
            if line0:
                line0 = line0.group()
                line0 = re.sub('# BEGIN DATA PATH.*', '', line0)
                line0 = re.sub('\^Z', '\n', line0)
            else:
                line0 = re.sub('\^Z', '\n', line0)

            p = re.compile('# END DATA PATH.*')
            line1 = p.search(line)
            if line1:
                line1 = line1.group()
                if '^Z' in line1:
                    line1 = re.sub('# END DATA PATH['', ' ']*\^Z', '', line1)
                    line1 = re.sub('\^Z', '\n', line1)
                else:
                    line1 = ''
            else:
                line1 = ''

            if line0:
                out_line = '%s\n%s\n%s' %(line0, S, line1)
            else:
                out_line = '%s\n%s' %(S, line1)

        else:
            out_line = S

        with open(config_fname, 'w') as config_fobj:
            print(out_line, file=config_fobj)


    def load_telescope(self):

        ldf = 1 # load from file
        config_fname = self.cfn

        if not os.path.isfile(config_fname):
            ldf = 0
        else:
            with open(config_fname) as config_fobj:
                line = config_fobj.read()
            n = len(line)
            line = line[0:n-1]
            line = re.sub('\n', '^Z', line)

            line = self.key_matcher('TELESCOPE PARAMETERS', line)
            if line:
                line = [e.split()[0:2] for e in line]

                tp0 = dict(line)
                tp1 = {'Telescope': 'tele_name',
                        'Diameter': 'diam',
                        'Obs_lon': 'obs_lon',
                        'Obs_lat': 'obs_lat',
                        'Frequency': 'freq',
                        'Tcal': 'tcal'}

                for key in list(tp0.keys()):
                    if key == 'Telescope':
                        setattr(self, tp1[key], tp0[key])
                    else:
                        setattr(self, tp1[key], float(tp0[key]))
                ldf = 1

            else: ldf = 0

        if not ldf: self.tele_name = 'NSRT'

        return ldf


    def write_telescope(self):

        config_fname = self.cfn

        if os.path.isfile(config_fname):
            with open(config_fname) as config_fobj:
                line = config_fobj.read()
            n = len(line)
            line = line[0:n-1]
            line = re.sub('\n', '^Z', line)

            p = re.compile('.*BEGIN TELESCOPE PARAMETERS['', ' ']*\^Z')
            line0 = p.search(line)
            if line0:
                line0 = line0.group()
                line0 = re.sub('['', ' ']*\^Z# BEGIN TELESCOPE PARAMETERS.*',
                        '', line0)
                line0 = re.sub('\^Z', '\n', line0)
            else:
                line0 = re.sub('\^Z', '\n', line)

            p = re.compile('# END TELESCOPE PARAMETERS.*')
            line1 = p.search(line)
            if line1:
                line1 = line1.group()
                if '^Z' in line1:
                    line1 = re.sub('# END TELESCOPE PARAMETERS['', ' ']*\^Z',
                           '', line1)
                    line1 = re.sub('\^Z', '\n', line1)
                else:
                    line1 = ''
            else:
                line1 = ''

        names = ['Diameter', 'Obs_lon', 'Obs_lat', 'Frequency', 'Tcal']
        cls = ['diam', 'obs_lon', 'obs_lat', 'freq', 'tcal']
        units = ['[m]', '[deg.]', '[deg.]', '[MHz]', '[K]']

        S = '# BEGIN TELESCOPE PARAMETERS'
        S = '%s%s' %(S, '\n  Telescope %11s' %self.tele_name)

        for i in range(5):
            S = '%s\n  %-10s %10s %6s' \
                        %(S, names[i], str(getattr(self, cls[i])), units[i])
        S = '%s%s' %(S, '\n# END TELESCOPE PARAMETERS')

        if line1:
            out_line = '%s\n%s\n%s' %(line0, S, line1)
        else:
            out_line = '%s\n%s' %(line0, S)

        with open(config_fname, 'w') as config_fobj:
            print(out_line, file=config_fobj)


    def load_quality_rules(self):

        ldf = 1 # load from file
        config_fname = self.cfn

        if not os.path.isfile(config_fname):
            ldf = 0
        else:
            config_fobj = open(config_fname)
            line = config_fobj.read()
            n = len(line)
            line = line[0:n-1]
            line = re.sub('\n', '^Z', line)

            line = self.key_matcher('FLAG OPTIONS', line)
            if line:
                line = [e.split()[0:2] for e in line]

                tp0 = dict(line)
                tp1 = {'Beam_Width': 'beam',
                        'D_Off': 'doff',
                        'D_Width': 'dwidth',
                        'Rela_Err': 'rerr',
                        'D_Symmetry': 'dsymm',
                        'D_Avg': 'davg'}

                for key in list(tp0.keys()):
                    setattr(self, tp1[key], float(tp0[key]))
                ldf = 1

            else: ldf = 0

        if not ldf:
            self.beam = 3E2/self.freq/self.diam*180/np.pi*3600*1.12
            self.doff = 0.15*self.beam
            self.dwidth = 0.1
            self.rerr = 0.1
            self.dsymm = 0.1
            self.davg = 0.15

        return ldf


    def write_quality_rules(self):

        config_fname = self.cfn

        if os.path.isfile(config_fname):
            with open(config_fname) as config_fobj:
                line = config_fobj.read()
            n = len(line)
            line = line[0:n-1]
            line = re.sub('\n', '^Z', line)

            p = re.compile('.*BEGIN FLAG OPTIONS['', ' ']*\^Z')
            line0 = p.search(line)
            if line0:
                line0 = line0.group()
                #line0 = re.sub('['', ' ']*\^Z# BEGIN FLAG OPTIONS['', ' ']*\^Z', '', line0)
                line0 = re.sub('['', ' ']*\^Z# BEGIN FLAG OPTIONS.*', '', line0)
                line0 = re.sub('\^Z', '\n', line0)
            else:
                line0 = re.sub('\^Z', '\n', line)

            p = re.compile('# END FLAG OPTIONS.*')
            line1 = p.search(line)
            if line1:
                line1 = line1.group()
                if '^Z' in line1:
                    line1 = re.sub('# END FLAG OPTIONS['', ' ']*\^Z', '', line1)
                    line1 = re.sub('\^Z', '\n', line1)
                else:
                    line1 = ''
            else:
                line1 = ''

        fpk = ['Beam_Width', 'D_Off', 'D_Width', 'Rela_Err', \
                'D_Symmetry', 'D_Avg']
        fp = {'Beam_Width': 'beam',
                'D_Off': 'doff',
                'D_Width': 'dwidth',
                'Rela_Err': 'rerr',
                'D_Symmetry': 'dsymm',
                'D_Avg': 'davg'}
        S = '# BEGIN FLAG OPTIONS'
        for key in fpk:
            S = '%s\n  %-10s %6.2f' %(S, key, getattr(self, fp[key]))
        S = '%s\n# END FLAG OPTIONS' %S

        if line1:
            out_line = '%s\n%s\n%s' %(line0, S, line1)
        else:
            out_line = '%s\n%s' %(line0, S)

        with open(config_fname, 'w') as config_fobj:
            print(out_line, file=config_fobj)


    def load_sources(self):

        ldf = 1 # load from file
        config_fname = self.cfn

        if not os.path.isfile(config_fname):
            ldf = 0

        else:
            with open(config_fname) as config_fobj:
                line = config_fobj.read()
            n = len(line)
            line = line[0:n-1]
            line = re.sub('\n', '^Z', line)

            pri_line = self.key_matcher('PRI.CAL.', line)
            #sec_line = self.key_matcher('SEC.CAL.', line)
            gain_line = self.key_matcher('GAIN CAL.', line)
            time_line = self.key_matcher('TIME CAL.', line)
            target_line = self.key_matcher('TARGETS', line)

            if pri_line != []:
                self.pri_cal = [e.split()[0] for e in pri_line]
            else:
                self.pri_cal = []
                ldf = 0

#            if sec_line != []:
#                self.sec_cal = [e.split()[0] for e in sec_line]
#            else:
#                self.sec_cal = []

            if gain_line != []:
                self.gain_cal = [e.split()[0] for e in gain_line]
            else:
                self.gain_cal = []

            if time_line != []:
                self.time_cal = [e.split()[0] for e in time_line]
            else:
                self.time_cal = []

            if target_line != []:
                self.target = [e.split()[0] for e in target_line]
            else:
                self.target = []

            self.all_src = self.pri_cal + self.gain_cal + \
                    self.time_cal + self.target
#            self.all_src = self.pri_cal + self.sec_cal + self.gain_cal + \
#                    self.time_cal + self.target

            self.all_src.sort()
            self.all_src = list(set(self.all_src))

        return ldf


    def write_sources(self):

        config_fname = self.cfn

        if os.path.isfile(config_fname):
            with open(config_fname) as config_fobj:
                line = config_fobj.read()
            n = len(line)
            line = line[0:n-1]
            line = re.sub('\n', '^Z', line)

            p = re.compile('.*BEGIN PRI.CAL.['', ' ']*\^Z')
            line0 = p.search(line)
            if line0:
                line0 = line0.group()
                #line0 = re.sub('['', ' ']*\^Z# BEGIN PRI.CAL.['', ' ']*\^Z', '', line0)
                line0 = re.sub('['', ' ']*\^Z# BEGIN PRI.CAL.*', '', line0)
                line0 = re.sub('\^Z', '\n', line0)
            else: line0 = re.sub('\^Z', '\n', line)

        p = re.compile('# END TARGETS.*')
        line1 = p.search(line)
        if line1:
            line1 = line1.group()
            if '^Z' in line1:
                line1 = re.sub('# END TARGETS['', ' ']*\^Z', '', line1)
                line1 = re.sub('\^Z', '\n', line1)
            else:
                line1 = ''
        else:
            line1 = ''

        if self.pri_cal != []:
            S1 = ('\n  ').join([e for e in self.pri_cal])
            S1 = '# %s\n  %s\n# %s' %('BEGIN PRI.CAL.', S1, 'END PRI.CAL.')
        else:
            S1 = '# %s\n# %s' %('BEGIN PRI.CAL.', 'END PRI.CAL.')

#        if self.sec_cal != []:
#            S2 = ('\n  ').join([e for e in self.sec_cal])
#            S2 = '# %s\n  %s\n# %s' %('BEGIN SEC.CAL.', S2, 'END SEC.CAL.')
#        else:
#            S2 = '# %s\n# %s' %('BEGIN SEC.CAL.', 'END SEC.CAL.')

        if self.gain_cal != []:
            S3 = ('\n  ').join([e for e in self.gain_cal])
            S3 = '# %s\n  %s\n# %s' %('BEGIN GAIN CAL.', S3, 'END GAIN CAL.')
        else:
            S3 = '# %s\n# %s' %('BEGIN GAIN CAL.', 'END GAIN CAL.')

        if self.time_cal != []:
            S4 = ('\n  ').join([e for e in self.time_cal])
            S4 = '# %s\n  %s\n# %s' %('BEGIN TIME CAL.', S4, 'END TIME CAL.')
        else:
            S4 = '# %s\n# %s' %('BEGIN TIME CAL.', 'END TIME CAL.')

        if self.target != []:
            S5 = ('\n  ').join([e for e in self.target])
            S5 = '# %s\n  %s\n# %s' %('BEGIN TARGETS', S5, 'END TARGETS')
        else:
            S5 = '# %s\n# %s' %('BEGIN TARGETS', 'END TARGETS')

        if line1:
#            out_line = '%s\n%s\n%s\n%s\n%s\n%s\n%s' \
#                    %(line0, S1, S2, S3, S4, S5, line1)
            out_line = '%s\n%s\n%s\n%s\n%s\n%s' \
                    %(line0, S1, S3, S4, S5, line1)
        else:
            #out_line = '%s\n%s\n%s\n%s\n%s\n%s' %(line0, S1, S2, S3, S4, S5)
            out_line = '%s\n%s\n%s\n%s\n%s' %(line0, S1, S3, S4, S5)

        with open(config_fname, 'w') as config_fobj:
            print(out_line, file=config_fobj)


    def load_opacity(self):

        ldf = 1 # load from file
        config_fname = self.cfn

        if not os.path.isfile(config_fname):
            ldf = 0
        else:
            with open(config_fname) as config_fobj:
                line = config_fobj.read()
            n = len(line)
            line = line[0:n-1]
            line = re.sub('\n', '^Z', line)

            line = self.key_matcher('AIR OPACITY', line)
            if line != []:
                line = [e.split()[0:2] for e in line]

                tp0 = dict(line)
                tp1 = {'Tground': 'Tgr',
                        'T0': 'Trec',
                        'Tau0': 'Tau0'}

                for key in list(tp0.keys()):
                    setattr(self, tp1[key], float(tp0[key]))
                ldf = 1

            else:
                self.Tgr = 0.0
                self.Trec = 0.0
                self.Tau0 = 0.0
                ldf = 0

        return ldf


    def write_opacity(self):

        config_fname = self.cfn

        if os.path.isfile(config_fname):
            with open(config_fname) as config_fobj:
                line = config_fobj.read()
            n = len(line)
            line = line[0:n-1]
            line = re.sub('\n', '^Z', line)

            p = re.compile('.*AIR OPACITY['', ' ']*\^Z')
            line0 = p.search(line)
            if line0:
                line0 = line0.group()
                line0 = re.sub('['', ' ']*\^Z# BEGIN AIR OPACITY.*', '', line0)
                line0 = re.sub('\^Z', '\n', line0)
            else:
                line0 = re.sub('\^Z', '\n', line)

            p = re.compile('# END AIR OPACITY.*')
            line1 = p.search(line)
            if line1:
                line1 = line1.group()
                if '^Z' in line1:
                    line1 = re.sub('# END AIR OPACITY['', ' ']*\^Z', '', line1)
                    line1 = re.sub('\^Z', '\n', line1)
                else:
                    line1 = ''
            else:
                line1 = ''

            S = '\n  Tground %10.1f' %self.Tgr
            S = '%s\n  T0 %15.2f' %(S, self.Trec)
            S = '%s\n  Tau0 %13.4f' %(S, self.Tau0)
            S = '# BEGIN AIR OPACITY%s' %S
            S = '%s\n# END AIR OPACITY' %S

        if line1:
            out_line = '%s\n%s\n%s' %(line0, S, line1)
        else:
            out_line = '%s\n%s' %(line0, S)

        with open(config_fname, 'w') as config_fobj:
            print(out_line, file=config_fobj)


    def key_matcher(self, key, line):

        out_line = re.sub('.*BEGIN %s['', ' ']*\^Z' %key, '', line)
        out_line = re.sub('['', ' ']*\^Z#.END %s.*' %key, '', out_line)
        out_line = out_line.split('^Z')
        if out_line[0][0] == '#':
            out_line = []
        return out_line


class RawData(object):

    def __init__(self, conf, fname, **kwargs):

        args = {'chan': '0.5*(L+R)',
                'nchan': 4,
                'nphase': 4,
                'despike': 0,
                'cut': [0.1, 0.1],
                'window': [1.35, 1.35],
                'npeak': 1,
                'bo': 1,
                'smodel': [0.0, 0.0, 0.0, 0.0, 0.0]}

        for key in list(args.keys()):
            if key in kwargs:
                setattr(self, key, kwargs[key])
            else:
                setattr(self, key, args[key])

        self.conf = conf # an instance of Config class
        self.fname = os.path.join(self.conf.fits_path, fname)
        self.fitsobj = pyfits.open(self.fname)

        self.src = self.fitsobj[0].header['object']
        self.scann = self.fitsobj[1].header['scannum']
        self.slon = self.fitsobj[1].header['crval1']*np.pi/180.0
        self.slat = self.fitsobj[1].header['crval2']*np.pi/180.0
        self.nsub = (len(self.fitsobj)-3)/3

        self.sig_onln = self.fitsobj[2].header['sig_onln']
        self.ref_onln = self.fitsobj[2].header['ref_onln']
        self.sig_pol = self.fitsobj[2].header['sig_pol']
        self.ref_pol = self.fitsobj[2].header['ref_pol']

    @staticmethod
    def extract_sig(sig_mode, arg1, arg2):
        """
        Define signal/reference modes, for data extraction.
        """

        return sig_mode.replace(')', '-1]').\
            replace('(', '%s[' %arg1).\
            replace(',', '-1+%s,' %arg2)


    @staticmethod
    def total_int(sig, cal):
        """
        extract total intensity

        Returns
        -------
        source and noise temperature
        """

        n = len(sig)
        if n > 1:
            x = np.arange(n)
            fit_val = ml.fit_poly(x, cal, 1, 0)[0]
            cal_corr = np.polyval(fit_val, x)

            return (sig-cal)/(2*cal_corr), cal_corr/2.

        else:
            return (sig-cal)/(2*cal), cal/2.


    @staticmethod
    def total_pol(chan3, cal3, chan4, cal4):
        """
        polarization parameters u and q

        Returns
        -------
        polarization strength and polarization angle
        """

        chan3_corr = (chan3 - cal3)/4.
        chan4_corr = (chan4 - cal4)/4.
        cal3_corr, cal4_corr = cal3/2., cal4/2.

        n = len(chan3_corr)
        x = np.arange(n)
        fit3 = ml.fit_poly(x, cal3_corr, 1, 0)[0]
        fit4 = ml.fit_poly(x, cal4_corr, 1, 0)[0]
        cal3_corr = np.polyval(fit3, x)
        cal4_corr = np.polyval(fit4, x)
        del fit3, fit4

        chan_corr = (chan3_corr**2. + chan4_corr**2.)**0.5
        theta = np.arctan(cal3_corr/cal4_corr)
        psin = np.sin(theta)
        pcos = np.cos(theta)
        s3 = chan3*pcos - chan4*psin
        s4 = chan3*psin + chan4*pcos

        return s3/chan_corr, s4/chan_corr, chan_corr, 90*np.mean(theta)/np.pi


    def load_data_par(self, subn):

        """
        function that manages the 'datapar' extension

        """
        fits_obj = self.fitsobj[3*(subn-1)+4]  # datapar

        try:
            self.scandir = fits_obj.header['scandir']

        except:
            if subn < self.nsub/2:
                self.scandir = 'ALON'
            else:
                self.scandir = 'ALAT'

        self.azi = fits_obj.data.field('azimuth')
        sublen = len(self.azi)
        if sublen < 24:
            return False

        self.azi = np.mean(self.azi[1:sublen-1])
        self.elv = np.mean(fits_obj.data.field('elevatio')[1:sublen-1])
        self.prg = np.mean(fits_obj.data.field('parangle')[1:sublen-1])

        self.mjd_array = fits_obj.data.field('mjd')[1:sublen-1]
        self.mjd = np.mean(self.mjd_array)

        if self.scandir== 'ALON':
            self.offset_array = 3600.*fits_obj.data.field('longoff')[1:sublen-1]

        elif self.scandir== 'ALAT':
            self.offset_array = 3600.*fits_obj.data.field('latoff')[1:sublen-1]

        else:
            print('Bad Data Format! Please Check!')
            return

        ut = al.mjd2utc(self.mjd)
        tmp, tmp, self.sun_ang = \
        al.ra_dec2az_el(self.conf.obs_lon, self.conf.obs_lat,
                self.slon, self.slat, ut, 1)

        del fits_obj, sublen, ut, tmp

        return True


    def load_array_data(self, subn):

        """
        function that manages the 'arraydata' extension
        """

        fits_obj = self.fitsobj[3*(subn-1)+5].data  # arraydata
        time = fits_obj.field('mjd')
        sig = fits_obj.field('data')

        sublen = len(sig)/4  # subscan length, 4 phases
        # channel 1,2,3,4 and Tcal of theses channels
        self.obs_mjd = np.zeros(sublen)
        self.chan_val = \
            np.zeros(self.nchan*self.nphase*sublen).\
                reshape(self.nchan, self.nphase, sublen)

        ### create empty arrays to restore intensity and polarization data
        self.int_array = np.zeros(sublen)
        self.cal_array = np.zeros(sublen)

        self.pol1_array = np.zeros(sublen)
        self.pol2_array = np.zeros(sublen)
        self.pol1_cal_array = np.zeros(sublen)
        self.pol2_cal_array = np.zeros(sublen)

        time = time.reshape(sublen, self.nphase)
        self.obs_mjd = np.array([np.mean(e) for e in time])

        sig = sig.reshape(sublen,
                          self.nphase,
                          self.nchan).transpose(2, 1, 0)

        for i in range(self.nchan/2):
            self.chan_val[i][0] = \
                eval(self.extract_sig(self.sig_onln, 'sig', 'i'))
            self.chan_val[i][1] = \
                eval(self.extract_sig(self.ref_onln, 'sig', 'i'))
            self.chan_val[i+2][0] = \
                eval(self.extract_sig(self.sig_pol, 'sig', 'i+2'))
            self.chan_val[i+2][1] = \
                eval(self.extract_sig(self.ref_pol, 'sig', 'i+2'))


    def align_data(self):

        """
        function that aligns the 'datapar' and 'arraydata' extensions

        """

        # possible trim of the data
        flt = (self.obs_mjd > min(self.mjd_array))*\
               (self.obs_mjd < max(self.mjd_array))
        self.obs_mjd = self.obs_mjd[flt]
        sublen = len(flt)
        self.sublen = len(self.obs_mjd)

        flt = np.array((flt.tolist()*self.nchan*self.nphase)).\
            reshape(self.nchan,
                    self.nphase,
                    sublen)

        self.chan_val = self.chan_val[flt].reshape(self.nchan,
                                                    self.nphase,
                                                    self.sublen)

        self.scan_offsets = np.interp(self.obs_mjd,
                                     self.mjd_array,
                                     self.offset_array)
        self.scan_amp = np.array([])
        self.cal_array_align = np.array([])
        self.pol1_array_align = np.array([])
        self.pol2_array_align = np.array([])
        self.pol_cal_array_align = np.array([])
        self.pa = 0

        sig_sky, sig_ref = [0, 0], [0, 0]

        for i in range(self.nchan/2):
            sig_sky[i], sig_ref[i] = \
                self.total_int(self.chan_val[i][0], self.chan_val[i][1])
        self.pol1_array_align, self.pol2_array_align, \
        self.pol_cal_array_align, self.pa = self.total_pol(self.chan_val[2][0],
                                                 self.chan_val[2][1],
                                                 self.chan_val[3][0],
                                                 self.chan_val[3][1])

        p_sky = re.sub('L', 'sig_sky[0]', self.chan)
        p_sky = re.sub('R', 'sig_sky[1]', p_sky)
        p_ref = re.sub('L', 'sig_ref[0]', self.chan)
        p_ref = re.sub('R', 'sig_ref[1]', p_ref)
        self.scan_amp = eval(p_sky)
        self.cal_array_align = eval(p_ref)

#        if self.despike:
#            self.scan_amp = ml.DeSpike(self.scan_amp, 5, self.despike)
#            self.cal_array_align = \
#                    ml.DeSpike(self.cal_array_align, 5, self.despike)
#            self.pol1_array_align = \
#                    ml.DeSpike(self.pol1_array_align, 5, self.despike)
#            self.pol2_array_align = \
#                    ml.DeSpike(self.pol2_array_align, 5, self.despike)
#            self.pol_cal_array_align = \
#                    ml.DeSpike(self.pol_cal_array_align, 5, self.despike)

        self.pol_cal = np.mean(self.pol_cal_array_align)
        self.int_cal = np.mean(self.cal_array_align)
        self.sublen = len(self.scan_amp)


    def fit_sub(self, subn):

        status = self.load_data_par(subn)
        if not status:
            return status
        self.load_array_data(subn)
        self.align_data()

        self.scan_amp*=self.conf.tcal

        ###### first fit  #############
        span = max(self.scan_offsets) - min(self.scan_offsets)
        dcut = self.scan_offsets >= min(self.scan_offsets) + span*self.cut[0]
        dcut*= self.scan_offsets <= max(self.scan_offsets) - span*self.cut[1]

        del span

        exp_len = self.npeak*3 + self.bo + 1
        true_len = len(self.smodel)
        if exp_len > true_len:
            diff_len = exp_len - true_len
            rp0 = np.append(self.smodel[:self.npeak*3], np.zeros(diff_len))
            rp0 = np.append(rp0, self.smodel[self.npeak*3:])

        elif exp_len == true_len:
            rp0 = np.array(self.smodel)*1.0

        x = self.scan_offsets[dcut]
        y = self.scan_amp[dcut]
        if len(x) < 2+self.bo+self.npeak*3:
            dcut = np.arange(int(0.2*len(self.scan_amp)),
                            int(0.8*len(self.scan_amp)))
            x = self.scan_offsets[dcut]
            y = self.scan_amp[dcut]
        span = max(x) - min(x)
        if 0 in rp0[3*self.npeak:]:
            rp0[3*self.npeak:] = ml.fit_poly(x, y, 0, self.bo)[0]

        gs = ml.Gauss(self.npeak, self.bo)
        for i in range(exp_len-self.bo-1):
            yclean = y - gs.fGauss(rp0, x)
            if rp0[i] == 0 and i%3 == 0:
                rp0[i] = (max(yclean)-min(yclean))*0.5
            elif rp0[i] == 0 and i%3 == 1:
                rp0[i] = x[yclean.argmax()]
            elif rp0[i] == 0 and i%3 == 2:
                rp0[i] = 3E2/self.conf.freq/self.conf.diam*180/np.pi*3600*1.12
        self.smodel = rp0*1.0

#        if np.array(self.smodel).sum() == 0.0:
#            x = self.scan_offsets[dcut]
#            y = self.scan_amp[dcut]
#            span = max(x) - min(x)
#            rp0 = np.zeros(3*self.npeak+self.bo+1)
#            rp0[3*self.npeak:] = ml.fit_poly(x, y, 0, self.bo)[0]
#
#            for i in xrange(self.npeak):
#                gs = ml.Gauss((len(rp0)-2)/3, self.bo)
#                rp0[i*3+2] = \
#                        3E2/self.conf.freq/self.conf.diam*180/np.pi*3600*1.12
#                yclean = y - gs.fGauss(rp0, x)
#                rp0[i*3] = (max(yclean) - min(yclean))*0.5
#
#                yclean /= np.fabs(x)+rp0[i*3+2]*2.0
#                rp0[i*3+1] = x[yclean.argmax()]
#            self.smodel = rp0*1.0

        del x, y, span, rp0, gs, yclean

        gauss = ml.Gauss(self.npeak, self.bo)

        self.pm, self.pe = gauss.fit_window(self.scan_offsets[dcut],
                self.scan_amp[dcut], self.smodel, self.conf.beam,
                self.window, self.despike)

        self.bpm = self.pm[3*self.npeak:3*self.npeak+self.bo+1]
        self.tsys = self.bpm[-1]
#        self.smodel = np.zeros(3*self.npeak+self.bo+1)

        return status


    def fit_data(self):

        ffit = open(os.path.join(self.conf.fit_path,
            '%04d.fit' %self.scann), 'w')
        fdat = open(os.path.join(self.conf.fit_path,
            '%04d.dat' %self.scann), 'w')

        print('%04d   %s' % (self.scann, self.src))

        print('# %04d %s' % (self.scann, self.src), file=ffit)
        print('#', file=ffit)
        print('# %7s %8s %7s %6s %7s %12s %9s %9s %11s ' \
                        '%12s %5s %7s %7s %11s %6s %8s %11s' \
                        % ('Source', 'Scan', 'direct', 'Subsc', 'MJD', 'Azi',
                           'Elv', 'Amp', 'err','Offset', 'err', 'FWHM', 'err',
                           'Parangle', 'Tsys', 'Cal', 'SunAng'), file=ffit)
        print('#', file=ffit)

        print('#! %04d %s' %(self.scann, self.src), file=fdat)
        print('#---------------', file=fdat)

        for subn in range(1, 1+self.nsub):
            status = self.fit_sub(subn)
            if not status:
                continue
            #self.fit_sub(7)

            tx, ty = self.scan_offsets, self.scan_amp
            sublen = len(tx)
            blstr = ty - np.polyval(self.bpm, tx)
            residual = ml.Gauss(1, self.bo).fResidual(self.pm, tx, ty)
            self.rms = np.std(residual)

            print('  %6s %4d %12.5e %7.2f %8.2f %8.3f' \
                    % (self.scandir, subn, self.pm[0],
                        self.pm[1], self.pm[2], self.tsys))

            print('  %-11s %04d %6s %4d %15.8f ' \
                    '%9.3f %8.3f %12.5e %12.5e %7.2f ' \
                    '%6.2f %8.2f %7.2f %8.3f %8.3f %10.3f %9.3f' \
                    % (self.src, self.scann, self.scandir,
                        subn, self.mjd, self.azi, self.elv,
                        self.pm[0], self.pe[0],
                        self.pm[1], self.pe[1],
                        self.pm[2], self.pe[2],
                        self.prg, self.tsys,
                        self.int_cal, self.sun_ang), file=ffit)

            print('#! SUBSCAN %02d %s %d' \
                    %(subn, self.scandir, sublen), file=fdat)

            print('#%16s %9s %8s %8s %8s %12s %9s %9s %9s %8s' \
                    %('MJD', 'Azim.', 'Elv.', 'par.Ang.', 'Tsys',
                      'dTsys', 'TcalI', 'TcalP', 'RMS', 'Sun.Ang'), file=fdat)
            print('#!%15.8f %9.3f %8.3f %8.3f %8.3f ' \
                           '%+12.3E %9.3f %9.3f %9.2E %8.3f' \
                           %(self.mjd, self.azi, self.elv,
                             self.prg, self.tsys, self.bpm[-2],
                             self.int_cal,self.pol_cal, self.rms,
                             self.sun_ang), file=fdat)

            print('# %14s %23s %26s %13s %19s %21s' \
                    %('offsets', 'amplitudes', 'baseline_subtracted',
                        'residual', 'pol_u', 'pol_q'), file=fdat)
            print('#'+'-'*125, file=fdat)

            for j in range(sublen):
                print('%20.9f %20.9f %20.9f %20.9f %20.9f %20.9f'\
                    %(tx[j], ty[j], blstr[j], residual[j],
                      self.pol1_array_align[j], self.pol2_array_align[j]), file=fdat)

        ffit.close()
        fdat.close()


if __name__ == '__main__':
    pass

