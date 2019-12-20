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

"""
GUI for View.
"""

__author__ = 'Jun Liu'
__copyright__ = 'Copyright (c) 2015-2020 Jun Liu <liuj@xao.ac.cn, jliu@mpifr-bonn.mpg.de>'
__version__ = 2.1

import wx
import os
import sys
import re
from ppgplot import *
import numpy as np
from pubsub import pub
import Misc
import ObsData
import Calibration


class Fitting(wx.Dialog):
    def __init__(self, parent, *args, **kwds):

        kwds['style'] = wx.DEFAULT_DIALOG_STYLE
        wx.Dialog.__init__(self, parent, *args, **kwds)
        self.cfg = parent.drp.cfg
        self.drp = parent.drp
        self.turn_gray = parent.turn_gray

        hpbw = 3E2/self.cfg.freq/self.cfg.diam*180/np.pi*3600*1.12
        self.tc_fits = wx.TextCtrl(self, wx.ID_ANY, self.cfg.fits_path)
        self.tc_fit = wx.TextCtrl(self, wx.ID_ANY, self.cfg.fit_path)
        self.tc_amp = wx.TextCtrl(self, wx.ID_ANY, '0.0')
        self.tc_off = wx.TextCtrl(self, wx.ID_ANY, '0.1')
        self.tc_hpbw = wx.TextCtrl(self, wx.ID_ANY, '%0.1f' %hpbw)
        self.tc_tsys = wx.TextCtrl(self, wx.ID_ANY, '20.0')
        self.tc_tcal = wx.TextCtrl(self, wx.ID_ANY, '%0.1f' %self.cfg.tcal)
        self.tc_npeak = wx.TextCtrl(self, wx.ID_ANY, '1')
        self.tc_bo = wx.TextCtrl(self, wx.ID_ANY, '1')
        self.tc_chan = wx.TextCtrl(self, wx.ID_ANY, '0.5*(L+R)')
        self.tc_despike = wx.TextCtrl(self, wx.ID_ANY, '0')
        self.tc_cut0 = wx.TextCtrl(self, wx.ID_ANY, '0.1')
        self.tc_cut1 = wx.TextCtrl(self, wx.ID_ANY, '0.1')
        # self.tc_win0 = wx.TextCtrl(self, wx.ID_ANY, '1.35')
        # self.tc_win1 = wx.TextCtrl(self, wx.ID_ANY, '1.35')
        self.tc_win0 = wx.TextCtrl(self, wx.ID_ANY, '1.25')
        self.tc_win1 = wx.TextCtrl(self, wx.ID_ANY, '1.25')
        self.tc_scan0 = wx.TextCtrl(self, wx.ID_ANY, '0')
        self.tc_scan1 = wx.TextCtrl(self, wx.ID_ANY, '9999')

        self.Bind(wx.EVT_CLOSE, self.on_cancel)

        self.__do_layout()

    def __do_layout(self):

        self.SetTitle('Fitting')
        sz0 = wx.BoxSizer(wx.VERTICAL)

        sz_pth = Misc.wx_hsbs(self)
        sz_mod = Misc.wx_hsbs(self)
        sz_cfg = Misc.wx_hsbs(self)
        sz_bt = Misc.wx_hsbs(self)

        sz_cfg2 = wx.BoxSizer(wx.VERTICAL)
        sz_scan = wx.BoxSizer(wx.HORIZONTAL)
        sz_win = wx.BoxSizer(wx.HORIZONTAL)
        sz_cut = wx.BoxSizer(wx.HORIZONTAL)
        sz_cfg1 = wx.BoxSizer(wx.VERTICAL)
        sz_mod2 = wx.BoxSizer(wx.VERTICAL)
        sz_mod1 = wx.BoxSizer(wx.VERTICAL)
        sz_pth2 = wx.BoxSizer(wx.VERTICAL)
        sz_pth1 = wx.BoxSizer(wx.VERTICAL)
        sz_pth0 = wx.BoxSizer(wx.VERTICAL)
        sz_pth0.Add(Misc.wx_st(self, 'Raw'), 1, wx.ALIGN_RIGHT, 0)
        sz_pth0.Add(Misc.wx_st(self, 'Fit'), 1, wx.ALIGN_RIGHT, 0)
        sz_pth.Add(sz_pth0, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 5)
        sz_pth1.Add(self.tc_fits, 1, wx.EXPAND, 0)
        sz_pth1.Add(self.tc_fit, 1, wx.EXPAND | wx.ALIGN_CENTER_VERTICAL, 0)
        sz_pth.Add(sz_pth1, 1, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz_pth2.Add(Misc.wx_bt(self, 'Browse', self.on_browse_fits, \
                style=wx.BU_EXACTFIT), 1, 0, 0)
        sz_pth2.Add(Misc.wx_bt(self, 'Browse', self.on_browse_fit, \
                style=wx.BU_EXACTFIT), 1, 0, 0)
        sz_pth.Add(sz_pth2, 0, wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz0.Add(sz_pth, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz_mod1.Add(Misc.wx_st(self, 'Amplitude'), 1, wx.ALIGN_RIGHT, 0)
        sz_mod1.Add(Misc.wx_st(self, 'Offset'), 1, wx.ALIGN_RIGHT, 0)
        sz_mod1.Add(Misc.wx_st(self, 'HPBW'), 1, wx.ALIGN_RIGHT, 0)
        sz_mod1.Add(Misc.wx_st(self, 'Tsys'), 1, wx.ALIGN_RIGHT, 0)
        sz_mod1.Add(Misc.wx_st(self, 'Tcal'), 1, wx.ALIGN_RIGHT, 0)
        sz_mod.Add(sz_mod1, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 5)
        sz_mod2.Add(self.tc_amp, 1, wx.LEFT | wx.EXPAND, 30)
        sz_mod2.Add(self.tc_off, 1, wx.LEFT | wx.EXPAND, 30)
        sz_mod2.Add(self.tc_hpbw, 1, wx.LEFT | wx.EXPAND, 30)
        sz_mod2.Add(self.tc_tsys, 1, wx.LEFT | wx.EXPAND, 30)
        sz_mod2.Add(self.tc_tcal, 1, wx.LEFT | wx.EXPAND, 30)
        sz_mod.Add(sz_mod2, 2, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz0.Add(sz_mod, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz_cfg1.Add(Misc.wx_st(self, 'Num. of Peaks'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg1.Add(Misc.wx_st(self, 'Baseline Order'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg1.Add(Misc.wx_st(self, 'Channel'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg1.Add(Misc.wx_st(self, 'De-Noise'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg1.Add(Misc.wx_st(self, 'Data Cut'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg1.Add(Misc.wx_st(self, 'Window'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg1.Add(Misc.wx_st(self, 'Scan Num.'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg.Add(sz_cfg1, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 5)
        sz_cfg2.Add(self.tc_npeak, 1, wx.LEFT | wx.EXPAND, 30)
        sz_cfg2.Add(self.tc_bo, 1, wx.LEFT | wx.EXPAND, 30)
        sz_cfg2.Add(self.tc_chan, 1, wx.LEFT | wx.EXPAND, 30)
        sz_cfg2.Add(self.tc_despike, 1, wx.LEFT | wx.EXPAND, 30)
        sz_cut.Add(self.tc_cut0, 1, wx.EXPAND, 0)
        sz_cut.Add(self.tc_cut1, 1, wx.EXPAND, 0)
        sz_cfg2.Add(sz_cut, 1, wx.LEFT | wx.EXPAND, 30)
        sz_win.Add(self.tc_win0, 1, wx.EXPAND, 0)
        sz_win.Add(self.tc_win1, 1, wx.EXPAND, 0)
        sz_cfg2.Add(sz_win, 1, wx.LEFT | wx.EXPAND, 30)
        sz_scan.Add(self.tc_scan0, 1, wx.EXPAND, 0)
        sz_scan.Add(self.tc_scan1, 1, wx.EXPAND, 0)
        sz_cfg2.Add(sz_scan, 1, wx.LEFT | wx.EXPAND, 30)
        sz_cfg.Add(sz_cfg2, 1, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz0.Add(sz_cfg, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)

        sz_bt.Add(Misc.wx_bt(self, 'Cancel', self.on_cancel,
                style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Reset', self.on_reset, \
                style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Go', self.on_go, style=wx.BU_EXACTFIT), \
                1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz0.Add(sz_bt, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        self.SetSizer(sz0)
        sz0.Fit(self)
        self.Layout()


    def on_browse_fits(self, evt):

        dialog = wx.DirDialog(self, 'Set Raw FITS Location', \
                defaultPath = self.tc_fits.Value+'/..',
                style = wx.DD_DEFAULT_STYLE)
        if dialog.ShowModal() == wx.ID_OK:
            self.tc_fits.SetValue(dialog.GetPath())
            self.cfg.fits_path = self.tc_fits.Value
            self.cfg.write_path()
            dialog.Destroy()

    def on_browse_fit(self, evt):

        dialog = wx.DirDialog(self, 'Set resultant fit Location', \
                defaultPath = './', style = wx.DD_DEFAULT_STYLE)
        if dialog.ShowModal() == wx.ID_OK:
            self.tc_fit.SetValue(dialog.GetPath())
            self.cfg.fit_path = self.tc_fit.Value
            self.cfg.write_path()
            dialog.Destroy()

    def on_reset(self, evt):

        hpbw = 3E2/self.cfg.freq/self.cfg.diam*180/np.pi*3600*1.12
        self.tc_fits.SetValue('./fits')
        self.tc_fit.SetValue('./fit')
        self.tc_amp.SetValue('0.0')
        self.tc_off.SetValue('0.1')
        self.tc_hpbw.SetValue('%0.1f' %hpbw)
        self.tc_tsys.SetValue('20.0')
        self.tc_tcal.SetValue('%0.1f' %self.cfg.tcal)

        self.tc_npeak.SetValue('1')
        self.tc_bo.SetValue('1')
        self.tc_chan.SetValue('0.5*(L+R)')
        self.tc_despike.SetValue('0')
        self.tc_cut0.SetValue('0.1')
        self.tc_cut1.SetValue('0.1')
        # self.tc_win0.SetValue('1.35')
        # self.tc_win1.SetValue('1.35')
        self.tc_win0.SetValue('1.25')
        self.tc_win1.SetValue('1.25')
        self.tc_scan0.SetValue('0')
        self.tc_scan1.SetValue('9999')

    def on_cancel(self, evt):

        wx.CallAfter(pub.sendMessage, 'bt_color', msg='fit')
        self.Destroy()

    def on_go(self, evt):

        self.cfg.write_path()
        smodel= '[%s,%s,%s,0,%s]' %(self.tc_amp.Value, self.tc_off.Value,
                self.tc_hpbw.Value, self.tc_tsys.Value)
        scan = '[%d,%d]' %(int(self.tc_scan0.Value), int(self.tc_scan1.Value))
        chan = self.tc_chan.Value
        despike = float(self.tc_despike.Value)
        cut = '[%f,%f]' %(float(self.tc_cut0.Value), float(self.tc_cut1.Value))
        window = '[%f,%f]' %(float(self.tc_win0.Value),
                float(self.tc_win1.Value))
        npeak = int(self.tc_npeak.Value)
        bo = int(self.tc_bo.Value)

        cl = '-smodel %s -scan %s -chan %s -denoise %f '\
                '-cut %s -window %s -npeak %d -bo %d' \
                %(smodel, scan, chan, despike, cut, window, npeak, bo)
        Misc.NewThread(self.drp.do_fit, cl).start()


class Pointing(wx.Dialog):
    def __init__(self, parent, *args, **kwds):
        kwds["style"] = wx.DEFAULT_DIALOG_STYLE
        wx.Dialog.__init__(self, parent, *args, **kwds)
        self.cfg = parent.drp.cfg
        self.drp = parent.drp
        self.clink = wx.Choice(self, wx.ID_ANY, choices=["0", "1"])
        self.cflag = wx.Choice(self, wx.ID_ANY, choices=["0", "1"])
        self.tc_beam = wx.TextCtrl(self, wx.ID_ANY, '%0.1f' %self.cfg.beam)
        self.tc_off = wx.TextCtrl(self, wx.ID_ANY, '%0.1f' %self.cfg.doff)
        self.tc_width = wx.TextCtrl(self, wx.ID_ANY, str(self.cfg.dwidth))
        self.tc_rerr = wx.TextCtrl(self, wx.ID_ANY, str(self.cfg.rerr))
        self.tc_dsym = wx.TextCtrl(self, wx.ID_ANY, str(self.cfg.dsymm))
        self.tc_davg = wx.TextCtrl(self, wx.ID_ANY, str(self.cfg.davg))

        self.Bind(wx.EVT_CLOSE, self.on_cancel)

        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.SetTitle("Pointing")
        self.clink.SetSelection(1)
        self.cflag.SetSelection(1)

    def __do_layout(self):
        sz0 = wx.BoxSizer(wx.VERTICAL)
        sz1 = Misc.wx_hsbs(self)
        sz2 = Misc.wx_hsbs(self)
        sz_bt = Misc.wx_hsbs(self)

        sz21 = wx.BoxSizer(wx.VERTICAL)
        sz20 = wx.BoxSizer(wx.VERTICAL)
        za11 = wx.BoxSizer(wx.VERTICAL)
        za10 = wx.BoxSizer(wx.VERTICAL)
        za10.Add(Misc.wx_st(self, 'Update Result'), 1, wx.TOP | wx.ALIGN_RIGHT, 3)
        za10.Add(Misc.wx_st(self, 'Quality Control'), 1, wx.TOP | wx.ALIGN_RIGHT, 3)
        sz1.Add(za10, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 5)
        za11.Add(self.clink, 0, wx.EXPAND | wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 5)
        za11.Add(self.cflag, 0, wx.EXPAND | wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 5)
        sz1.Add(za11, 1, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz0.Add(sz1, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz20.Add(Misc.wx_bt(self, 'Beam_Width', self.on_beam,
            style=wx.BU_EXACTFIT), 1, wx.ALIGN_RIGHT, 0)
        sz20.Add(Misc.wx_st(self, 'D_Off'), 1, wx.ALIGN_RIGHT | wx.TOP, 6)
        sz20.Add(Misc.wx_st(self, 'D_Width'), 1, wx.ALIGN_RIGHT | wx.TOP, 2)
        sz20.Add(Misc.wx_st(self, 'Rela_Err'), 1, wx.ALIGN_RIGHT | wx.TOP, 2)
        sz20.Add(Misc.wx_st(self, 'D_Symmetry'), 1, wx.ALIGN_RIGHT | wx.TOP, 0)
        sz20.Add(Misc.wx_st(self, 'D_Avg'), 1, wx.ALIGN_RIGHT | wx.TOP, 0)
        sz2.Add(sz20, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz21.Add(self.tc_beam, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_off, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_width, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_rerr, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_dsym, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_davg, 1, wx.LEFT | wx.EXPAND, 10)
        sz2.Add(sz21, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz0.Add(sz2, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Cancel', self.on_cancel,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Reset', self.on_reset,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Go', self.on_go,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz0.Add(sz_bt, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        self.SetSizer(sz0)
        sz0.Fit(self)
        self.Layout()

    def on_beam(self, evt):

        val = 3E2/self.cfg.freq/self.cfg.diam*180/np.pi*3600*1.12
        self.tc_beam.SetValue('%0.1f' %val)

    def on_cancel(self, evt):

        wx.CallAfter(pub.sendMessage, 'bt_color', msg='point')
        self.Destroy()

    def on_reset(self, evt):

        self.clink.SetSelection(0)
        self.cflag.SetSelection(1)

        val = 3E2/self.cfg.freq/self.cfg.diam*180/np.pi*3600*1.12
        self.tc_beam.SetValue('%0.1f' %val)
        self.tc_off.SetValue('100.0')
        self.tc_width.SetValue('0.1')
        self.tc_rerr.SetValue('0.1')
        self.tc_dsym.SetValue('0.1')
        self.tc_davg.SetValue('0.15')


    def on_go(self, evt):

        link = self.clink.Selection
        flag = self.cflag.Selection

        self.cfg.beam = float(self.tc_beam.Value)
        self.cfg.doff = float(self.tc_off.Value)
        self.cfg.dwidth = float(self.tc_width.Value)
        self.cfg.rerr = float(self.tc_rerr.Value)
        self.cfg.dsymm = float(self.tc_dsym.Value)
        self.cfg.davg = float(self.tc_davg.Value)
        self.cfg.write_quality_rules()

        cl = '-link %d -flag %d' %(link, flag)
        Misc.NewThread(self.drp.do_pointing, cl).start()


class Tau(wx.Dialog):
    def __init__(self, parent, *args, **kwds):
        kwds["style"] = wx.DEFAULT_DIALOG_STYLE
        wx.Dialog.__init__(self, parent, *args, **kwds)
        self.cfg = parent.drp.cfg
        self.drp = parent.drp
        self.tc_tgr = wx.TextCtrl(self, wx.ID_ANY, "0.0")

        self.Bind(wx.EVT_CLOSE, self.on_cancel)

        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.SetTitle("Tau")

    def __do_layout(self):
        sz1 = wx.BoxSizer(wx.VERTICAL)
        sz2 = Misc.wx_hsbs(self)
        sz3 = Misc.wx_hsbs(self)
        sz6 = wx.BoxSizer(wx.VERTICAL)
        sz5 = wx.BoxSizer(wx.VERTICAL)
        sz4 = wx.BoxSizer(wx.VERTICAL)
        sz4.Add(Misc.wx_st(self, 'Ground Temp.'), 1,
            wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        sz2.Add(sz4, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 5)
        sz5.Add(self.tc_tgr, 0, wx.EXPAND, 0)
        sz2.Add(sz5, 1, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz6.Add(Misc.wx_st(self, 'Deg. Cent.'), 1, wx.ALIGN_CENTER_VERTICAL, 0)
        sz2.Add(sz6, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 5)
        sz1.Add(sz2, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz3.Add(Misc.wx_bt(self, 'Cancel', self.on_cancel), 1,
            wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz3.Add(Misc.wx_bt(self, 'Reset', self.on_reset), 1,
            wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz3.Add(Misc.wx_bt(self, 'Go', self.on_go), 1,
            wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz1.Add(sz3, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        self.SetSizer(sz1)
        sz1.Fit(self)
        self.Layout()

    def on_cancel(self, evt):

        wx.CallAfter(pub.sendMessage, 'bt_color', msg='tau')
        self.Destroy()

    def on_reset(self, evt):

        self.tc_tgr.SetValue('0.0')

    def on_go(self, evt):

        self.cfg.Tgr = float(self.tc_tgr.Value)

        def func(argv):
            self.drp.do_get_tau(argv)
            self.drp.do_tau_corr(argv)

        Misc.NewThread(func, '-tgr %f' %self.cfg.Tgr).start()



class GainElv(wx.Dialog):
    def __init__(self, parent, *args, **kwds):
        kwds["style"] = wx.DEFAULT_DIALOG_STYLE
        wx.Dialog.__init__(self, parent, *args, **kwds)
        self.drp = parent.drp
        self.tc_inpth = wx.TextCtrl(self, wx.ID_ANY, "1.tau")
        self.tc_order = wx.TextCtrl(self, wx.ID_ANY, "2")
        self.tc_weight = wx.TextCtrl(self, wx.ID_ANY, "0")
        self.tc_limit = wx.TextCtrl(self, wx.ID_ANY, "3")
        self.tc_error = wx.TextCtrl(self, wx.ID_ANY, "0")

        self.Bind(wx.EVT_CLOSE, self.on_cancel)

        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.SetTitle("Gain Elevation")

    def __do_layout(self):
        sz0 = wx.BoxSizer(wx.VERTICAL)
        sz_bt = Misc.wx_hsbs(self)
        sz2 = Misc.wx_hsbs(self)
        sz21 = wx.BoxSizer(wx.VERTICAL)
        sz20 = wx.BoxSizer(wx.VERTICAL)
        sz20.Add(Misc.wx_st(self, 'Input Data'), 1, wx.ALIGN_RIGHT, 0)
        sz20.Add(Misc.wx_st(self, 'Order'), 1, wx.ALIGN_RIGHT, 0)
        sz20.Add(Misc.wx_st(self, 'Weighting'), 1, wx.ALIGN_RIGHT, 0)
        sz20.Add(Misc.wx_st(self, 'Limit'), 1, wx.ALIGN_RIGHT, 0)
        sz20.Add(Misc.wx_st(self, 'Error'), 1, wx.ALIGN_RIGHT, 0)
        sz2.Add(sz20, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 5)
        sz21.Add(self.tc_inpth, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_order, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_weight, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_limit, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_error, 1, wx.LEFT | wx.EXPAND, 10)
        sz2.Add(sz21, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz0.Add(sz2, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Cancel', self.on_cancel,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Reset', self.on_reset,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Go', self.on_go,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz0.Add(sz_bt, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)

        self.SetSizer(sz0)
        sz0.Fit(self)
        self.Layout()

    def on_cancel(self, evt):

        wx.CallAfter(pub.sendMessage, 'bt_color', msg='gain')
        self.Destroy()

    def on_reset(self, evt):

        self.tc_inpth.SetValue('1.tau')
        self.tc_order.SetValue('2')
        self.tc_weight.SetValue('0')
        self.tc_limit.SetValue('3')
        self.tc_error.SetValue('0')


    def on_go(self, evt):

        inpth = self.tc_inpth.Value
        order = int(self.tc_order.Value)
        weight = int(self.tc_weight.Value)
        limit = float(self.tc_limit.Value)
        error = float(self.tc_error.Value)

        def func(argv):
            self.drp.do_get_gain_curve(argv)
            self.drp.do_gain_corr(argv)

        cl = '-inpth %s -order %d -weight %d -limit %f -err %d' \
                %(inpth, order, weight, limit, error)
        Misc.NewThread(func, cl).start()


class GainTime(wx.Dialog):
    def __init__(self, parent, *args, **kwds):
        kwds["style"] = wx.DEFAULT_DIALOG_STYLE
        wx.Dialog.__init__(self, parent, *args, **kwds)
        self.drp = parent.drp
        self.tc_inpth = wx.TextCtrl(self, wx.ID_ANY, "2.gain")
        self.tc_bs = wx.TextCtrl(self, wx.ID_ANY, "0.08")
        self.tc_error = wx.TextCtrl(self, wx.ID_ANY, "0")

        self.Bind(wx.EVT_CLOSE, self.on_cancel)

        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.SetTitle("Gain Time")

    def __do_layout(self):
        sz0 = wx.BoxSizer(wx.VERTICAL)
        sz_bt = Misc.wx_hsbs(self)
        sz2 = Misc.wx_hsbs(self)
        sz21 = wx.BoxSizer(wx.VERTICAL)
        sz20 = wx.BoxSizer(wx.VERTICAL)
        sz20.Add(Misc.wx_st(self, 'Input Data'), 1, wx.ALIGN_RIGHT, 0)
        sz20.Add(Misc.wx_st(self, 'Bin Size'), 1, wx.ALIGN_RIGHT, 0)
        sz20.Add(Misc.wx_st(self, 'Error'), 1, wx.ALIGN_RIGHT, 0)
        sz2.Add(sz20, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 5)
        sz21.Add(self.tc_inpth, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_bs, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_error, 1, wx.LEFT | wx.EXPAND, 10)
        sz2.Add(sz21, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz0.Add(sz2, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Cancel', self.on_cancel,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Reset', self.on_reset,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Go', self.on_go,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz0.Add(sz_bt, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)

        self.SetSizer(sz0)
        sz0.Fit(self)
        self.Layout()

    def on_cancel(self, evt):

        wx.CallAfter(pub.sendMessage, 'bt_color', msg='time')
        self.Destroy()

    def on_reset(self, evt):

        self.tc_inpth.SetValue('2.gain')
        self.tc_bs.SetValue('0.08')
        self.tc_error.SetValue('0')


    def on_go(self, evt):

        inpth = self.tc_inpth.Value
        bs = float(self.tc_bs.Value)
        error = float(self.tc_error.Value)

        def func(argv):
            self.drp.do_get_time_curve(argv)
            self.drp.do_time_corr(argv)

        cl = '-inpth %s -bs %f -err %d' %(inpth, bs, error)
        Misc.NewThread(func, cl).start()


class AbsFlux(wx.Dialog):
    def __init__(self, parent, *args, **kwds):
        kwds["style"] = wx.DEFAULT_DIALOG_STYLE
        wx.Dialog.__init__(self, parent, *args, **kwds)
        self.drp = parent.drp
        self.tc_inpth = wx.TextCtrl(self, wx.ID_ANY, '3.time')
        self.tc_error = wx.TextCtrl(self, wx.ID_ANY, '1.0')

        self.Bind(wx.EVT_CLOSE, self.on_cancel)

        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.SetTitle("Abs. Flux")

    def __do_layout(self):
        sz0 = wx.BoxSizer(wx.VERTICAL)
        sz_bt = Misc.wx_hsbs(self)
        sz2 = Misc.wx_hsbs(self)
        sz21 = wx.BoxSizer(wx.VERTICAL)
        sz20 = wx.BoxSizer(wx.VERTICAL)
        sz20.Add(Misc.wx_st(self, 'Input Data'), 1, wx.ALIGN_RIGHT, 0)
        sz20.Add(Misc.wx_st(self, 'Error'), 1, wx.ALIGN_RIGHT, 0)
        sz2.Add(sz20, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 5)
        sz21.Add(self.tc_inpth, 1, wx.LEFT | wx.EXPAND, 10)
        sz21.Add(self.tc_error, 1, wx.LEFT | wx.EXPAND, 10)
        sz2.Add(sz21, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz0.Add(sz2, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Cancel', self.on_cancel,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Reset', self.on_reset,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Go', self.on_go,
            style=wx.BU_EXACTFIT), 1, wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)
        sz0.Add(sz_bt, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)

        self.SetSizer(sz0)
        sz0.Fit(self)
        self.Layout()

    def on_cancel(self, evt):

        wx.CallAfter(pub.sendMessage, 'bt_color', msg='flux')
        self.Destroy()

    def on_reset(self, evt):

        self.tc_inpth.SetValue('3.time')
        self.tc_error.SetValue('1.0')

    def on_go(self, evt):

        inpth = self.tc_inpth.Value
        error = float(self.tc_error.Value)

        def func(argv):
            self.drp.do_get_flux_factor(argv)
            self.drp.do_abs_flux_corr(argv)

        cl = '-inpth %s -err %f' %(inpth, error)
        Misc.NewThread(func, cl).start()


class View(wx.Dialog):

    def __init__(self, parent, *args, **kwds):

        kwds['style'] = wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        wx.Dialog.__init__(self, parent, *args, **kwds)
        self.cfg = parent.cfg
        self.drp = parent.drp
        self.GFS = parent.GFS
        self.src_dist = sys.path[0] + '/'

        self.c_pg = wx.Choice(self, -1, choices=[])
        self.sz_pg_staticbox = wx.StaticBox(self, -1, '')
        self.bt_prical = wx.Button(self, -1, 'Pri. Cal.', style=wx.BU_EXACTFIT)
        self.lb_prical = wx.ListBox(self, -1, choices=[], style=wx.LB_MULTIPLE)
        self.bt_g2p = wx.BitmapButton(self, -1,
                wx.Bitmap(self.src_dist+'arrow_u.gif', wx.BITMAP_TYPE_ANY))
        self.bt_p2g = wx.BitmapButton(self, -1,
                wx.Bitmap(self.src_dist+'arrow_d.gif', wx.BITMAP_TYPE_ANY))
        self.bt_gaincal = wx.Button(self, -1, 'Gain Cal.', style=wx.BU_EXACTFIT)
        self.lb_gcal = wx.ListBox(self, -1, choices=[], style=wx.LB_MULTIPLE)
        self.bt_p2tt = wx.BitmapButton(self, -1,
                wx.Bitmap(self.src_dist+'arrow_l.gif', wx.BITMAP_TYPE_ANY))
        self.bt_tt2p = wx.BitmapButton(self, -1,
                wx.Bitmap(self.src_dist+'arrow_r.gif', wx.BITMAP_TYPE_ANY))
        self.bt_tt2g = wx.BitmapButton(self, -1,
                wx.Bitmap(self.src_dist+'arrow_l.gif', wx.BITMAP_TYPE_ANY))
        self.bt_g2tt = wx.BitmapButton(self, -1,
                wx.Bitmap(self.src_dist+'arrow_r.gif', wx.BITMAP_TYPE_ANY))
        self.bt_tt2t = wx.BitmapButton(self, -1,
                wx.Bitmap(self.src_dist+'arrow_l.gif', wx.BITMAP_TYPE_ANY))
        self.bt_t2tt = wx.BitmapButton(self, -1,
                wx.Bitmap(self.src_dist+'arrow_r.gif', wx.BITMAP_TYPE_ANY))
        self.bt_targets = wx.Button(self, -1, 'Targets', style=wx.BU_EXACTFIT)
        self.lb_targets = wx.ListBox(self, -1, choices=[], style=wx.LB_MULTIPLE)
        self.bt_t2g = wx.BitmapButton(self, -1,
                wx.Bitmap(self.src_dist+'arrow_u.gif', wx.BITMAP_TYPE_ANY))
        self.bt_g2t = wx.BitmapButton(self, -1,
                wx.Bitmap(self.src_dist+'arrow_d.gif', wx.BITMAP_TYPE_ANY))
        self.bt_timecal = wx.Button(self, -1, 'Time Cal.', style=wx.BU_EXACTFIT)
        self.lb_tcal = wx.ListBox(self, -1, choices=[], style=wx.LB_MULTIPLE)
        self.sz_src_staticbox = wx.StaticBox(self, -1, '')
        self.c_line = wx.Choice(self, -1, choices=['0', '1'])
        self.c_show = wx.Choice(self, -1, choices=['0', '1'])
        self.c_proc = wx.Choice(self, -1, choices=['0.point'])
        self.c_xaixs = wx.Choice(self, -1,
                choices=['mjd', 'scan', 'azi', 'elv', 'utc', 'lst', 'sunang'])
        self.c_yaixs = wx.Choice(self, -1,
            choices=['amp', 'azi', 'elv', 'cal', 'tsys', 'pa', 'err'])
        self.c_yerr = wx.Choice(self, -1, choices=['err', 'none'])
        self.sz_cfg_staticbox = wx.StaticBox(self, -1, '')
        self.lb_scans = wx.ListBox(self, -1, choices=[])
        self.sz_scans_staticbox = wx.StaticBox(self, -1, '')
        self.sz_bt_staticbox = wx.StaticBox(self, -1, '')

        self.Bind(wx.EVT_CLOSE, self.on_cancel)

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_BUTTON, self.on_prical, self.bt_prical)
        self.Bind(wx.EVT_BUTTON, self.on_gcal, self.bt_gaincal)
        self.Bind(wx.EVT_BUTTON, self.on_tcal, self.bt_timecal)
        self.Bind(wx.EVT_BUTTON, self.on_targets, self.bt_targets)
        self.Bind(wx.EVT_BUTTON, self.on_g2p, self.bt_g2p)
        self.Bind(wx.EVT_BUTTON, self.on_p2g, self.bt_p2g)
        self.Bind(wx.EVT_BUTTON, self.on_tt2p, self.bt_p2tt)
        self.Bind(wx.EVT_BUTTON, self.on_p2tt, self.bt_tt2p)
        self.Bind(wx.EVT_BUTTON, self.on_g2t, self.bt_g2t)
        self.Bind(wx.EVT_BUTTON, self.on_t2g, self.bt_t2g)
        self.Bind(wx.EVT_BUTTON, self.on_g2tt, self.bt_g2tt)
        self.Bind(wx.EVT_BUTTON, self.on_tt2g, self.bt_tt2g)
        self.Bind(wx.EVT_BUTTON, self.on_t2tt, self.bt_t2tt)
        self.Bind(wx.EVT_BUTTON, self.on_tt2t, self.bt_tt2t)
        self.Bind(wx.EVT_LISTBOX_DCLICK, self.dclick_scan, self.lb_scans)

        self.chk_pid()
        self.slct_srcs = []

    def __set_properties(self):
        self.SetTitle('View')
#        self.SetSize((self.GFS*600, 210*(self.GFS-1)+800))
        self.SetMinSize((self.GFS*540, 210*(self.GFS-1)+960))
#        self.SetMaxSize((self.GFS*600, -1))

        self.c_line.SetSelection(0)
        self.c_show.SetSelection(0)
        self.c_proc.SetSelection(0)
        self.c_xaixs.SetSelection(0)
        self.c_yaixs.SetSelection(0)
        self.c_yerr.SetSelection(0)


    def __do_layout(self):
        sz0 = wx.BoxSizer(wx.HORIZONTAL)
        sz2 = wx.BoxSizer(wx.VERTICAL)
        self.sz_bt_staticbox.Lower()
        sz_bt = wx.StaticBoxSizer(self.sz_bt_staticbox, wx.HORIZONTAL)
        self.sz_scans_staticbox.Lower()
        sz_scans = wx.StaticBoxSizer(self.sz_scans_staticbox, wx.VERTICAL)
        sz1 = wx.BoxSizer(wx.VERTICAL)
        self.sz_cfg_staticbox.Lower()
        sz_cfg = wx.StaticBoxSizer(self.sz_cfg_staticbox, wx.HORIZONTAL)
        sz_go = wx.BoxSizer(wx.VERTICAL)
        sz_cfg2 = wx.BoxSizer(wx.VERTICAL)
        sz_cfg1 = wx.BoxSizer(wx.VERTICAL)
        self.sz_src_staticbox.Lower()
        sz_src = wx.StaticBoxSizer(self.sz_src_staticbox, wx.HORIZONTAL)
        sz_src2 = wx.BoxSizer(wx.VERTICAL)
        sz11_copy = wx.BoxSizer(wx.HORIZONTAL)
        sz_src1 = wx.BoxSizer(wx.VERTICAL)
        sz_src11 = wx.BoxSizer(wx.VERTICAL)
        sz_src12 = wx.BoxSizer(wx.VERTICAL)
        sz_src13 = wx.BoxSizer(wx.VERTICAL)
        sz_src0 = wx.BoxSizer(wx.VERTICAL)
        sz11 = wx.BoxSizer(wx.HORIZONTAL)
        self.sz_pg_staticbox.Lower()
        sz_pg = wx.StaticBoxSizer(self.sz_pg_staticbox, wx.HORIZONTAL)
        gsz_pg = wx.GridSizer(3, 2, 1, 3)
        gsz_pg.Add(Misc.wx_st(self, 'Select PGWIN'), 1,
                wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 0)
        gsz_pg.Add(self.c_pg, 0, wx.EXPAND, 0)
        gsz_pg.Add(Misc.wx_bt(self, 'Open XW', self.on_open_xw,
            style=wx.BU_EXACTFIT), 1, wx.EXPAND, 0)
        gsz_pg.Add(Misc.wx_bt(self, 'Open CPS', self.on_open_cps,
            style=wx.BU_EXACTFIT), 1, wx.EXPAND, 0)
        gsz_pg.Add(Misc.wx_bt(self, 'Close Current', self.on_cls_crt,
            style=wx.BU_EXACTFIT), 1, wx.BOTTOM | wx.EXPAND, 5)
        gsz_pg.Add(Misc.wx_bt(self, 'Close All', self.on_cls_all,
            style=wx.BU_EXACTFIT), 1, wx.BOTTOM | wx.EXPAND, 5)
        sz_pg.Add(gsz_pg, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz1.Add(sz_pg, 0, wx.LEFT | wx.EXPAND, 5)
        sz_src0.Add(self.bt_prical, 0, wx.BOTTOM | wx.EXPAND, 3)
        sz_src0.Add(self.lb_prical, 1, wx.EXPAND, 0)
        sz11.Add(wx.Panel(self, -1), 1, wx.EXPAND, 0)
        sz11.Add(self.bt_g2p, 0, wx.ALL | wx.EXPAND, 2)
        sz11.Add(self.bt_p2g, 0, wx.ALL, 2)
        sz11.Add(wx.Panel(self, -1), 1, wx.EXPAND, 0)
        sz_src0.Add(sz11, 0, wx.TOP | wx.BOTTOM | wx.EXPAND, 1)
        sz_src0.Add(self.bt_gaincal, 0, wx.BOTTOM | wx.EXPAND, 3)
        sz_src0.Add(self.lb_gcal, 1, wx.EXPAND, 0)
        sz11_copy.Add(wx.Panel(self, -1), 1, wx.EXPAND, 0)
        sz11_copy.Add(self.bt_t2g, 0, wx.ALL | wx.EXPAND, 2)
        sz11_copy.Add(self.bt_g2t, 0, wx.ALL, 2)
        sz11_copy.Add(wx.Panel(self, -1), 1, wx.EXPAND, 0)
        sz_src0.Add(sz11_copy, 0, wx.TOP | wx.BOTTOM | wx.EXPAND, 1)
        sz_src0.Add(self.bt_timecal, 0, wx.BOTTOM | wx.EXPAND, 3)
        sz_src0.Add(self.lb_tcal, 1, wx.EXPAND, 0)
        sz_src.Add(sz_src0, 1, wx.LEFT | wx.BOTTOM | wx.EXPAND, 5)
        sz_src1.Add(sz_src11, 1, wx.EXPAND, 0)
        sz_src1.Add(sz_src12, 1, wx.EXPAND, 0)
        sz_src1.Add(sz_src13, 1, wx.EXPAND, 0)
        sz_src11.Add(wx.Panel(self, -1), 1, wx.EXPAND, 0)
        sz_src11.Add(self.bt_p2tt, 0, wx.LEFT | wx.TOP | wx.BOTTOM, 2)
        sz_src11.Add(self.bt_tt2p, 0, wx.LEFT | wx.TOP | wx.BOTTOM, 2)
        sz_src11.Add(wx.Panel(self, -1), 1, wx.EXPAND, 0)
#
        sz_src12.Add(wx.Panel(self, -1), 1, wx.EXPAND, 0)
        sz_src12.Add(self.bt_tt2g, 0, wx.LEFT | wx.TOP | wx.BOTTOM, 2)
        sz_src12.Add(self.bt_g2tt, 0, wx.LEFT | wx.TOP | wx.BOTTOM, 2)
        sz_src12.Add(wx.Panel(self, -1), 1, wx.EXPAND, 0)
#
        sz_src13.Add(wx.Panel(self, -1), 1, wx.EXPAND, 0)
        sz_src13.Add(self.bt_tt2t, 0, wx.LEFT | wx.TOP | wx.BOTTOM, 2)
        sz_src13.Add(self.bt_t2tt, 0, wx.LEFT | wx.TOP | wx.BOTTOM, 2)
        sz_src13.Add(wx.Panel(self, -1), 1, wx.EXPAND, 0)
        sz_src.Add(sz_src1, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 3)
        sz_src2.Add(self.bt_targets, 0, wx.BOTTOM | wx.EXPAND, 3)
        sz_src2.Add(self.lb_targets, 1, wx.EXPAND, 0)
        sz_src.Add(sz_src2, 1, wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz1.Add(sz_src, 1, wx.LEFT | wx.EXPAND, 5)
        sz_cfg1.Add(Misc.wx_st(self, 'Line'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg1.Add(Misc.wx_st(self, 'Show'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg1.Add(Misc.wx_st(self, 'Proc'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg1.Add(Misc.wx_st(self, 'X-aixs'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg1.Add(Misc.wx_st(self, 'Y-aixs'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg1.Add(Misc.wx_st(self, 'Y-Err'), 1, wx.ALIGN_RIGHT, 0)
        sz_cfg.Add(sz_cfg1, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 5)
        sz_cfg2.Add(self.c_line, 0, wx.EXPAND, 0)
        sz_cfg2.Add(self.c_show, 0, wx.EXPAND, 0)
        sz_cfg2.Add(self.c_proc, 0, wx.EXPAND, 0)
        sz_cfg2.Add(self.c_xaixs, 0, wx.EXPAND, 0)
        sz_cfg2.Add(self.c_yaixs, 0, wx.EXPAND, 0)
        sz_cfg2.Add(self.c_yerr, 0, wx.EXPAND, 0)
        sz_cfg.Add(sz_cfg2, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz_go.Add(Misc.wx_bt(self, 'Go', self.on_go), 1, wx.EXPAND, 0)
        sz_cfg.Add(sz_go, 1, wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz1.Add(sz_cfg, 0, wx.LEFT | wx.BOTTOM | wx.EXPAND, 5)
        sz0.Add(sz1, 1, wx.EXPAND, 0)
        sz_scans.Add(Misc.wx_st(self, 'List of Scans'), 0,
                wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 5)
        sz_scans.Add(self.lb_scans, 1, wx.EXPAND, 0)
        sz2.Add(sz_scans, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Previous', self.on_prev), 1,
                wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz_bt.Add(Misc.wx_bt(self, 'Next', self.on_next), 1,
                wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz2.Add(sz_bt, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz0.Add(sz2, 0, wx.EXPAND, 0)
        self.SetSizer(sz0)
        sz0.Fit(self)
        self.Layout()

        self.load_srcs()
        self.load_scans()
        self.load_procs()

    def load_srcs(self):

        self.lb_prical.SetItems(self.cfg.pri_cal)
        self.lb_gcal.SetItems(self.cfg.gain_cal)
        self.lb_tcal.SetItems(self.cfg.time_cal)
        self.lb_targets.SetItems(self.cfg.target)

    def load_scans(self):

        self.clf = Calibration.ListFlux(self.cfg, '0.point', '0.point')
        self.clf.load_listflux()
        self.scans = self.clf.data[2]
        self.srcs = self.clf.src
        cont = ['%04d %12s' %(a, b) for (a,b) in zip(self.scans, self.srcs)]
        self.lb_scans.SetItems(cont)

    def load_procs(self):

        procs = os.listdir('./')
        procs = [e for e in procs if (e[0].isdigit() and os.path.isdir(e))]
        procs.sort()
        self.c_proc.SetItems(procs)
        self.c_proc.Select(0)

    def chk_pid(self):

        if self.drp.all_pid:
            self.c_pg.SetItems([str(e) for e in self.drp.all_pid])
            pid = pgqid()
            idx = self.drp.all_pid.index(pid)
            self.c_pg.Select(idx)

        else:
            self.c_pg.SetItems([])


    def on_open_xw(self, evt):

        self.drp.do_openX('xw')

        self.c_pg.SetItems([str(e) for e in self.drp.all_pid])
        self.c_pg.Select(len(self.drp.all_pid)-1)

    def on_open_cps(self, evt):

        self.drp.do_openX('cps')

        self.c_pg.SetItems([str(e) for e in self.drp.all_pid])
        self.c_pg.Select(len(self.drp.all_pid)-1)

    def on_cls_crt(self, evt):

        idx = self.c_pg.Selection
        pid = self.drp.all_pid[idx]
        pgslct(pid)
        pgask(0)
        pgclos(pid)

        self.drp.all_pid.remove(pid)
        if pid == self.drp.gs_pid: self.drp.gs_pid = 0

        if self.drp.all_pid:
            new_pid = self.drp.all_pid[-1]
            pgslct(new_pid)

        self.c_pg.SetItems([str(e) for e in self.drp.all_pid])
        self.c_pg.Select(len(self.drp.all_pid)-1)

    def on_cls_all(self, evt):

        n = len(self.drp.all_pid)
        for i in range(n):
            self.on_cls_crt(evt)


    def on_prical(self, evt):

        if self.bt_prical.BackgroundColour == self.BackgroundColour:
            self.bt_prical.SetBackgroundColour((0, 255, 0))
            for i in range(len(self.cfg.pri_cal)):
                self.lb_prical.Select(i)
        else:
            self.bt_prical.SetBackgroundColour(self.BackgroundColour)
            self.slct_srcs = [e for e in self.slct_srcs if \
                                e not in self.cfg.pri_cal]
            self.slct_srcs.sort()
            self.lb_prical.Select(-1)


    def on_gcal(self, evt):

        if self.bt_gaincal.BackgroundColour == self.BackgroundColour:
            self.bt_gaincal.SetBackgroundColour((0, 255, 0))
            for i in range(len(self.cfg.gain_cal)):
                self.lb_gcal.Select(i)
        else:
            self.bt_gaincal.SetBackgroundColour(self.BackgroundColour)
            self.slct_srcs = [e for e in self.slct_srcs if \
                                e not in self.cfg.gain_cal]
            self.slct_srcs.sort()
            self.lb_gcal.Select(-1)

    def on_tcal(self, evt):

        if self.bt_timecal.BackgroundColour == self.BackgroundColour:
            self.bt_timecal.SetBackgroundColour((0, 255, 0))
            for i in range(len(self.cfg.time_cal)):
                self.lb_tcal.Select(i)
        else:
            self.bt_timecal.SetBackgroundColour(self.BackgroundColour)
            self.slct_srcs = [e for e in self.slct_srcs if \
                                e not in self.cfg.time_cal]
            self.slct_srcs.sort()
            self.lb_tcal.Select(-1)

    def on_targets(self, evt):

        if self.bt_targets.BackgroundColour == self.BackgroundColour:
            self.bt_targets.SetBackgroundColour((0, 255, 0))
            for i in range(len(self.cfg.target)):
                self.lb_targets.Select(i)
        else:
            self.bt_targets.SetBackgroundColour(self.BackgroundColour)
            self.slct_srcs = [e for e in self.slct_srcs if \
                                e not in self.cfg.target]
            self.slct_srcs.sort()
            self.lb_targets.Select(-1)

    def on_g2p(self, evt):

        idxs = self.lb_gcal.Selections
        srcs = [self.cfg.gain_cal[e] for e in idxs]
        self.cfg.pri_cal += srcs
        self.cfg.pri_cal = list(set(self.cfg.pri_cal))
        self.cfg.pri_cal.sort()

        self.lb_prical.SetItems(self.cfg.pri_cal)
        self.cfg.write_sources()


    def on_p2g(self, evt):

        idxs = self.lb_prical.Selections
        srcs = [self.cfg.pri_cal[e] for e in idxs]
        self.cfg.gain_cal += srcs
        self.cfg.gain_cal = list(set(self.cfg.gain_cal))
        self.cfg.gain_cal.sort()

        self.lb_gcal.SetItems(self.cfg.gain_cal)
        self.cfg.write_sources()


    def on_t2g(self, evt):

        idxs = self.lb_tcal.Selections
        srcs = [self.cfg.time_cal[e] for e in idxs]
        self.cfg.gain_cal += srcs
        self.cfg.gain_cal = list(set(self.cfg.gain_cal))
        self.cfg.gain_cal.sort()

        self.lb_gcal.SetItems(self.cfg.gain_cal)
        self.cfg.write_sources()

    def on_g2t(self, evt):

        idxs = self.lb_gcal.Selections
        srcs = [self.cfg.gain_cal[e] for e in idxs]
        self.cfg.time_cal += srcs
        self.cfg.time_cal = list(set(self.cfg.time_cal))
        self.cfg.time_cal.sort()

        self.lb_tcal.SetItems(self.cfg.time_cal)
        self.cfg.write_sources()

    def on_tt2p(self, evt):

        idxs = self.lb_targets.Selections
        srcs = [self.cfg.target[e] for e in idxs]
        self.cfg.pri_cal += srcs
        self.cfg.pri_cal = list(set(self.cfg.pri_cal))
        self.cfg.pri_cal.sort()
        self.cfg.target = [e for e in self.cfg.target if e not in srcs]
        self.cfg.target.sort()

        self.lb_prical.SetItems(self.cfg.pri_cal)
        self.lb_targets.SetItems(self.cfg.target)
        self.cfg.write_sources()


    def on_p2tt(self, evt):

        idxs = self.lb_prical.Selections
        srcs = [self.cfg.pri_cal[e] for e in idxs]
        self.cfg.target += srcs
        self.cfg.target = list(set(self.cfg.target))
        self.cfg.target.sort()
        self.cfg.pri_cal = [e for e in self.cfg.pri_cal if e not in srcs]
        self.cfg.pri_cal.sort()

        self.lb_prical.SetItems(self.cfg.pri_cal)
        self.lb_targets.SetItems(self.cfg.target)
        self.cfg.write_sources()

    def on_t2tt(self, evt):

        idxs = self.lb_tcal.Selections
        srcs = [self.cfg.time_cal[e] for e in idxs]
        self.cfg.target += srcs
        self.cfg.target = list(set(self.cfg.target))
        self.cfg.target.sort()
        self.cfg.time_cal = [e for e in self.cfg.time_cal if e not in srcs]
        self.cfg.time_cal.sort()

        self.lb_tcal.SetItems(self.cfg.time_cal)
        self.lb_targets.SetItems(self.cfg.target)
        self.cfg.write_sources()


    def on_tt2t(self, evt):

        idxs = self.lb_targets.Selections
        srcs = [self.cfg.target[e] for e in idxs]
        self.cfg.time_cal += srcs
        self.cfg.time_cal = list(set(self.cfg.time_cal))
        self.cfg.time_cal.sort()
        self.cfg.target = [e for e in self.cfg.target if e not in srcs]
        self.cfg.target.sort()

        self.lb_tcal.SetItems(self.cfg.time_cal)
        self.lb_targets.SetItems(self.cfg.target)
        self.cfg.write_sources()


    def on_g2tt(self, evt):

        idxs = self.lb_gcal.Selections
        srcs = [self.cfg.gain_cal[e] for e in idxs]
        self.cfg.target += srcs
        self.cfg.target = list(set(self.cfg.target))
        self.cfg.target.sort()
        self.cfg.gain_cal = [e for e in self.cfg.gain_cal if e not in srcs]
        self.cfg.gain_cal.sort()

        self.lb_gcal.SetItems(self.cfg.gain_cal)
        self.lb_targets.SetItems(self.cfg.target)
        self.cfg.write_sources()


    def on_tt2g(self, evt):

        idxs = self.lb_targets.Selections
        srcs = [self.cfg.target[e] for e in idxs]
        self.cfg.gain_cal += srcs
        self.cfg.gain_cal = list(set(self.cfg.gain_cal))
        self.cfg.gain_cal.sort()
        self.cfg.target = [e for e in self.cfg.target if e not in srcs]
        self.cfg.target.sort()

        self.lb_gcal.SetItems(self.cfg.gain_cal)
        self.lb_targets.SetItems(self.cfg.target)
        self.cfg.write_sources()


    def dclick_scan(self, evt):

        pid = pgqid()
        if not self.drp.gs_pid:
            self.drp.gs_pid = pid

        idx = self.lb_scans.Selection
        scannum = self.scans[idx]
        self.drp.do_chk_fit(str(scannum))


    def on_prev(self, evt):

        pid = pgqid()
        if not self.drp.gs_pid:
            self.drp.gs_pid  = pid

        pgslct(self.drp.gs_pid)

        idx = self.lb_scans.Selection
        idx -= 1
        self.lb_scans.Select(idx)
        scannum = self.scans[idx]
        self.drp.do_chk_fit(str(scannum))


    def on_next(self, evt):

        pid = pgqid()
        if not self.drp.gs_pid:
            self.drp.gs_pid  = pid

        pgslct(self.drp.gs_pid)

        idx = self.lb_scans.Selection
        idx += 1
        self.lb_scans.Select(idx)
        scannum = self.scans[idx]
        self.drp.do_chk_fit(str(scannum))

    def on_go(self, evt):

        idx = self.c_pg.Selection
        pid = self.drp.all_pid[idx]
        pgslct(pid)
        pgask(0)
        self.chk_pid()

        line = int(self.c_line.Selection)
        show = int(self.c_show.Selection)
        proc = self.c_proc.Items[self.c_proc.Selection]
        self.load_procs()
        idx = self.c_proc.Items.index(proc)
        self.c_proc.Select(idx)

        x = self.c_xaixs.Items[self.c_xaixs.Selection]
        y = self.c_yaixs.Items[self.c_yaixs.Selection]
        yerr = self.c_yerr.Items[self.c_yerr.Selection]
        if yerr == 'none':  sets = [proc, x, y]
        else: sets = [proc, x, y, yerr]

        self.slct_srcs = []
        if self.bt_prical.BackgroundColour == self.BackgroundColour:
            idx = self.lb_prical.Selections
            for i in idx:
                self.slct_srcs.append(self.lb_prical.Items[i])
        else:
            self.slct_srcs += self.cfg.pri_cal

        if self.bt_gaincal.BackgroundColour == self.BackgroundColour:
            idx = self.lb_gcal.Selections
            for i in idx:
                self.slct_srcs.append(self.lb_gcal.Items[i])

        else:
            self.slct_srcs += self.cfg.gain_cal

        if self.bt_timecal.BackgroundColour == self.BackgroundColour:
            idx = self.lb_tcal.Selections
            for i in idx:
                self.slct_srcs.append(self.lb_tcal.Items[i])
        else:
            self.slct_srcs += self.cfg.time_cal

        if self.bt_targets.BackgroundColour == self.BackgroundColour:
            idx = self.lb_targets.Selections
            for i in idx:
                self.slct_srcs.append(self.lb_targets.Items[i])
        else:
            self.slct_srcs += self.cfg.target

        self.slct_srcs = list(set(self.slct_srcs))
        self.slct_srcs.sort()

        cl = '-show %d -line %d -sets %s %s' \
                %(show, line, (' ').join(sets), (',').join(self.slct_srcs))

        Misc.NewThread(self.drp.do_view, cl).start()


    def on_cancel(self, evt):

        wx.CallAfter(pub.sendMessage, 'bt_color', msg='view')
        self.Destroy()

class Pipe(wx.Dialog):

    def __init__(self, parent, *args, **kwds):

        kwds["style"] = wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        wx.Dialog.__init__(self, parent, *args, **kwds)
        self.drp = parent.drp
        self.cfg = parent.drp.cfg
        self.GFS = self.Parent.GFS
        self.SetSize((self.GFS*760, 120*(self.GFS-1)+230))
        self.SetMinSize((self.GFS*700, 120*(self.GFS-1)+230))
        self.panel = wx.Panel(self, -1)

        self.SetTitle('Pipeline')
        self.p1 = wx.TextCtrl(self.panel, -1, "", style=wx.TE_MULTILINE)
        self.p2 = wx.TextCtrl(self.panel, -1, "", style=wx.TE_MULTILINE)
        self.p1.SetFont(wx.Font(self.GFS*10, wx.MODERN, wx.NORMAL, wx.NORMAL, 0, ""))
        self.p2.SetFont(wx.Font(self.GFS*10, wx.MODERN, wx.NORMAL, wx.NORMAL, 0, ""))

        self.box = wx.StaticBoxSizer(wx.StaticBox(self.panel, -1, ""), wx.VERTICAL)
        self.sizerU = wx.BoxSizer(wx.HORIZONTAL)
        self.sizerB = wx.GridSizer(rows=1, cols=5, hgap=2, vgap=1)
        self.box.Add(self.sizerU, 1, wx.EXPAND|wx.TOP, 3)
        self.box.Add(self.sizerB, 0, wx.EXPAND)

        self.sizerU.Add(self.p1, 1, wx.EXPAND|wx.RIGHT, 2)
        self.sizerU.Add(self.p2, 1, wx.EXPAND|wx.LEFT, 2)
        self.createButtonBar(self.panel, self.sizerB)
        self.panel.SetSizer(self.box)

        self.LoadPipe()

        self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)

        self.Bind(wx.EVT_CLOSE, self.on_cancel)


    def createButtonBar(self, parent, sizer):

        buttonData = (("Procedure 1", self.OnP1),
                ("Reset", self.OnReset),
                ("Go", self.OnGo),
                ("Refresh", self.OnRefresh),
                ("Procedure 2", self.OnP2))

        for eachLabel, eachHandler in buttonData:
            Button = self.buildOneButton(parent, eachLabel, eachHandler)
            sizer.Add(Button, 1, wx.EXPAND|wx.TOP, 2)


    def buildOneButton(self, parent, label, handler, pos=(0,0)):

        button = wx.Button(parent, -1, label, pos, size=(self.GFS*70,self.GFS*29))
        self.Bind(wx.EVT_BUTTON, handler, button)
        return button


    def OnCloseWindow(self, evt):

        self.Destroy()


    def LoadPipe(self):

        if not os.path.isfile('pipeline'):
            import Calibration
            clf = Calibration.ListFlux(self.cfg, '0.point', '0.point')
            self.drp.creat_pipe_script(clf)
            del Calibration, clf

        with open('pipeline') as f:
            trs = f.read()[:-1].split('\n')
        trs = [e.split() for e in trs]
        trs = [(' ').join(e) for e in trs if len(e)>0]
        line1 = [e for e in trs if 'procedure 1' in e][0]
        line2 = [e for e in trs if 'procedure 2' in e][0]
        idx1 = trs.index(line1)
        idx2 = trs.index(line2)
        P1 = trs[idx1+2:idx2]
        P2 = trs[idx2+2:]
        self.defP1 = P1
        self.defP2 = P2
        P1 = ('\n').join(P1)
        P2 = ('\n').join(P2)

        self.p1.SetValue(P1)
        self.p2.SetValue(P2)

    def OnP1(self, evt):

        self.LoadPipe()
        self.defP1 = [re.sub('#', '', e) for e in self.defP1]
        P1 = ('\n').join(self.defP1)
        self.p1.SetValue(P1)

        self.defP2 = [e.startswith('#') and e or '#'+e for e in self.defP2]
        P2 = ('\n').join(self.defP2)
        self.p2.SetValue(P2)

        f = open("pipeline", "w")
        print(('# pipeline script procedure 1:'
        '\n# -----------------------------------'), file=f)
        print(self.p1.GetValue(), file=f)
        print(('# pipeline script procedure 2:'
        '\n# -----------------------------------'), file=f)
        print(self.p2.GetValue(), file=f)
        f.close()

    def OnReset(self, evt):

        if os.path.isfile('pipeline'): os.remove('pipeline')
        self.LoadPipe()
        #self.p1.SetValue(self.defP1)
        #self.p2.SetValue(self.defP2)
        self.OnRefresh(evt)

    def OnRefresh(self, evt):
        f = open("pipeline", "w")
        print(('# pipeline script procedure 1:'
        '\n# -----------------------------------'), file=f)
        print(self.p1.GetValue(), file=f)
        print(('# pipeline script procedure 2:'
        '\n# -----------------------------------'), file=f)
        print(self.p2.GetValue(), file=f)
        f.close()

    def OnP2(self, evt):

        self.LoadPipe()
        self.defP2 = [re.sub('#', '', e) for e in self.defP2]
        P2 = ('\n').join(self.defP2)
        self.p2.SetValue(P2)

        self.defP1 = [e.startswith('#') and e or '#'+e for e in self.defP1]
        P1 = ('\n').join(self.defP1)
        self.p1.SetValue(P1)

        f = open("pipeline", "w")
        print(('# pipeline script procedure 1:'
        '\n# -----------------------------------'), file=f)
        print(self.p1.GetValue(), file=f)
        print(('# pipeline script procedure 2:'
        '\n# -----------------------------------'), file=f)
        print(self.p2.GetValue(), file=f)
        f.close()

    def on_cancel(self, evt):

        wx.CallAfter(pub.sendMessage, 'bt_color', msg='pipe')
        self.Destroy()


    def OnGo(self, evt):

        f = open("pipeline", "w")
        print(('# pipeline script procedure 1:'
        '\n# -----------------------------------'), file=f)
        print(self.p1.GetValue(), file=f)
        print(('# pipeline script procedure 2:'
        '\n# -----------------------------------'), file=f)
        print(self.p2.GetValue(), file=f)
        f.close()

        Misc.NewThread(self.drp.do_pipe, 0).start()


class Option(wx.Dialog):
    def __init__(self, parent, *args, **kwds):

        kwds['style'] = wx.DEFAULT_DIALOG_STYLE
        wx.Dialog.__init__(self, parent, *args, **kwds)
        self.cfg = parent.drp.cfg
        self.GFS = parent.GFS
        self.c_tele = wx.Choice(self, -1,
                choices=['NSRT', 'Effelsberg', 'Customize'])
        self.tc_dia = wx.TextCtrl(self, -1, '')
        self.tc_lon = wx.TextCtrl(self, -1, '')
        self.tc_lat = wx.TextCtrl(self, -1, '')
        self.tc_freq = wx.TextCtrl(self, -1, '')
        self.tc_tcal = wx.TextCtrl(self, -1, '')
        self.sz1_staticbox = wx.StaticBox(self, -1, '')

        self.Bind(wx.EVT_CLOSE, self.on_cancel)

        self.__set_properties()
        self.__do_layout()

        self.load_cfg()

    def __set_properties(self):
        self.SetTitle('Option')
        self.c_tele.SetSelection(0)

    def __do_layout(self):

        sz1 = wx.BoxSizer(wx.VERTICAL)
        sz2 = Misc.wx_hsbs(self)
        sz3 = Misc.wx_hsbs(self)
        sz5_copy_1 = wx.BoxSizer(wx.VERTICAL)
        sz5_copy = wx.BoxSizer(wx.VERTICAL)
        sz5 = wx.BoxSizer(wx.VERTICAL)
        sz5.Add(Misc.wx_st(self, 'Telescope'), 1,
                wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        sz5.Add(Misc.wx_st(self, 'Diameter'), 1,
                wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        sz5.Add(Misc.wx_st(self, 'Obs_lon'), 1,
                wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        sz5.Add(Misc.wx_st(self, 'Obs_lat'), 1,
                wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        sz5.Add(Misc.wx_st(self, 'Frequency'), 1,
                wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        sz5.Add(Misc.wx_st(self, 'Tcal'), 1,
                wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        sz3.Add(sz5, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 10)
        sz5_copy.Add(self.c_tele, 0, wx.EXPAND, 0)
        sz5_copy.Add(self.tc_dia, 0, wx.EXPAND, 0)
        sz5_copy.Add(self.tc_lon, 0, wx.EXPAND, 0)
        sz5_copy.Add(self.tc_lat, 0, wx.EXPAND, 0)
        sz5_copy.Add(self.tc_freq, 0, wx.EXPAND, 0)
        sz5_copy.Add(self.tc_tcal, 0, wx.BOTTOM | wx.EXPAND, 5)
        sz3.Add(sz5_copy, 1, wx.EXPAND, 0)
        sz5_copy_1.Add(wx.Panel(self), 1, wx.EXPAND, 0)
        sz5_copy_1.Add(Misc.wx_st(self, 'meter'), 1,
                wx.ALIGN_CENTER_VERTICAL, 0)
        sz5_copy_1.Add(Misc.wx_st(self, 'degree'), 1, wx.ALIGN_CENTER_VERTICAL, 0)
        sz5_copy_1.Add(Misc.wx_st(self, 'degree'), 1, wx.ALIGN_CENTER_VERTICAL, 0)
        sz5_copy_1.Add(Misc.wx_st(self, 'MHz'), 1, wx.ALIGN_CENTER_VERTICAL, 0)
        sz5_copy_1.Add(Misc.wx_st(self, 'K'), 1, wx.ALIGN_CENTER_VERTICAL, 0)
        sz3.Add(sz5_copy_1, 0, wx.LEFT | wx.RIGHT | wx.TOP | wx.EXPAND, 10)
        sz1.Add(sz3, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        sz2.Add(Misc.wx_bt(self, 'Cancel', self.on_cancel), 1,
                wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz2.Add(Misc.wx_bt(self, 'Reset', self.on_reset), 1,
                wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz2.Add(Misc.wx_bt(self, 'OK', self.on_ok), 1,
                wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)
        sz1.Add(sz2, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 5)
        self.SetSizer(sz1)
        sz1.Fit(self)
        self.Layout()

    def load_cfg(self):

        if not self.cfg.load_telescope():
            self.cfg.write_telescope()

        self.c_tele.Select(0)
        self.tc_dia.SetValue('%0.1f' %self.cfg.diam)
        self.tc_lon.SetValue('%f' %self.cfg.obs_lon)
        self.tc_lat.SetValue('%f' %self.cfg.obs_lat)
        self.tc_freq.SetValue('%0.1f' %self.cfg.freq)
        self.tc_tcal.SetValue('%0.2f' %self.cfg.tcal)

    def on_cancel(self, evt):

        wx.CallAfter(pub.sendMessage, 'bt_color', msg='option')
        self.Destroy()

    def on_reset(self, evt):
        self.cfg.tele_name = 'NSRT'
        self.cfg.diam = 26.0
        self.cfg.obs_lon = 87.178108
        self.cfg.obs_lat = 43.4709389
        self.cfg.freq = 4800.0
        self.cfg.tcal = 1.7

        self.load_cfg()

    def on_ok(self, evt):

        idx = self.c_tele.Selection
        self.cfg.tele_name = self.c_tele.Items[idx]
        self.cfg.diam = float(self.tc_dia.Value)
        self.cfg.obs_lon = float(self.tc_lon.Value)
        self.cfg.obs_lat = float(self.tc_lat.Value)
        self.cfg.freq = float(self.tc_freq.Value)
        self.cfg.tcal = float(self.tc_tcal.Value)

        self.cfg.write_telescope()
        self.Destroy()
