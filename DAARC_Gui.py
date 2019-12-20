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


import wx
import sys
import threading
from pubsub import pub
import DAARC
import Misc
import Gui_Toolkit as gtk
from perform import sys_pfm

class RedirectText(object):

    def __init__(self, aWxTextCtrl):
        self.out = aWxTextCtrl

    def write(self, string):
        if threading.currentThread().getName() is 'MainThread':
            self.out.WriteText(string)
            self.out.Update()
        else:
            wx.CallAfter(self.out.WriteText, string)
            wx.CallAfter(self.out.Update)


class MsgWin(wx.Frame):

    def __init__(self, parent):
        wx.Frame.__init__(self, parent, -1, "Message Window", size=(800,250),
        style = wx.DEFAULT_FRAME_STYLE^wx.CLOSE_BOX)
        self.panel = wx.Panel(self, -1)
        #self.GFS = wx.SystemSettings_GetFont(wx.SYS_OEM_FIXED_FONT).GetPointSize()/10.0
        self.GFS = 1.0

        self.msg = wx.TextCtrl(self.panel, -1, "",
        style=wx.TE_MULTILINE|wx.HSCROLL|wx.TE_READONLY)
        self.msg.SetBackgroundColour(wx.Colour(0, 0, 0))
        self.msg.SetForegroundColour(wx.Colour(255, 255, 255))
        self.msg.SetFont(wx.Font(self.GFS*10, wx.MODERN, wx.NORMAL, wx.NORMAL, 0, ""))

        self.sizer = wx.BoxSizer(wx.VERTICAL)

        self.sizer.Add(self.msg, 1, wx.EXPAND|wx.ALL, 5)
        self.panel.SetSizer(self.sizer)


class MainWin(wx.Frame):

    def __init__(self, *args, **kwds):
        kwds['style'] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)

        self.GFS = 1.0
        self.SetMinSize((self.GFS*240*1.2, self.GFS*120*1.2))
        self.SetMaxSize((self.GFS*240*1.5, self.GFS*120*1.5))
        self.SetTitle("DrCont")
        sz = wx.GridSizer(rows=4, cols=3, hgap=2, vgap=1) # creat a grid sizer
        self.bt_fitting = Misc.wx_bt(self, 'Fitting', self.on_fitting)
        self.bt_point = Misc.wx_bt(self, 'Pointing', self.on_pointing)
        self.bt_tau = Misc.wx_bt(self, 'Tau', self.on_tau)
        self.bt_gain = Misc.wx_bt(self, 'Gain Elv.', self.on_gain)
        self.bt_time = Misc.wx_bt(self, 'Gain Time', self.on_time)
        self.bt_flux = Misc.wx_bt(self, 'Abs. Flux', self.on_flux)
        self.bt_view = Misc.wx_bt(self, 'View', self.on_view)
        self.bt_pipe = Misc.wx_bt(self, 'Pipeline', self.on_pipe)
        self.bt_fdvar = Misc.wx_bt(self, 'Find Var.', self.on_fdvar)
        self.bt_option = Misc.wx_bt(self, 'Option', self.on_option)
        self.bt_perform = Misc.wx_bt(self, 'Performance', self.on_perform)
        self.bt_help = Misc.wx_bt(self, 'Help', self.on_help)
        sz.Add(self.bt_fitting, 1, wx.EXPAND)
        sz.Add(self.bt_point, 1, wx.EXPAND)
        sz.Add(self.bt_tau, 1, wx.EXPAND)
        sz.Add(self.bt_gain, 1, wx.EXPAND)
        sz.Add(self.bt_time, 1, wx.EXPAND)
        sz.Add(self.bt_flux, 1, wx.EXPAND)
        sz.Add(self.bt_view, 1, wx.EXPAND)
        sz.Add(self.bt_pipe, 1, wx.EXPAND)
        sz.Add(self.bt_fdvar, 1, wx.EXPAND)
        sz.Add(self.bt_option, 1, wx.EXPAND)
        sz.Add(self.bt_perform, 1, wx.EXPAND)
        sz.Add(self.bt_help, 1, wx.EXPAND)

        self.SetSizer(sz)

        self.drp = DAARC.DRP()
        self.cfg = self.drp.cfg

        output = MsgWin(self)
        redir = RedirectText(output.msg)
        sys.stdout = redir
        output.Show()

        pub.subscribe(self.turn_gray, 'bt_color')


    def turn_gray(self, msg):

        if msg == 'fit':
            self.bt_fitting.SetBackgroundColour(self.BackgroundColour)
        if msg == 'point':
            self.bt_point.SetBackgroundColour(self.BackgroundColour)
        if msg == 'tau':
            self.bt_tau.SetBackgroundColour(self.BackgroundColour)
        if msg == 'gain':
            self.bt_gain.SetBackgroundColour(self.BackgroundColour)
        if msg == 'time':
            self.bt_time.SetBackgroundColour(self.BackgroundColour)
        if msg == 'flux':
            self.bt_flux.SetBackgroundColour(self.BackgroundColour)
        if msg == 'view':
            self.bt_view.SetBackgroundColour(self.BackgroundColour)
        if msg == 'pipe':
            self.bt_pipe.SetBackgroundColour(self.BackgroundColour)
        if msg == 'option':
            self.bt_option.SetBackgroundColour(self.BackgroundColour)


    def on_fitting(self, evt):

        if self.bt_fitting.BackgroundColour == self.BackgroundColour:
            dlg = gtk.Fitting(self, -1)
            val = dlg.Show()
            self.bt_fitting.SetBackgroundColour((0, 255, 0))

    def on_pointing(self, evt):

        if self.bt_point.BackgroundColour == self.BackgroundColour:
            dlg = gtk.Pointing(self, -1)
            #val = dlg.ShowModal()
            #dlg.CenterOnScreen()
            #dlg.Destroy()
            val = dlg.Show()
            self.bt_point.SetBackgroundColour((0, 255, 0))

    def on_tau(self, evt):

        if self.bt_tau.BackgroundColour == self.BackgroundColour:
            dlg = gtk.Tau(self, -1)
            val = dlg.Show()
            self.bt_tau.SetBackgroundColour((0, 255, 0))

    def on_gain(self, evt):

        if self.bt_gain.BackgroundColour == self.BackgroundColour:
            dlg = gtk.GainElv(self, -1)
            val = dlg.Show()
            self.bt_gain.SetBackgroundColour((0, 255, 0))

    def on_time(self, evt):

        if self.bt_time.BackgroundColour == self.BackgroundColour:
            dlg = gtk.GainTime(self, -1)
            val = dlg.Show()
            self.bt_time.SetBackgroundColour((0, 255, 0))

    def on_flux(self, evt):

        if self.bt_flux.BackgroundColour == self.BackgroundColour:
            dlg = gtk.AbsFlux(self, -1)
            val = dlg.Show()
            self.bt_flux.SetBackgroundColour((0, 255, 0))

    def on_view(self, evt):

        if self.bt_view.BackgroundColour == self.BackgroundColour:
            dlg = gtk.View(self, -1)
            val = dlg.Show()
            self.bt_view.SetBackgroundColour((0, 255, 0))

    def on_pipe(self, evt):

        if self.bt_pipe.BackgroundColour == self.BackgroundColour:
            dlg = gtk.Pipe(self, -1)
            val = dlg.Show()
            self.bt_pipe.SetBackgroundColour((0, 255, 0))

    def on_fdvar(self, evt):

        self.drp.do_fdvar(None)

    def on_option(self, evt):

        if self.bt_option.BackgroundColour == self.BackgroundColour:
            dlg = gtk.Option(self, -1)
            val = dlg.Show()
            self.bt_option.SetBackgroundColour((0, 255, 0))

    def on_perform(self, evt):

        sys_pfm(self.cfg.diam, self.cfg.freq, self.cfg.tcal, self.cfg.Tgr)


    def on_help(self, evt):

        evt.Skip()


def main():

    class App(wx.App):

        def OnInit(self):
            MainWin(None, -1).Show()
            return True

    app = App()
    app.MainLoop()


if __name__ == "__main__":

    main()
