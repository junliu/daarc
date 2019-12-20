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


import subprocess
import threading
import os
import wx

### FileIO functions ###
def comment_remover(line):

    if line.strip()[0]!='%' and line.strip()[0]!='#':
        return True

    else:
        return False

def type_guess(string):

    string = string.strip()

    try:
        dtype = type(eval(string))
        try: dtype(string)
        except: dtype = str
    except: dtype = str

    return dtype

def type_specify(string):

    if string == 'F':
        dtype = float
    elif string == 'D':
        dtype = int
    elif string == 'S':
        dtype = str
    else:
        print('Wrong data format, use i, f, or s')
        return

    return dtype


def readcol(fname, cols=None, fmt=None,
                    start=0, stop=0,
                    comment='#', flag=True):
    """
    readcol(fname, cols=None, fmt=None,
            start=0, stop=0,
            comment='#', flag=True):

    """

    with open(fname) as fobj:
        lines = fobj.readlines()
    n = len(lines)
    if stop == 0: stop = n+1
    lines = lines[start:stop]
    if flag:
        lines = [line for line in lines if line.strip()[0]!=comment]

    if not cols:
        tmp = lines[0].split()
        if flag: ncol = len(tmp)
        else:
            if tmp[0] != comment: ncol = len(tmp) + 1
            else: ncol = len(tmp)
        cols = list(range(ncol))
    else:
        if flag: ncol = len(cols)
        else:
            cols = [e+1 for e in cols]
            cols.insert(0, 0)
            ncol = len(cols)

    n = len(lines)

    if not fmt:
        dfmt = []
        line = lines[0].split()
        if len(line) < ncol: line.insert(0, '#')
        for i in cols:
            dfmt.append(type_guess(line[i]))
    else:
        fmt = [e for e in fmt if e!=',' and e!='%']
        fmt = [e.strip().upper() for e in fmt]
        fmt = [e for e in fmt if e!='']
        dfmt = [type_specify(s) for s in fmt]

    if ncol == len(dfmt) + 1: dfmt.insert(0, str)

    data = [list(range(n)) for _ in range(ncol)]
    if flag:
        for i in range(n):
            line = lines[i].split()
            for j in range(ncol):
                col = cols[j]
                data[j][i] = dfmt[j](line[col])
    else:
        for i in range(n):
            line = lines[i].split()
            if line[0] != comment: line.insert(0, ' ')
            for j in range(ncol):
                col = cols[j]
                data[j][i] = dfmt[j](line[col])

    return data



def writecol(fname, hdr, tail, data, fmt):

    ncol = len(data)
    try:
        nline = len(data[0])
    except:
        nline = ncol*1
        ncol = 1
        data = [data]

    dlist = []

    fobj = open(fname, 'w')
    print(hdr, file=fobj)
    for i in range(nline):
        dlist = []
        for j in range(ncol):
            dlist.append(data[j][i])
        dtuple = tuple(dlist)
        print(fmt %dtuple, file=fobj)

    if tail:
        print(tail, file=fobj)

    fobj.close()

### end FileIO functions ###

def array_flat(array, norm=0, ref_array=None):
    """
    array = [array1, array2, array3 ....]
    """
    out_array = np.array([])
    n = len(array)
    locate = np.zeros(n)
    if (type(norm) == int) and (norm == 0):
        norm = np.ones(n)
    elif (type(norm) == int) and (norm == 1):
        norm = np.array([np.mean(e) for e in array])

    elif (type(norm) == int) and (norm > 1):
        order = copy.deepcopy(norm)
        norm = np.zeros(n)

        for i in range(n):
            par = np.polyfit(ref_array[i], array[i], order)
            func = lambda x: -np.polyval(par, x)
            mini = minimize(func, [48.0])
            norm[i] = -mini.fun
            loc = mini.x
            if len(loc) > 1:
                idx = np.argmin(np.fabs(loc-48))[0][0]
                locate[i] = loc[idx]
            else: locate[i] = loc[0]

    for i in range(n):
        if array[i].dtype =='S1':
            out_array = np.append(out_array, array[i])
        else:
            out_array = np.append(out_array, array[i]/norm[i])

    return out_array, norm, locate

###

def get_terminal_size():

#    c = commands.getoutput('resize').split('\n')
#    print c
#    c0 = c[0].split('=')[-1][:-1]
#    c1 = c[1].split('=')[-1][:-1]
#
#    return int(c0), int(c1)

    rows, columns = os.popen('stty size', 'r').read().split()
    return int(rows), int(columns)


class NewThread(threading.Thread):
    def __init__(self, func, args):
        threading.Thread.__init__(self)
        self.func = func
        self.args = args

    def run(self):
        if self.args == None: self.func()
        else: self.func(self.args)


### build wx objects
def wx_st(parent, label):
    return wx.StaticText(parent, -1, label)


def wx_bt(parent, label, handler, **kwds):

    button = wx.Button(parent, -1, label, **kwds)
    parent.Bind(wx.EVT_BUTTON, handler, button)

    return button


def wx_hsbs(parent, **kwds):

    return wx.StaticBoxSizer(wx.StaticBox(parent, -1, ''), wx.HORIZONTAL)


def wx_vsbs(parent, **kwds):

    return wx.StaticBoxSizer(wx.StaticBox(parent, -1, ''), wx.VERTICAL)


def wx_lpt(line):
    # layout properties
    p = {'a' : '240',
        'at': '0',
        'v' : '8',
        'h' : '4',
        'b': '128',
        't' : '64',
        'e': '8192',
        'l' : '16',
        'r' : '32',
        'al' : '0',
        'ar' : '512',
        'ach' : '256',
        'acv' : '2048'}
    line = line.split('|')
    line = [p[e] for e in line]
    line = ('|').join(line)

    return eval(line)

### end
