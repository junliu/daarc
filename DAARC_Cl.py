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


import Misc
import DAARC
from cmd import *

if __name__ == '__main__':

    termW = Misc.get_terminal_size()[1]
    strlen = int((termW + 50)/2)
    welcome  = '\n' + \
    '%s' %('*************************************************************\n'.rjust(strlen)) + \
    '%s' %('*    Data Reduction package for radio Continuum (DrCont)    *\n'.rjust(strlen)) + \
    '%s' %('*       observation performed with cross-scan mode.         *\n'.rjust(strlen)) + \
    '%s' %('*                                                           *\n'.rjust(strlen)) + \
    '%s' %('*                          Jun LIU                          *\n'.rjust(strlen)) + \
    '%s' %('*                       liuj@xao.ac.cn                      *\n'.rjust(strlen)) + \
    '%s' %('*           Xinjiang Astronomical Observatory, CAS          *\n'.rjust(strlen)) + \
    '%s' %('*                                                           *\n'.rjust(strlen)) + \
    '%s' %('*************************************************************\n'.rjust(strlen))
    app = DAARC.DRP()
    app.cmdloop(welcome)

