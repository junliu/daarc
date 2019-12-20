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

import ephem
import numpy as np


def utc2mjd(ut):
    """
    Convert UTC day to Julian day.
    Returns float.
    """

    return utc2jd(ut) - 2400000.5


def utc2jd(ut):
    """
    Convert UTC day to Julian day.
    Returns float.
    """

    date = ephem.date(ut)
    return ephem.julian_date(date)


def jd2utc(jd):
    """
    Convert Julian day to UTC day.
    Returns UTC year, month, day.
    float.
    """

    return ephem.date(jd-2415020.0)


def mjd2utc(mjd):
    """
    Convert Modified Julian day to UTC day.
    Returns UTC year, month, day.
    float.
    """

    return jd2utc(mjd+2400000.5)


def mjd2lst(mjd, lon):
    """
    Local apparent sidereal time.
    Returns float.
    """

    obs = ephem.Observer()
    obs.lon = str(lon)
    obs.lat = 0
    obs.date = ephem.date(mjd+2400000.5-2415020.0)
    time = ephem.degrees(str(obs.sidereal_time()))
    return float(time)*180.0/ephem.pi


def ra_dec2az_el(obslon, obslat, ra, dec, ut, sun=0):
    obs = ephem.Observer()
    obs.long = obslon
    obs.lat = obslat
    obs.date = ut

    src = ephem.readdb('NAME,f|M|F7, %s, %s, 0, 2000' %(ra, dec))
    src.compute(obs)

    if sun:
        Sun = ephem.Sun()
        Sun.compute(obs)
        sep = ephem.separation([Sun.az, Sun.alt], [src.az, src.alt])
    else: sep = 0

    return src.az*180.0/ephem.pi, src.alt*180.0/ephem.pi, sep*180.0/ephem.pi


def cal_flux(srcname, frq):
    """
    Calculate the flux density of the calibrator at a given frq.
    Here we make 3c286 the default calibrator.
    Returns float.
    """
    # parameters are from Effeslberg radiocal
    cal_flux = {'3c286':[1.11578, 0.494671, -0.15199],
                '3c48': [2.74926, -0.155745, -0.10576],
                '3c123': [2.525, 0.2460, -0.1638],
                '3c147': [2.47011, 0.0469864, -0.128600],
                '3c161': [1.62941, 0.500078, -0.195202],
                '3c218': [4.729, -1.025, 0.01300],
                '3c227': [6.757, -2.801, 0.296900],
                '3c249.1': [2.537, -0.565, -0.0404],
                'vira': [4.484, -0.603, -0.0280],
                '3c295': [1.69210, 0.642636, -0.238765],
                '3c309.1': [2.617, -0.437, -0.0373],
                '3c348': [3.852, -0.361, -0.1053],
                '3c353': [3.148, -0.157, -0.0911],
                'cyga': [8.360, -1.565, 0]
    }

    if srcname.lower() == 'ngc7027':
        _frq = frq/1e3

        return 0.703958*_frq**2*(1-np.exp(-11.0752*_frq**-2.1))

    else:
        poly = cal_flux[srcname.lower()]
        _frq = np.log10(frq)
        return 10**(poly[0] + _frq*poly[1] + _frq**2*poly[2])


if __name__ == '__main__':
    print(ra_dec2az_el('87:10.67439', '43:28.26913', '01 37 41.2994', \
            '33 09 35.134', '2017/03/29', 0))

