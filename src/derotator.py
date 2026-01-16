#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Joerg Weingrill"
__copyright__ = "Copyright 2023 Leibniz-Insitute for Astrophysics Potsdam (AIP)"
__credits__ = ["Joerg Weingrill"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Joerg Weingrill"
__email__ = "jweingrill@aip.de"
__status__ = "Development"
__date__ = "7/7/23"

import numpy as np
from skyfield import api



if __name__ == "__main__":
    ts = api.load.timescale()
    t = ts.utc(2023, 9, 13, 20)
    latitude = 28.301195
    longitude = -16.509209
    izana = api.wgs84.latlon(latitude_degrees=latitude, longitude_degrees=longitude)
    observer = izana.at(t)
    # Altair JNow:	19h 51m 56s	 08° 55' 45"
    # J2000:	19h 50m 46s	 08° 52' 02"
    # AzAlt:	 98° 15' 59"	 33° 23' 26"
    pleiades = position_of_radec(ra_hours, dec_degrees, t=t, center=earth)

    sinpa = np.sin(azimuth)*np.cos(np.radians(latitude))/cos(declination)

    pass
