"""ecliptic.py
This module calculates the path of the ecliptic for plotting"""
import streamlit
#from datetime import datetime
from skyfield.api import Star, load, wgs84, Angle
from skyfield.data import hipparcos, stellarium
from skyfield.projections import build_stereographic_projection
import numpy as np


# Create ts object
# We will just pick an arbitrary year, since accuracy isn't paramount
# and its essentially the same from year to year

#@streamlit.cache_data
def make_ecliptic_path(year=2010):
#year = 2010
    t_sun = ts.utc(year, 1, range(365))

    obs = loc.at(t_sun).observe(eph_sun)
    x_ecliptic, y_ecliptic = projection(obs)

    return x_ecliptic, y_ecliptic
