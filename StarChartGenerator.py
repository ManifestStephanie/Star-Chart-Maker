#!/usr/bin/env python
# coding: utf-8

# Allsky chart code based on https://github.com/skyfielders/python-skyfield/discussions/636

import math
import re
from datetime import datetime, time

import numpy as np
import streamlit as st
from geopy.geocoders import Nominatim
from matplotlib import pyplot as plt
from matplotlib.backends.backend_agg import RendererAgg
from matplotlib.collections import LineCollection
from pytz import timezone
from skyfield.api import Star, load, wgs84, Angle
from skyfield.data import hipparcos, stellarium
from skyfield.projections import build_stereographic_projection
from timezonefinder import TimezoneFinder
import ecliptic

# Enter location
loc_text = st.text_input('Address', value="Space Needle", placeholder="345 Blueberry Ln Lexington KY")

#@st.cache_data
def get_location(loc_text=loc_text):
    geolocator = Nominatim(user_agent="Star_Chart_Generator")
    nom_location = geolocator.geocode(loc_text)
    long = nom_location.longitude
    lat = nom_location.latitude
    location = wgs84.latlon(lat, long)

    # Calculate the correct time zone
    obj_tz = TimezoneFinder()
    zone_name = obj_tz.timezone_at(lat=lat, lng=long)
    zone = timezone(zone_name)
    # st.write(f"zone: {zone}")
    return zone, zone_name, location, long, lat

zone, zone_name, location, long, lat = get_location()
st.write(f"Chart Long, Lat: {long}, {lat}")

# OpenStreetMaps requisite attribution
st.write("Powered by [OpenStreetMap](http://www.openstreetmap.org/copyright)")
#seattle = wgs84.latlon(47.61352679507131, -122.30535433025425, elevation_m=100)
#geolocator = Nominatim(user_agent="Star_Chart_Generator")


# Load timescale
ts = load.timescale()
ecliptic.ts = ts

# Enter date and time
date = st.date_input('Date')
#st.write(f"Date: {date}")
user_time = st.time_input('Time (Local)', value=datetime.combine(date, time(hour=22, minute=0)), step=3600)
date_time = datetime.combine(date, user_time)
date_time = zone.localize(date_time)
st.write(date_time)
#time_local = date_time.astimezone(zone)
#st.write(f"Local time: {time_local.strftime('%H:%M')}")
t = ts.from_datetime(date_time)


# 180 = South, 0 = North
degrees = 0.0
# Set zenith
zenith = location.at(t).from_altaz(alt_degrees=90, az_degrees=degrees) 

@st.cache_resource
def load_eph():
    return load('de421.bsp')
eph = load_eph()

earth = eph['earth']

eph_obj = ('sun', 'moon', 'mercury', 'venus', 'mars', 'jupiter barycenter', 'saturn barycenter',
           'uranus barycenter', 'neptune barycenter')
ecliptic.eph_sun = eph['sun']
radius_km = (695700, 1737.4, 2437.7, 6051.8, 3389.5, 69911, 58232, 25362, 24622)
radius_km_dict = dict(zip(eph_obj, radius_km))


# The Hipparcos mission provides our star catalog.
@st.cache_data
def load_stars():
    with load.open(hipparcos.URL) as f:
        return hipparcos.load_dataframe(f)
stars = load_stars()

# The constellation outlines come from Stellarium. We make a list
# of the stars at which each edge starts, and the star at which each edge
# ends

with load.open('constellationship.fab') as f:
    constellations = stellarium.parse_constellations(f)

edges = [edge for name, edges in constellations for edge in edges]
edges_star1 = [star1 for star1, star2 in edges]
edges_star2 = [star2 for star1, star2 in edges]


# Center the chart on the zenith
projection = build_stereographic_projection(zenith)
ecliptic.projection = projection
field_of_view_degrees = 180.0

# Allow user to select limiting magnitude
limiting_magnitude = st.slider(label="Limiting Magnitude of Stars", min_value=1.0, max_value=7.0, value=5.5, step=0.5)
#limiting_magnitude = 6.0


# Now that we have constructed our projection, compute the x and y
# coordinates that each star and the jupiter will have on the plot
loc = earth + location
ecliptic.loc = loc

star_positions = earth.at(t).observe(Star.from_dataframe(stars))
#star_alt = star_positions.apparent().altaz()[0].degrees
#star_positions = star_positions[ star_alt > 0 ]
stars['x'], stars['y'] = projection(star_positions)

x = {}
y = {}
altaz_dist = {}
dist_km = {}
apparent_diam_arcsec = {}

for obj in eph_obj:
    obs = loc.at(t).observe(eph[obj])
    x[obj], y[obj] = projection(obs)
    altaz_dist[obj] = obs.apparent().altaz()
    dist_km[obj] = altaz_dist[obj][-1].km

    apparent_diam_arcsec[obj] = Angle(radians=np.arcsin(radius_km_dict[obj] / dist_km[obj]) * 2.0).arcseconds()


# Create a True/False mask marking the stars bright enough to be
# included in our plot. And go ahead and compute how large their
# markers will be on the plot

bright_stars = (stars.magnitude <= limiting_magnitude)
magnitude = stars['magnitude'][bright_stars]
marker_size = (0.5 + limiting_magnitude - magnitude) ** 2.0


# The constellation lines will each begin at the x, y of one star and end
# at the x, y of another. We have to "rollaxis" the resulting coordinate
# array into the shape that matplotlib expects.

xy1 = stars[['x', 'y']].loc[edges_star1].values
xy2 = stars[['x', 'y']].loc[edges_star2].values
lines_xy = np.rollaxis(np.array([xy1, xy2]), 1)


# Time to build the figure!

_lock = RendererAgg.lock

with _lock:

    angle = np.pi - field_of_view_degrees / 360.0 * np.pi
    limit = np.sin(angle / (1.0 - np.cos(angle)))



    fig, ax = plt.subplots(figsize=[12, 12], facecolor='white')

    # Draw horizon as a dashed line
    # 24 points horizon

    horizon = []
    h0 = projection(location.at(t).from_altaz(alt_degrees=0, az_degrees=0.0))
    for i in range(1, 73):
        delta = 5.0
        current = i*delta
        h1 = projection(location.at(t).from_altaz(alt_degrees=0, az_degrees=current))
        horizon.append([h0, h1])
        h0 = h1

    ax.add_collection(LineCollection(horizon, colors='#00f2', linewidths=1, linestyle='dashed',
                                     zorder=0, alpha=0.5))

    # Draw the constellation lines

    ax.add_collection(LineCollection(lines_xy, colors=(0, 0, 1, .4)))

    # Draw the stars

    ax.scatter(stars['x'][bright_stars], stars['y'][bright_stars],
               s=marker_size, color='k', zorder=0)

    # Draw ecliptic
    x_ecliptic, y_ecliptic = ecliptic.make_ecliptic_path()
    ax.plot(x_ecliptic, y_ecliptic, linestyle='--')

    # Draw the solar system positions

    color = {'sun': 'yellow',
             'moon': 'gray',
             'mercury': 'brown',
             'venus': 'blue',
             'mars': 'red',
             'jupiter barycenter': 'darkorange',
            'saturn barycenter': 'gold',
            'uranus barycenter': 'teal',
            'neptune barycenter': 'blue'} 
    offset = 0.016

    for obj in eph_obj:
        if altaz_dist[obj][0].degrees > 0:
            # Make planet marker size relative to the log of the apparent diam
            size = math.log(apparent_diam_arcsec[obj])*25
            ax.scatter(x[obj], y[obj], s=size, marker='o', c=color[obj], edgecolors='black', zorder=3)
            #time = t.astimezone(zone).strftime('%m-%d')
            name = re.split("\s", obj)[0].capitalize()
            text = ax.text(x[obj] + offset, y[obj] - offset, name, color='black',
                           ha='left', va='top', fontsize=10, weight='bold')

    #for xi, yi, tstr in zip(x[obj], y[obj], eph_obj):
    #    time = t.utc_strftime('%m-%d')
    #    tstr = tstr.lstrip('0')
    #    text = ax.text(xi + offset, yi - offset, tstr, color=jupiter_color,
    #                   ha= 'left', va='top', fontsize=10, weight='bold')

    # Finally, title the plot and set some final parameters

    angle = np.pi - field_of_view_degrees / 360.0 * np.pi
    limit = np.sin(angle / (1.0 - np.cos(angle)))
    
    # Horizon circle
    horiz_circle = plt.Circle((0, 0), limit, color='white', zorder=-1)
    ax.add_patch(horiz_circle)
    ax.set_facecolor("linen")
    

    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_aspect(1.0)
    ax.set_title('Star Chart for {} {}'.format(
        t.astimezone(zone).strftime('%Y %B %d %H:%M'),
        zone
    ))

    st.pyplot(fig)


# In[ ]:


# TODO: Plot ecliptic
# Add Moon phase
# Add major star names
# Display Lat and Long
# Add cardinal directions