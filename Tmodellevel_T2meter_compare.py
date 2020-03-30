#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""for dO3/dT research 2-meter temperatures were used for to calculate dO3/dT;
however, this isn't what the CTM uses (it uses 3d temperature fields and, 
specificially, the temperature in the lowest model box). The following code 
compares MERRA-2 2-meter temperatures to instantaneous temperatures at the 
lowest model box (985 hPa) at 18z for JJA 2008. 
"""
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import Normalize
import sys
sys.path.append('/Users/ghkerr/phd/GMI/')
import commensurability 
# Load MERRA-2 2-meter temperatures
t2m, t10m, u2m, u10m, v2m, v10m, ps, merra_lat, merra_lon, times_all = \
commensurability.load_MERRA2([2008])
# Find 18z 2-meter temperatures from hourly output
where18z = []
for time in times_all:
    if time.hour == 18.:
        time_where = np.where(np.array(times_all) == time)
        where18z.append(time_where[0][0])
t2m18z = t2m[where18z]        
# Load MERRA-2 air temperatures (should have the same spatial extent as the 
# above 2-meter temperatures, so don't reload the dimensional coordinates)
tair = Dataset('/Users/ghkerr/phd/meteorology/data/T_jja2008.nc', 'r')
tair = tair.variables['T']
# Select 3-hourly air temperatures at 18Z (1080 minutes since beginning of 
# day in 'time' variable). This is closest to the 1-2pm mean temperatures 
# from the overpass2 output. Select model level 72 (985 hPa)
tair18z = tair[:, 6, 1]
# Create map
m = Basemap(projection = 'merc', llcrnrlon = -126., llcrnrlat = 24., 
            urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'h', 
            area_thresh = 1000)
x, y = np.meshgrid(merra_lon, merra_lat)
x, y = m(x, y)    
fig = plt.figure()  
# Difference in means
ax1 = plt.subplot2grid((2, 1), (0, 0))
ax1.set_title('$\mu_{\mathregular{T model level}}$ - '+
              '$\mu_{\mathregular{T 2meter}}$', ha = 'left', 
              fontsize = 16, x = 0.03, y = 1.03)
vmin = -3.; vmax = 3.
cmap = plt.get_cmap('bwr', 10)
norm = Normalize(vmin = vmin, vmax = vmax)
clevs = np.linspace(vmin, vmax, 11, endpoint = True)
m.contourf(x, y, np.mean(tair18z, axis=0) - np.mean(t2m18z, axis=0), 
           clevs, cmap = cmap, norm = norm, extend = 'both')
m.drawstates(color = 'k', linewidth = 0.5)
m.drawcountries(color = 'k', linewidth = 1.0)
m.drawcoastlines(color = 'k', linewidth = 1.0)    
fill_oceans(ax1, m)  
cb = plt.colorbar()  
cb.set_label(label = '[K]', size = 12)
# Difference in variability
ax2 = plt.subplot2grid((2, 1), (1, 0))
ax2.set_title('$\sigma_{\mathregular{T model level}}$ - '+
              '$\sigma_{\mathregular{T 2meter}}$', ha = 'left', 
              fontsize = 16, x = 0.03, y = 1.03)
vmin = -1.; vmax = 1.
cmap = plt.get_cmap('bwr', 10)
norm = Normalize(vmin = vmin, vmax = vmax)
clevs = np.linspace(vmin, vmax, 11, endpoint = True)

m.contourf(x, y, np.std(tair18z, axis=0) - np.std(t2m18z, axis=0), 
           clevs, cmap = cmap, norm = norm, extend = 'both')
m.drawstates(color = 'k', linewidth = 0.5)
m.drawcountries(color = 'k', linewidth = 1.0)
m.drawcoastlines(color = 'k', linewidth = 1.0)    
fill_oceans(ax2, m)    
cb = plt.colorbar()  
cb.set_label(label = '[K$^{\mathregular{2}}$]', size = 12)
plt.savefig('/Users/ghkerr/phd/GMI/figs/'+'Tmodellevel_T2meter_compare.eps',
            dpi=300)