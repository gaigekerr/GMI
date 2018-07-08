#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAME
    aqs_utils.py
PURPOSE
    Script deals with trace gas measurements from AQS 
PROGRAMMER
    Gaige Hunter Kerr
REVISION HISTORY
    24062018 -- initial version created and function 'boxplot_tracegas_byenv' 
                added 
"""
# # # # # # # # # # # # #
# change font
import matplotlib 
prop = matplotlib.font_manager.FontProperties(fname = 
    '/Users/ghkerr/Library/Fonts/cmunbmr.ttf')
matplotlib.rcParams['font.family'] = prop.get_name()
prop = matplotlib.font_manager.FontProperties(fname = 
    '/Users/ghkerr/Library/Fonts/cmunbbx.ttf')
matplotlib.rcParams['mathtext.bf'] = prop.get_name()
# for unicode minus/negative sign implementation
matplotlib.rcParams['axes.unicode_minus'] = False
# # # # # # # # # # # # #
def boxplot_tracegas_byenv(co_all, co_u, co_su, co_r, mr2_co, no2_all, 
                           no2_u, no2_su, no2_r, mr2_no2, o3_all_castnet, 
                           o3_all, o3_u, o3_su, o3_r, mr2_o3, years, region):
    """function plots boxplots of regionally-averaged CO, NO2, and O3 
    concentrations for every summer day in measuring period. Boxplots 
    represent (1) an average of all AQS sites within the search radius of 
    CASTNet sites, (2) an average of urban AQS sites, (3) an average of 
    suburban sites, (4) an average of rural sites, and (5) GMI CTM. 
  
    Parameters
    ----------        
    co_all : numpy.ndarray
        AQS CO measurements averaged over all sites within search radius of 
        CASTNet sites, units of parts per million, [years in measuring 
        period * days in months in 'sampling_months',]  
    co_u : numpy.ndarray
        AQS CO measurements averaged over all urban sites within search radius 
        of CASTNet sites, units of parts per million, [years in measuring 
        period * days in months in 'sampling_months',]      
    co_su : numpy.ndarray
        AQS CO measurements averaged over all suburban sites within search 
        radius of CASTNet sites, units of parts per million, [years in 
        measuring period * days in months in 'sampling_months',]          
    co_r : numpy.ndarray
        AQS CO measurements averaged over all rural sites within search radius 
        of CASTNet sites, units of parts per million, [years in measuring 
        period * days in months in 'sampling_months',]          
    mr2_co : numpy.ndarray
        GMI CTM CO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        parts per million, [years in measuring period * days in months in 
        'sampling_months',]     
    no2_all : numpy.ndarray
        Same as 'co_all' but for NO2, units of parts per billion
    no2_u : numpy.ndarray
        Same as 'co_u' but for NO2, units of parts per billion    
    no2_su : numpy.ndarray
        Same as 'co_su' but for NO2, units of parts per billion    
    no2_r : numpy.ndarray
        Same as 'co_r' but for NO2, units of parts per billion    
    mr2_no2 : numpy.ndarray
        Same as 'mr2_co' but for NO2, units of parts per billion        
    o3_all_castnet : numpy.ndarray
        Regionally-averaged CASTNet O3 observations in region, units of parts 
        per billion, [years in measuring period * days in months in 
        'sampling_months',]      
    o3_all : numpy.ndarray 
        Same as 'co_all' but for O3, units of parts per million    
    o3_u : numpy.ndarray
        Same as 'co_all' but for O3, units of parts per million    
    o3_su : numpy.ndarray
        Same as 'co_all' but for O3, units of parts per million    
    o3_r : numpy.ndarray
        Same as 'co_all' but for O3, units of parts per million    
    mr2_o3 : numpy.ndarray
        Same as 'mr2_co' but for O3, units of parts per million        
    years : list
        Years of interest    
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function

    Returns
    ----------
    None    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    # dictionary of data for CO boxplots
    data = {}
    data['a'] = np.hstack(co_all)
    data['b'] = np.hstack(co_u)
    data['c'] = np.hstack(co_su)
    data['d'] = np.hstack(co_r)
    data['e'] = np.hstack(mr2_co)
    color_dict = {'AQS$_{\mathregular{ all}}$' : '#66c2a5',
                  'AQS$_{\mathregular{ urban}}$' : '#fc8d62',
                  'AQS$_{\mathregular{ suburban}}$' : '#8da0cb', 
                  'AQS$_{\mathregular{ rural}}$' : '#e78ac3',
                  'MR2' : 'w'}
    controls = ['AQS$_{\mathregular{ all}}$', 'AQS$_{\mathregular{ urban}}$',
                'AQS$_{\mathregular{ suburban}}$', 'AQS$_{\mathregular{ rural}}$',
                'MR2']
    # initialize figure, axis
    ax = plt.subplot2grid((1, 1), (0, 0))
    boxplot_dict = ax.boxplot(
            [data[x] for x in ['a', 'b', 'c', 'd', 'e']],
            positions = [1, 1.5, 2, 2.5, 3],
            labels = controls, 
            patch_artist = True,
            widths = 0.25, showfliers = True)
    for b in boxplot_dict['medians']:
        b.set_color('k')
        b.set_linewidth(2.)
    i = 0
    for b in boxplot_dict['boxes']:
        lab = ax.get_xticklabels()[i].get_text()
        b.set_facecolor(color_dict[lab]) 
        b.set_linewidth(2.)
        i += 1
    for b in boxplot_dict['fliers']:
        b.set_marker('.')
        b.set_markersize(1)
        b.set_markerfacecolor('k')
    for b in boxplot_dict['whiskers']:
        b.set_linewidth(2.)
    for b in boxplot_dict['caps']:
        b.set_linewidth(2.)
    ax.set_xlim([0.75, 3.25])
    ax.set_ylabel('CO [ppmv]')
    del data
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'boxplot_co_byenv_%d-%d_%s.eps' %(years[0], years[-1], 
                                                  region), dpi = 300)
    # dictionary of data for NO2 boxplots
    data = {}
    data['a'] = np.hstack(no2_all)
    data['b'] = np.hstack(no2_u)
    data['c'] = np.hstack(no2_su)
    data['d'] = np.hstack(no2_r)
    data['e'] = np.hstack(mr2_no2)
    ax = plt.subplot2grid((1, 1), (0, 0))
    boxplot_dict = ax.boxplot(
            [data[x] for x in ['a', 'b', 'c', 'd', 'e']],
            positions = [1, 1.5, 2, 2.5, 3],
            labels = controls, 
            patch_artist = True,
            widths = 0.25, showfliers = True)
    for b in boxplot_dict['medians']:
        b.set_color('k')
        b.set_linewidth(2.)
    i = 0
    for b in boxplot_dict['boxes']:
        lab = ax.get_xticklabels()[i].get_text()
        b.set_facecolor(color_dict[lab]) 
        b.set_linewidth(2.)
        i += 1
    for b in boxplot_dict['fliers']:
        b.set_marker('.')
        b.set_markersize(1)
        b.set_markerfacecolor('k')
    for b in boxplot_dict['whiskers']:
        b.set_linewidth(2.)
    for b in boxplot_dict['caps']:
        b.set_linewidth(2.)
    ax.set_xlim([0.75, 3.25])
    ax.set_ylabel('NO$_{2}$ [ppbv]')
    del data
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'boxplot_no2_byenv_%d-%d_%s.eps' %(years[0], years[-1], 
                                                   region), dpi = 300)
    # dictionary of data for O3 boxplots
    data = {}
    data['aa'] = np.hstack(o3_all_castnet)
    data['a'] = np.hstack(o3_all) * 1000.
    data['b'] = np.hstack(o3_u) * 1000.
    data['c'] = np.hstack(o3_su) * 1000.
    data['d'] = np.hstack(o3_r) * 1000.
    data['e'] = np.hstack(mr2_o3)
    color_dict = {'CASTNet' : '#a6d854',
                  'AQS$_{\mathregular{ all}}$' : '#66c2a5',
                  'AQS$_{\mathregular{ urban}}$' : '#fc8d62',
                  'AQS$_{\mathregular{ suburban}}$' : '#8da0cb', 
                  'AQS$_{\mathregular{ rural}}$' : '#e78ac3',
                  'MR2' : 'w'}    
    controls = ['CASTNet', 'AQS$_{\mathregular{ all}}$', 
                'AQS$_{\mathregular{ urban}}$', 'AQS$_{\mathregular{ suburban}}$', 
                'AQS$_{\mathregular{ rural}}$', 'MR2']
    ax = plt.subplot2grid((1, 1), (0, 0))
    boxplot_dict = ax.boxplot(
            [data[x] for x in ['aa', 'a', 'b', 'c', 'd', 'e']],
            positions = [1, 1.5, 2, 2.5, 3, 3.5],
            labels = controls, 
            patch_artist = True,
            widths = 0.25, showfliers = True)
    for b in boxplot_dict['medians']:
        b.set_color('k')
        b.set_linewidth(2.)
    i = 0
    for b in boxplot_dict['boxes']:
        lab = ax.get_xticklabels()[i].get_text()
        b.set_facecolor(color_dict[lab]) 
        b.set_linewidth(2.)
        i += 1
    for b in boxplot_dict['fliers']:
        b.set_marker('.')
        b.set_markersize(1)
        b.set_markerfacecolor('k')
    for b in boxplot_dict['whiskers']:
        b.set_linewidth(2.)
    for b in boxplot_dict['caps']:
        b.set_linewidth(2.)
    ax.set_xlim([0.75, 3.75])
    ax.set_ylabel('O$_{3}$ [ppbv]')
    del data
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'boxplot_o3_byenv_%d-%d_%s.eps' %(years[0], years[-1], 
                                                  region), dpi = 300)
    return 
# # # # # # # # # # # # # 
def map_gmiaqsbias(obs, gmi, lat, lon, species, filename):
    """function replaces CTM trace gas concentrations with NaNs if no 
    observation for a particular day/grid cell exist and then finds the time-
    averaged summer (JJA) bias between the CTM and observations. This bias 
    is plotted as a map over the focus region. Values greater than (less than)
    8 (-8) are set to a constant value. 
    
    Parameters
    ----------        
    obs : numpy.ndarray
        Daily 12Z AQS trace gas measurements for trace gas of interest (i.e., 
        NO2, CO, or O3) within CTM grid cells and environment of interest. If 
        no AQS stations are located within the bounds of a CTM grid cell, a nan 
        value is returned. If > 1 AQS station exists in a grid cell, the daily 
        values at these stations are averaged, [time, lat, lon]    
    gmi : numpy.ndarray
        Daily 12Z gridded CTM output for the species of interest at the lowest 
        model pressure level, same units as variable 'obs', [time, lat, lon]     
    lat : numpy.ndarray
        Latitude coordinates of focus region, units of degrees north, [lat,]
    lon : numpy.ndarray
        Longitude coordinates of focus region, units of degrees east, [lon,]
    species : str
        Used for labeling the colorbar of the map
    filename : str
        Used for the filename extension and plot's title, this variable should
        define which environment (i.e., rural, urban, sururban, all) 
        observations are derived from

    Returns
    ----------
    None      
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from matplotlib.colors import rgb2hex, Normalize
    from matplotlib.cm import ScalarMappable
    from matplotlib.patches import Polygon
    from matplotlib.colorbar import ColorbarBase
    # for all missing values in gridded observation (NaNs), replace 
    # correpsonding value in CTM grid with NaN
    nanwhere = np.argwhere(np.isnan(obs))
    gmi[nanwhere[:, 0], nanwhere[:, 1], nanwhere[:, 2]] == np.nan
    # find average difference between observations and model 
    diff = obs - gmi
    diff = np.nanmean(diff, axis = 0)
    # CTM resolution
    rlon = np.diff(lon).mean()
    rlat = np.diff(lat).mean()
    # initialize figure, axis
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    m = Basemap(projection = 'mill', llcrnrlon = -87., llcrnrlat = 33, 
                urcrnrlon = -66., urcrnrlat = 49., resolution = 'i', 
                area_thresh = 10000)   
    # from https://stackoverflow.com/questions/43646097/
    # add-colorbar-for-polygones-in-basemap-plot-of-us
    cmap = plt.cm.RdBu_r
    cmap.set_under((0.019607843137254902, 0.18823529411764706, 0.38039215686274508, 1.0))
    cmap.set_over((0.40392156862745099, 0.0, 0.12156862745098039, 1.0))
    vmin = -8.; vmax = 8.
    # alternatively, the following line could be used if hard coding is to be
    # avoided
    #vmin = np.nanmin(diff); vmax = np.nanmax(diff)
    norm = Normalize(vmin = vmin, vmax = vmax)
    # color mapper to covert values to colors
    mapper = ScalarMappable(norm = norm, cmap = cmap)
    # loop through GMI latitude and longitude coordinatesion
    for i, lat_gb in enumerate(lat):
        for j, lon_gb in enumerate(lon):
            # convert longitude from (0-360) to (-180 to 180)
            lon_gb = np.mod(lon_gb - 180.0, 360.0) - 180.0
            # determine bounding box centered on CTM lat/lon
            lon_left = lon_gb - (rlon/2.)
            lon_right = lon_gb + (rlon/2.)
            lat_top = lat_gb + (rlat/2.)
            lat_bottom = lat_gb - (rlat/2.) 
            # convert to map coordinates
            x1, y1 = m(lon_right, lat_bottom)
            x2, y2 = m(lon_right, lat_top)
            x3, y3 = m(lon_left, lat_top)
            x4, y4 = m(lon_left, lat_bottom)
            # add patch to map
            if np.isnan(diff[i, j]) != True:
                color = rgb2hex(mapper.to_rgba(diff[i, j]))
                p = Polygon([(x1, y1), (x2, y2), (x3, y3), (x4, y4)], 
                             facecolor = color,
                             edgecolor = 'none', linewidth = 0.0, zorder = 5, 
                             alpha = 1.0) 
                ax.add_patch(p) 
    m.drawmapboundary(fill_color='#EBF4FA')
    m.fillcontinents(color = '#CCCCCC', lake_color = '#EBF4FA')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    cax = fig.add_axes([0.25, 0.12, 0.53, 0.05])
    ColorbarBase(cax, cmap = cmap, norm = norm, orientation = 'horizontal', 
                 extend = 'both', label = '$\Delta$ %s' %species)
    plt.subplots_adjust(bottom = 0.2)
    ax.set_title('%s' %filename)
    # to check to ensure code is working:
    #x, y = np.meshgrid(lon, lat)
    #x, y = m(x, y)
    #cb = m.scatter(x, y, c = diff, s = 70, vmax = 8., vmin = -8., cmap = cmap)
    #plt.colorbar(cb)    
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 'map_gmiobsbias_%s.eps' 
                %filename, dpi = 300)
    return
# # # # # # # # # # # # #    
def map_aqsstations(rural_coords, urban_coords, suburban_coords, filename):
    """function plots the locations of urban, rural, and suburban AQS 
    stations measuring the species of interest. 
    
    Parameters
    ----------        
    rural_coords : list
        Latitude and longitude coordinates of AQS stations in rural 
        environments
    urban_coords : list
        Latitude and longitude coordinates of AQS stations in urban and center
        city environments
    suburban_coords : list
        Latitude and longitude coordinates of AQS stations in suburban 
        environments
    filename : str
        Used for the filename extension and plot's title, this variable should
        define which species' AQS stations are plotted

    Returns
    ----------
    None          
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    # initialize axis
    ax = plt.subplot2grid((1, 1), (0, 0))
    m = Basemap(projection = 'mill', llcrnrlon = -87., llcrnrlat = 33, 
                urcrnrlon = -66., urcrnrlat = 49., resolution = 'i', 
                area_thresh = 10000)   
    m.drawmapboundary(fill_color='#EBF4FA')
    m.fillcontinents(color = '#CCCCCC', lake_color = '#EBF4FA')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    # for rural stations
    sc_rural = np.vstack(np.array(rural_coords))
    xsc_rural, ysc_rural = m(sc_rural[:, 1], sc_rural[:, 0])
    m.scatter(xsc_rural, ysc_rural, marker = '.', color = '#4daf4a', 
              zorder = 10, label = 'Rural')
    # for suburban stations
    sc_suburban = np.vstack(np.array(suburban_coords))
    xsc_suburban, ysc_suburban = m(sc_suburban[:, 1], sc_suburban[:, 0])
    m.scatter(xsc_suburban, ysc_suburban, marker = '.', color = '#377eb8', 
              zorder = 10, label = 'Suburban')
    # for urban stations
    sc_urban = np.vstack(np.array(urban_coords))
    xsc_urban, ysc_urban = m(sc_urban[:, 1], sc_urban[:, 0])
    m.scatter(xsc_urban, ysc_urban, marker = '.', color = '#e41a1c', 
              zorder = 10, label = 'Urban')
    plt.legend(frameon = False)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 'map_aqsstations_%s.eps' 
                %filename, dpi = 300)
    return 
# # # # # # # # # # # # #    
import numpy as np
import sys
sys.path.append('/Users/ghkerr/phd/')
import pollutants_constants
sys.path.append('/Users/ghkerr/phd/GMI')
import commensurability
# values for Northeast U.S.
region = 'northeast'
castnet_sites_fr = ['ASH', 'WFM', 'WST', 'APT', 'SAR', 'CTH', 'WSP',
                    'ARE', 'BEL', 'PSU', 'LRL', 'PAR', 'PED', 'SHN', 
                    'CDR', 'VPI', 'MKG', 'KEF']
years = [2008, 2009, 2010]
sampling_months = [6, 7, 8]
sampling_hours = [15, 16, 17, 18, 19, 20]
# load data
# CTM/CASTNet
comm_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2', 
                                          years, sampling_months, 
                                          sampling_hours)  
# AQS at all sites
aqs_co, aqs_no2, aqs_o3, aqs_co_coords, aqs_no2_coords, aqs_o3_coords = \
commensurability.commensurate_aqstracegas(castnet_sites_fr, years, sampling_months, 
                                       sampling_hours)
# AQS at specific environments
(comm_co_su, comm_co_r, comm_co_u, coords_co_su, coords_co_r, coords_co_u, 
 comm_no2_su, comm_no2_r, comm_no2_u, coords_no2_su, coords_no2_r, 
 coords_no2_u, comm_o3_su, comm_o3_r, comm_o3_u, coords_o3_su, coords_o3_r, 
 coords_o3_u) = commensurability.commensurate_aqstracegas_siting(
 castnet_sites_fr, years, sampling_months, sampling_hours)
# average over sites 
# for CO
co_all = np.nanmean(aqs_co, axis = 1)
co_u = np.nanmean(comm_co_u, axis = 1)
co_su = np.nanmean(comm_co_su, axis = 1)
co_r = np.nanmean(comm_co_r, axis = 1)
mr2_co = np.nanmean(mr2_co, axis = 1)
# for NO2
no2_all = np.nanmean(aqs_no2, axis = 1)
no2_u = np.nanmean(comm_no2_u, axis = 1)
no2_su = np.nanmean(comm_no2_su, axis = 1)
no2_r = np.nanmean(comm_no2_r, axis = 1)
mr2_no2 = np.nanmean(mr2_no2, axis = 1)
# for O3
o3_all_castnet = np.nanmean(comm_castnet, axis = 1)
o3_all = np.nanmean(aqs_o3, axis = 1)
o3_u = np.nanmean(comm_o3_u, axis = 1)
o3_su = np.nanmean(comm_o3_su, axis = 1)
o3_r = np.nanmean(comm_o3_r, axis = 1)
mr2_o3 = np.nanmean(mr2_o3, axis = 1)
 plot boxplots of AQS O3, CO, and NO2 by environments
boxplot_tracegas_byenv(co_all, co_u, co_su, co_r, mr2_co, no2_all, 
                       no2_u, no2_su, no2_r, mr2_no2, o3_all_castnet, 
                       o3_all, o3_u, o3_su, o3_r, mr2_o3, years, region)
# open 12Z CTM output for focus region
(lat, lon, pressure, times, co, no, no2, o3) = \
commensurability.open_gridded_idailyCTM(years)
# load AQS observations
COdf, NO2df, O3df = commensurability.load_aqshourly(years)
# retrieve ALL AQS observations in CTM grid cells for NO2
no2_station_coordinates, colocated_no2obs = \
commensurability.commensurate_aqstracegas_gridded(NO2df, no2, times, lat, lon, 
                                                  [])
map_gmiaqsbias(colocated_no2obs, no2[:, 0] * 1e9, lat, lon, 'NO$_{2}$', 
               'All')
# retrieve RURAL AQS observations in CTM grid cells
no2_station_coordinates_rural, colocated_no2obs_rural = \
commensurability.commensurate_aqstracegas_gridded(NO2df, no2, times, lat, lon, 
                                                  ['RURAL'])
map_gmiaqsbias(colocated_no2obs_rural, no2[:, 0] * 1e9, lat, lon, 'NO$_{2}$', 
               'Rural')
# retrieve URBAN AQS observations in CTM grid cells
no2_station_coordinates_urban, colocated_no2obs_urban = \
commensurability.commensurate_aqstracegas_gridded(NO2df, no2, times, lat, lon, 
                                                  ['URBAN AND CENTER CITY'])
map_gmiaqsbias(colocated_no2obs_urban, no2[:, 0] * 1e9, lat, lon, 'NO$_{2}$', 
               'Urban')
# retrieve SUBURAN AQS observations in CTM grid cells
no2_station_coordinates_suburban, colocated_no2obs_suburban = \
commensurability.commensurate_aqstracegas_gridded(NO2df, no2, times, lat, lon, 
                                                  ['SUBURBAN'])
map_gmiaqsbias(colocated_no2obs_suburban, no2[:, 0] * 1e9, lat, lon, 
               'NO$_{2}$', 'Suburban')
# plot locations of rural, urban, and suburban AQS stations
map_aqsstations(no2_station_coordinates_rural, no2_station_coordinates_urban, 
                no2_station_coordinates_suburban, 'no2')
# retrieve ALL AQS observations in CTM grid cells for CO
co_station_coordinates, colocated_coobs = \
commensurability.commensurate_aqstracegas_gridded(COdf, co, times, lat, lon, 
                                                  [])
map_gmiaqsbias(colocated_coobs, co[:, 0] * 1e6, lat, lon, 'CO', 
               'All')
# retrieve RURAL AQS observations in CTM grid cells
co_station_coordinates_rural, colocated_coobs_rural = \
commensurability.commensurate_aqstracegas_gridded(COdf, co, times, lat, lon, 
                                                  ['RURAL'])
map_gmiaqsbias(colocated_coobs_rural, co[:, 0] * 1e6, lat, lon, 'CO', 
               'Rural')
# retrieve URBAN AQS observations in CTM grid cells
co_station_coordinates_urban, colocated_coobs_urban = \
commensurability.commensurate_aqstracegas_gridded(COdf, co, times, lat, lon, 
                                                  ['URBAN AND CENTER CITY'])
map_gmiaqsbias(colocated_coobs_urban, co[:, 0] * 1e6, lat, lon, 'CO', 
               'Urban')
# retrieve SUBURAN AQS observations in CTM grid cells
co_station_coordinates_suburban, colocated_coobs_suburban = \
commensurability.commensurate_aqstracegas_gridded(COdf, co, times, lat, lon, 
                                                  ['SUBURBAN'])
map_gmiaqsbias(colocated_coobs_suburban, co[:, 0] * 1e6, lat, lon, 
               'CO', 'Suburban')
# plot locations of rural, urban, and suburban AQS stations
map_aqsstations(co_station_coordinates_rural, co_station_coordinates_urban, 
                co_station_coordinates_suburban, 'co')