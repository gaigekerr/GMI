#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAME
    geos_chem_compare.py
PURPOSE
    Functions compare GOES-Chem output (from Katie Travis) with GMI CTM output
    for 1 August - 30 September 2013 period
PROGRAMMER
    Gaige Hunter Kerr
REVISION HISTORY
    10052018 -- initial version created
    16052018 -- functions 'diurnal_tracegases' and 'map_geosgmiaqscastnet' 
                added
"""
# global constants 
COLOR_GMI = '#66c2a5'
COLOR_GEOS = '#fc8d62'  
# # # # # # # # # # # # #
def timeseries_castneto3gmio3geoso3(castnet, mr2_o3, geos_o3):
    """plot timeseries of CASTNet, GMI CTM, and GEOS-Chem O3 for August - 
    September 2013. Show locations of weekends to see how GEOS-Chem handles 
    weekend effect differently than GMI CTM. 

    Parameters
    ----------    
    castnet : numpy.ndarray
        Regionally-averaged O3 observations from CASTNET, units of ppbv, 
        (61,)
    mr2_o3 : numpy.ndarray
        Regionally-averaged modeled O3 concentrations from GMI CTM HindcastMR2
        simulation, units of ppbv, (61,)
    geos_o3 : numpy.ndarray
        Regionally-averaged modeled O3 concentrations from GEOS-Chem, units of 
        ppbv, (61,)

    Returns
    ----------      
    None
    """    
    import numpy as np
    import datetime
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import matplotlib.font_manager as font_manager    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # set custom font
    path = pollutants_constants.FONTPATH_LIGHT
    prop = font_manager.FontProperties(fname=path)
    mpl.rcParams['font.family'] = prop.get_name()
    mpl.rcParams['axes.unicode_minus'] = False
    path = pollutants_constants.FONTPATH_BOLD
    prop = font_manager.FontProperties(fname = path)
    mpl.rcParams['mathtext.bf'] = prop.get_name()
    # find days in measuring period August - September 2013
    start = datetime.datetime.strptime("01-08-2013", "%d-%m-%Y")
    end = datetime.datetime.strptime("30-09-2013", "%d-%m-%Y")
    end = end + datetime.timedelta(days = 1)
    dates = [start + datetime.timedelta(days = x) for x in 
             range(0, (end-start).days)]
    # find which days are weekdays, datetime.weekday() returns the day of the 
    # week as an integer, where Monday is 0 and Sunday is 6
    weekdays = [x.weekday() for x in dates]
    where_weekend = np.where((np.array(weekdays) == 5) | 
                             (np.array(weekdays) == 6))[0]
    # initialize figure, axis 
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    # plot observed and modeled (GMI CTM, GEOS-Chem) O3
    ax.plot(castnet, '-k', lw = 2., label = 'CASTNet')
    ax.plot(mr2_o3, linestyle = '-', color = COLOR_GMI, lw = 1.5, 
            zorder = 10, label = 'GMI CTM')
    ax.plot(geos_o3, linestyle = '-', color = COLOR_GEOS, lw = 1.5, 
            zorder = 10, label = 'GEOS-Chem')
    # add legend 
    leg = ax.legend(bbox_to_anchor = (1.02, 0.65))
    leg.get_frame().set_linewidth(0.0)
    # shade locations of weekends; n.b. since shading 'axvline' starts at 
    # weekend index, add 0.5 so it centers shaded region in middle of the 
    # index
    for weekend_day in where_weekend:
        ax.axvspan(weekend_day - 0.5, weekend_day + 0.5, edgecolor = 'gainsboro', 
                   facecolor = 'gainsboro', alpha = 1., zorder = 1)
    ax.set_xlim([0, len(castnet) - 1])
    ax.set_xticks([0, 14, 31, 45])
    ax.set_xticklabels(['1 Aug 2013', '15 Aug 2013', '1 Sep 2013', '15 Sep 2013'])
    ax.set_ylabel('O$_{3}$ [ppbv]')
    plt.subplots_adjust(right = 0.75)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/GEOS-Chem/' + 
                'timeseries_castneto3gmio3geoso3.eps', dpi = 300)
    return 
# # # # # # # # # # # # #
def timeseries_tracegases(castnet, mr2_o3, geos_o3, mr2_no, geos_no, no2, 
                          mr2_no2, geos_no2, co, mr2_co, geos_co):
    """plot timeseries of regionally-averaged observed and modeled (GMI CTM and
    GEOS-Chem) O3, NO, NO2, and CO. Observations for O3 (NO2/CO) are from 
    CASTNet (AQS). 
    
    Parameters
    ----------    
    castnet : numpy.ndarray
        Regionally-averaged O3 observations from CASTNET, units of ppbv, 
        (61,)
    mr2_o3 : numpy.ndarray
        Regionally-averaged modeled O3 concentrations from GMI CTM HindcastMR2
        simulation, units of ppbv, (61,)
    geos_o3 : numpy.ndarray
        Regionally-averaged modeled O3 concentrations from GEOS-Chem, units of 
        ppbv, (61,)
    mr2_no : numpy.ndarray
        Regionally-averaged modeled NO concentrations from GMI CTM HindcastMR2
        simulation, units of ppbv, (61,)    
    geos_no : numpy.ndarray 
        Regionally-averaged modeled NO concentrations from GEOS-Chem, units of 
        ppbv, (61,)    
    no2 : numpy.ndarray
        Regionally-averaged NO2 observations from AQS, units of ppbv, (61,)    
    mr2_no2 : numpy.ndarray
    geos_no2 : numpy.ndarray 
    co : numpy.ndarray
        Regionally-averaged CO observations from AQS, units of ppmv, (61,)    
    mr2_co : numpy.ndarray
        Regionally-averaged modeled CO concentrations from GMI CTM HindcastMR2
        simulation, units of ppmv, (61,)        
    geos_co : numpy.ndarray 
        Regionally-averaged modeled CO concentrations from GEOS-Chem, units of 
        ppbv, (61,)    

    Returns
    ----------      
    None    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import matplotlib.font_manager as font_manager    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # set custom font
    path = pollutants_constants.FONTPATH_LIGHT
    prop = font_manager.FontProperties(fname=path)
    mpl.rcParams['font.family'] = prop.get_name()
    mpl.rcParams['axes.unicode_minus'] = False
    path = pollutants_constants.FONTPATH_BOLD
    prop = font_manager.FontProperties(fname = path)
    mpl.rcParams['mathtext.bf'] = prop.get_name()
    # initialize figure, axes
    fig = plt.figure()
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax3 = plt.subplot2grid((2, 2), (1, 0))
    ax4 = plt.subplot2grid((2, 2), (1, 1))
    # plot O3
    ax1.plot(castnet, '-k', lw = 2., label = 'CASTNet/AQS')
    ax1.plot(mr2_o3, linestyle = '-', color = COLOR_GMI, lw = 1.5, 
            zorder = 10, label = 'GMI CTM')
    ax1.plot(geos_o3, linestyle = '-', color = COLOR_GEOS, lw = 1.5, 
            zorder = 10, label = 'GEOS-Chem')
    ax1.set_xlim([0, len(castnet) - 1])
    ax1.set_xticks([0, 14, 31, 45])
    ax1.set_xticklabels([''])
    ax1.set_ylabel('O$_{3}$ [ppbv]')
    # plot NO
    ax2.plot(mr2_no, linestyle = '-', color = COLOR_GMI, lw = 1.5, 
            zorder = 10)
    ax2.plot(geos_no, linestyle = '-', color = COLOR_GEOS, lw = 1.5, 
             zorder = 10)
    ax2.set_xlim([0, len(castnet) - 1])
    ax2.set_xticks([0, 14, 31, 45])
    ax2.set_xticklabels([''])
    ax2.set_ylabel('NO [ppbv]')
    # plot NO2
    ax3.plot(no2, '-k', lw = 2.)
    ax3.plot(mr2_no2, linestyle = '-', color = COLOR_GMI, lw = 1.5, 
            zorder = 10)
    ax3.plot(geos_no2, linestyle = '-', color = COLOR_GEOS, lw = 1.5, 
             zorder = 10)
    ax3.set_xlim([0, len(castnet) - 1])
    ax3.set_xticks([0, 14, 31, 45])
    ax3.set_xticklabels(['1 Aug', '15 Aug', '1 Sep', '15 Sep'])
    ax3.set_ylabel('NO$_{2}$ [ppbv]')
    # plot CO
    ax4.plot(co, '-k', lw = 2.)
    ax4.plot(mr2_co, linestyle = '-', color = COLOR_GMI, lw = 1.5, 
            zorder = 10)
    ax4.plot(geos_co/1000., linestyle = '-', color = COLOR_GEOS, lw = 1.5, 
             zorder = 10)
    ax4.set_xticks([0, 14, 31, 45])
    ax4.set_xticklabels(['1 Aug', '15 Aug', '1 Sep', '15 Sep'])
    ax4.set_ylabel('CO [ppmv]')
    plt.tight_layout()
    plt.subplots_adjust(bottom = 0.2)
    # add legend
    leg = ax1.legend(bbox_to_anchor = (2.1, -1.35), ncol = 3)
    leg.get_frame().set_linewidth(0.0)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/GEOS-Chem/' + 
                'timeseries_tracegases.eps', dpi = 300)    
    return
# # # # # # # # # # # # #

def diurnal_tracegases(castnet_d, geos_o3_d, gmi_o3_d, geos_no_d, gmi_no_d, 
                       no2_d, geos_no2_d, gmi_no2_d, co_d, geos_co_d, gmi_co_d):
    """function plots diurnal curves of regionally-averaged observed and 
    modeled (GMI CTM and GEOS-Chem) O3, NO, NO2, and CO. Observations for O3 
    (NO2/CO) are from CASTNet (AQS). 
        
    Parameters
    ----------    
    castnet_d : numpy.ndarray
        Diurnal curve of regionally-averaged O3 observations from CASTNET, 
        units of ppbv, (24,)
    geos_o3_d : numpy.ndarray
        Diurnal curve of regionally-averaged O3 observations from GEOS-Chem, 
        units of ppbv, (12,)    
    gmi_o3_d : numpy.ndarray
        Diurnal curve of regionally-averaged O3 observations from GMI CTM,
        units of ppbv, (24,)       
    geos_no_d : numpy.ndarray
        Diurnal curve of regionally-averaged NO observations from GEOS-Chem, 
        units of ppbv, (12,)      
    gmi_no_d : numpy.ndarray
        Diurnal curve of regionally-averaged NO observations from GMI CTM,
        units of ppbv, (24,)           
    no2_d : numpy.ndarray
        Diurnal curve of regionally-averaged NO observations from AQS, units 
        of ppbv, (24,)        
    geos_no2_d : numpy.ndarray
        Diurnal curve of regionally-averaged NO2 observations from GEOS-Chem, 
        units of ppbv, (12,)          
    gmi_no2_d : numpy.ndarray
        Diurnal curve of regionally-averaged NO2 observations from GMI CTM,
        units of ppmv, (24,)      
    co_d : numpy.ndarray
        Diurnal curve of regionally-averaged CO observations from AQS, units 
        of ppmv, (24,)        
    geos_co_d : numpy.ndarray
        Diurnal curve of regionally-averaged CO observations from GEOS-Chem, 
        units of ppbv, (12,)         
    gmi_co_d : numpy.ndarray  
        Diurnal curve of regionally-averaged CO observations from GMI CTM,
        units of ppbv, (24,)      
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import matplotlib.font_manager as font_manager    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # set custom font
    path = pollutants_constants.FONTPATH_LIGHT
    prop = font_manager.FontProperties(fname=path)
    mpl.rcParams['font.family'] = prop.get_name()
    mpl.rcParams['axes.unicode_minus'] = False
    path = pollutants_constants.FONTPATH_BOLD
    prop = font_manager.FontProperties(fname = path)
    mpl.rcParams['mathtext.bf'] = prop.get_name()
    # shift CASTNet O3 observations ahead 4 hours so that they're in the same
    # time as CTMs (i.e. UTC)
    castnet_d = np.roll(castnet_d, 4)
    # initialize figure, axes
    fig = plt.figure()
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax3 = plt.subplot2grid((2, 2), (1, 0))
    ax4 = plt.subplot2grid((2, 2), (1, 1))
    # plot O3
    ax1.plot(np.arange(0, 23, 2), geos_o3_d, lw = 1.5, color = COLOR_GEOS, 
             label = 'GEOS-Chem') # since GEOS-Chem output is 2-hrly, only 
    # plot every two hours so it is same length as CASTNet
    # and GMI
    ax1.plot(gmi_o3_d, lw = 1.5, color = COLOR_GMI, label = 'GMI CTM')
    ax1.plot(castnet_d, lw = 2.0, color = 'k', label = 'CASTNet/AQS')
    ax1.set_xlim([0, 23])
    ax1.set_xticks(np.arange(0, 24, 6))
    ax1.set_xticklabels([''])
    ax1.set_ylabel('O$_{3}$ [ppbv]')
    # plot NO
    ax2.plot(np.arange(0, 23, 2), geos_no_d, lw = 1.5, color = COLOR_GEOS)
    ax2.plot(gmi_no_d, lw = 1.5, color = COLOR_GMI)
    ax2.set_xlim([0, 23])
    ax2.set_xticks(np.arange(0, 24, 6))
    ax2.set_xticklabels([''])
    ax2.set_ylabel('NO [ppbv]')
    # plot NO2
    ax3.plot(np.arange(0, 23, 2), geos_no2_d, lw = 1.5, color = COLOR_GEOS)
    ax3.plot(gmi_no2_d, lw = 1.5, color = COLOR_GMI)
    ax3.plot(no2_d, lw = 2.0, color = 'k')
    ax3.set_xlim([0, 23])
    ax3.set_xticks(np.arange(0, 24, 6))
    ax3.set_xticklabels(np.arange(0, 24, 6))
    ax3.set_xlabel('Hour [GMT]')
    ax3.set_ylabel('NO$_{2}$ [ppbv]')
    # plot CO
    ax4.plot(np.arange(0, 23, 2), geos_co_d/1000., lw = 1.5, color = COLOR_GEOS)
    ax4.plot(gmi_co_d, lw = 1.5, color = COLOR_GMI)
    ax4.plot(co_d, lw = 2.0, color = 'k')
    ax4.set_xlim([0, 23])
    ax4.set_xticks(np.arange(0, 24, 6))
    ax4.set_xticklabels(np.arange(0, 24, 6))
    ax4.set_xlabel('Hour [GMT]')
    ax4.set_ylabel('CO [ppmv]')
    plt.tight_layout()
    plt.subplots_adjust(bottom = 0.2)
    # add legend
    leg = ax1.legend(bbox_to_anchor = (2.1, -1.50), ncol = 3)
    leg.get_frame().set_linewidth(0.0)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/GEOS-Chem/' + 
                'diurnal_tracegases.eps', dpi = 300)  
    return 
# # # # # # # # # # # # #    
def map_geosgmiaqscastnet(castnet_sites_fr, sampling_months, sampling_hours, 
                          geos_lons, geos_lats, aqs_no2_lons, aqs_no2_lats, 
                          aqs_co_lons, aqs_co_lats):
    """for the CASTNet sites specified in variable 'castnet_sites_fr,' function
    opens the locations (latitude/longitude coordinates) of those CASTNet 
    sites and their co-located (or nearly co-located) MERRA grid cells, GMI 
    sites, and AQS CO and NO2 sites and plots their locations. 
    
    Parameters
    ----------    
    castnet_sites_fr : list
        CASTNET site names in focus region
    sampling_months : list 
        Months of interest
    sampling_hours : list 
        Hours during which trace gas concentrations are returned
    geos_lons : list
        Longitude of GEOS-Chem grid cells co-located with CASTNet sites       
    geos_lats : list
        Latitude of GEOS-Chem grid cells co-located with CASTNet sites    
    aqs_no2_lons : list
        Longitudes of AQS sites within the bounding box defined by the CASTNet 
        station's longitude and latitude +/- variable 'searchrad' which 
        measure NO2 during the years and months of interest
    aqs_no2_lats : list 
        Latitudes of AQS sites within the bounding box defined by the CASTNet 
        station's longitude and latitude +/- variable 'searchrad' which 
        measure NO2 during the years and months of interest    
    aqs_co_lons : list
        Longitudes of AQS sites within the bounding box defined by the CASTNet 
        station's longitude and latitude +/- variable 'searchrad' which 
        measure CO during the years and months of interest          
    aqs_co_lats : list
        Latitudes of AQS sites within the bounding box defined by the CASTNet 
        station's longitude and latitude +/- variable 'searchrad' which 
        measure CO during the years and months of interest    
    
    Returns
    ----------
    None    
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import matplotlib.font_manager as font_manager 
    import matplotlib as mpl    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/stagnation/')
    from plot_focus_region import plot_focus_regionOLD
    sys.path.append('/Users/ghkerr/phd/GMI/')
    import commensurability
    # set custom font
    path = pollutants_constants.FONTPATH_LIGHT
    prop = font_manager.FontProperties(fname=path)
    mpl.rcParams['font.family'] = prop.get_name()
    path = pollutants_constants.FONTPATH_BOLD
    prop = font_manager.FontProperties(fname = path)
    mpl.rcParams['mathtext.bf'] = prop.get_name()
    #
    # # # # open CASTNet and GMI sites
    castnet_lats, castnet_lons, gmi_lats, gmi_lons = \
    commensurability.commensurate_castnet_gmi_locations(castnet_sites_fr, 
        sampling_months, sampling_hours)
    # initialize figure, axis, and a new instance of basemap
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    m = Basemap(projection = 'merc', llcrnrlon = -84., llcrnrlat = 36, 
                urcrnrlon = -65., urcrnrlat = 49., resolution = 'c', 
                area_thresh = 10000)
    # plot GEOS-Chem grid cells 
    x_geos, y_geos = m(geos_lons, geos_lats)
    geos_loc = m.scatter(x_geos, y_geos, 25, color = '#fb9a99', marker = 's', 
                        edgecolor = 'none', linewidth = 0.85, zorder = 20, 
                        label = 'GEOS-Chem')
    # plot GMI sites 
    x_gmi, y_gmi = m(gmi_lons, gmi_lats)
    gmi_loc = m.scatter(x_gmi, y_gmi, 25, color = '#1f78b4', marker = 'x', 
                        edgecolor = 'none', linewidth = 0.85, zorder = 24, 
                        label = 'GMI CTM')
    # plot CASTNet stations
    x_castnet, y_castnet = m(castnet_lons, castnet_lats)
    castnet_loc = m.scatter(x_castnet, y_castnet, 20, color = 'k', marker = '^', 
                           edgecolor = 'k', linewidth = 0.5, zorder = 23, 
                           label = 'CASTNet')
    # plot AQS stations measuring NO2
    aqs_no2_coords = np.array([np.hstack(aqs_no2_lons), 
                               np.hstack(aqs_no2_lats)], dtype = float).T
    aqs_no2_coords = pd.DataFrame(aqs_no2_coords)
    aqs_no2_coords = aqs_no2_coords.drop_duplicates()
    x_aqs_no2, y_aqs_no2 = m(aqs_no2_coords[0].values, 
                             aqs_no2_coords[1].values)
    aqs_no2_loc = m.scatter(x_aqs_no2, y_aqs_no2, 8, color = '#33a02c', 
                            marker = 'o', edgecolor = 'none', linewidth = 0.5, 
                            zorder = 21, label = 'AQS NO$_{2}$')
    # plot AQS stations measuring CO
    aqs_co_coords = np.array([np.hstack(aqs_co_lons), 
                              np.hstack(aqs_co_lats)], dtype = float).T
    aqs_co_coords = pd.DataFrame(aqs_co_coords)
    aqs_co_coords = aqs_co_coords.drop_duplicates()
    x_aqs_co, y_aqs_co = m(aqs_co_coords[0].values, aqs_co_coords[1].values)
    aqs_co_loc = m.scatter(x_aqs_co, y_aqs_co, 8, color = '#b2df8a', marker = 'd', 
                           edgecolor = 'none', linewidth = 0.5, zorder = 22, 
                           label = 'AQS CO')
    # overlay shapefile of Northeast
    plot_focus_regionOLD(ax, m)
    plt.tight_layout()
    # add legend
    leg = ax.legend(handles = [gmi_loc, geos_loc, castnet_loc, aqs_no2_loc, 
                    aqs_co_loc], bbox_to_anchor = (0.8, 0.5), 
                    ncol = 1, fontsize = 12, scatterpoints = 1)
    leg.get_frame().set_linewidth(0.0)
    # remove axis frame from plot             
    ax.axis('off')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/GEOS-Chem/' + 
                'map_geosgmiaqscastnet.eps', dpi = 300)
    return     
# # # # # # # # # # # # #    
#import numpy as np
#import sys
#sys.path.append('/Users/ghkerr/phd/GMI/')
#import commensurability
#from calculate_regional_average import calculate_regional_average
#years = [2013]        
#sampling_months = [8, 9]
#sampling_hours = [15, 16, 17, 18, 19, 20]  
#castnet_sites_fr = ['ASH', 'WFM', 'WST', 'APT', 'SAR', 'CTH', 'WSP',
#                    'ARE', 'BEL', 'PSU', 'LRL', 'PAR', 'PED', 'SHN', 
#                    'CDR', 'VPI', 'MKG', 'KEF']
#
## # # # open HindcastMR2
#castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
#commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2', 
#                                          years, sampling_months, 
#                                          sampling_hours)
## # # # open diurnal curves of HindcastMR2 and average over region and 
## measuring period
#castnet_d, gmi_o3_d, gmi_no_d, gmi_no2_d, gmi_co_d = \
#commensurability.commensurate_castnet_gmi_diurnal(castnet_sites_fr,
#                                                  'HindcastMR2', years, 
#                                                  sampling_months)
#castnet_d, gmi_o3_d, gmi_no_d, gmi_no2_d, gmi_co_d = \
#calculate_regional_average(castnet_d[0], gmi_o3_d[0], gmi_no_d[0], 
#                           gmi_no2_d[0], gmi_co_d[0], (0, 1))
## # # # open diurnal curves of GEOS-Chem and average over region and 
## measuring period 
#temp, geos_o3_d, geos_no_d, geos_no2_d, geos_co_d = \
#commensurability.commensurate_geos_gmi_diurnal(castnet_sites_fr)
#temp, geos_o3_d, geos_no_d, geos_no2_d, geos_co_d = \
#calculate_regional_average(temp, geos_o3_d, geos_no_d, geos_no2_d, geos_co_d, (0, 1))
#del temp
## # # # open diurnal curves of AQS NO2 and CO and average over region and
## measuring period
#aqs_co_d, aqs_no2_d, aqs_no2_lons, aqs_no2_lats, aqs_co_lons, aqs_co_lats = \
#commensurability.commensurate_aqsno2co_diurnal(castnet, castnet_sites_fr, 
#                                               years, sampling_months)
#co_d = np.nanmean(aqs_co_d[0], axis = (0, 1))
#no2_d = np.nanmean(aqs_no2_d[0], axis = (0, 1))
## # # # average HindcastMR2 over region; n.b. didn't do this before 
## because function 'commensurate_aqsno2co_diurnal' needs unaveraged arrays
#castnet, mr2_o3, mr2_no, mr2_no2, mr2_co = \
#calculate_regional_average(castnet[0], mr2_o3[0], mr2_no[0], mr2_no2[0], mr2_co[0], 0)
## # # # open GEOS-Chem and average over region 
#temp, geos_o3, geos_no, geos_no2, geos_co, geos_lats, geos_lons = \
#commensurability.commensurate_geos_gmi(castnet_sites_fr, sampling_hours)
#temp, geos_o3, geos_no, geos_no2, geos_co = \
#calculate_regional_average(temp, geos_o3, geos_no, geos_no2, geos_co, 0)
#del temp
## # # # open AQS NO2 and CO and average over region
#aqs_co, aqs_no2, aqs_co_coords, aqs_no2_coords = \
#commensurability.commensurate_aqsno2co(castnet_sites_fr, years, 
#                                       sampling_months, sampling_hours)
#co = np.nanmean(aqs_co[0], axis = 0)
#no2 = np.nanmean(aqs_no2[0], axis = 0)
# # # # visualizations
#timeseries_castneto3gmio3geoso3(castnet, mr2_o3, geos_o3)
#timeseries_tracegases(castnet, mr2_o3, geos_o3, mr2_no, geos_no, no2, 
#                      mr2_no2, geos_no2, co, mr2_co, geos_co)
#diurnal_tracegases(castnet_d, geos_o3_d, gmi_o3_d, geos_no_d, gmi_no_d, 
#                   no2_d, geos_no2_d, gmi_no2_d, co_d, geos_co_d, gmi_co_d)
#map_geosgmiaqscastnet(castnet_sites_fr, sampling_months, sampling_hours, 
#                      geos_lons, geos_lats, aqs_no2_lons, aqs_no2_lats, 
#                      aqs_co_lons, aqs_co_lats)