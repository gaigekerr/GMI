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
import matplotlib as mpl
import matplotlib.font_manager as font_manager    
import sys
sys.path.append('/Users/ghkerr/phd/')
import pollutants_constants
# set custom font
path = pollutants_constants.FONTPATH_LIGHT
prop = font_manager.FontProperties(fname=path)
mpl.rcParams['font.family'] = prop.get_name()
path = pollutants_constants.FONTPATH_BOLD
prop = font_manager.FontProperties(fname = path)
mpl.rcParams['mathtext.bf'] = prop.get_name()
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







































#
#plt.plot(o3_u, '-k'); plt.plot(o3_su, '-r'); plt.plot(o3_r, '-b')
#plt.show()
#
#
#
#plt.plot(co_u, '-k'); plt.plot(co_su, '-r'); plt.plot(co_r, '-b'); plt.plot(mr2_co, '-g')
#plt.show()
#
#
#
#
#plt.plot(no2_u, '-k'); plt.plot(no2_su, '-r'); plt.plot(no2_r, '-b'); 
#plt.twinx()
#plt.plot(mr2_no2, '-g')
#plt.show()