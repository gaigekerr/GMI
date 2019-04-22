#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAME
    do3dt.py
PURPOSE
    Investigate the ozone-temperature sensitivty (dO3/dT) in observations 
    and GMI CTM across the continental U.S.
PROGRAMMER
    Gaige Hunter Kerr
REVISION HISTORY
    11092018 -- initial version created
    18092018 -- function 'find_conus_castnet' added
    02102018 -- functions 'calculate_gmi_r_do3dt2m' and 
                'map_ro3t2m_do3dt2m_conus_gmi' added
    04102018 -- change colorscheme in function 'map_ro3t2m_do3dt2m_conus_gmi'
    05102018 -- change function 'castnet_r_do3d2t' to return daily values 
                of O3 and T2m at each CASTNet site with data
    09102018 -- function 'fill_oceans' added; function 'map_std' changed to 
                'map_std_90ptile' to include plots showing the 90th percentile 
                as a metric of variability
    11102018 -- funciton 'castnet_r_do3d2t' edited to sample CASTNet stations 
                only at 1300 hours local time (time of overpass); function 
                'map_ro3t2m_do3dt2m_conus' edited to include both O3-T2m 
                sensitivity and correlation coefficient side-by-side
    12102018 -- function 'timeseries_mr2o3dato3t2m_atpoint' added
    14102018 -- function 'map_meanmr2o3meancastneto3_conus' added
    17102018 -- after adjusting CTM grid/domain size, timezonefinder was 
                giving errors so adjusted the delta_degree from 10 to 11 in 
                function 'calculate_gmi_r_do3dt2m'; function 
                'calculate_gmi_r_do3dt2m_regionmean' added
    18102018 -- functions 'calculate_castnet_r_do3dt2m_regionmean' and 
                'scatter_castnett2mo3_gmit2mo3_slopes' added
    19102018 -- functions 'timeseries_castneto3allgmio3' and 
                'scatterhist_castneto3allgmio3' added
    24102018 -- function 'map_allgmio3_percentcontribution' added
    25102018 -- edited 'map_allgmio3_percentcontribution' to calculated 
                Northeast-averaged values for differences on hot/cold days  
    27102018 -- printed regionally-averaged slopes in function 
                'scatterhist_castneto3allgmio3' removed because they were 
                showing the slope fit through regionally-averaged O3 and 
                2-meter temperature whereas we're interested in the slope 
                averaged over the region (not calculated with regionally-
                averaged quantities)
    28102018 -- function 'calculate_statistics_region' added
    29102018 -- function 'map_sensitivityratio_conus' edit to include a map
                of the O3 and 2-meter temperature correlation and the slope 
                from the Transport simulation in addition to the ratio of the 
                slopes from the Transport and + Chemistry simulations
    05112018 -- function 'timeseries_t2m_castneto3_cemsnox' added
    06112018 -- functions 'find_bias_nationwide_region' and 'map_r2o3t2m_conus'
                added
    08112018 -- add total least squares (ODR) to analysis in function 
                scatter_castnett2mo3_gmit2mo3_slopes
    11112018 -- fixed calculation of standard deviation in function 
                'calculate_statistics_region.' Before function had calculated 
                the regional averaged over the wrong dimensions (time dim) so
                the standard deviation wasn't representative of daily averages
                of the region. This section of code was migrated to function
                to 'calculate_gmi_r_do3dt2m_regionmean.'
    12112018 -- function 'timeseries_mr2o3dato3t2m_atpoint' edited to include
                timeseries of O3, T from a grid cell with high O3-T correlation
                and one with low O3-T correlation
    13112018 -- functions 'map_simulationschematic' and 
                'boxplot_cemsnox_castneto3_neus' added
    26112018 -- maximum and minimum biases and which CASTNet they occur at 
                recorded in function 'map_meanmr2o3meancastneto3_conus' 
    28112018 -- function 'scatter_inventorynonox_noxo3' added
    11122018 -- function 'timeseries_t2m_castneto3_cemsnox' edited and scatter-
                plot removed and model timeseries added
    20022019 -- function 'timeseries_map_hourlyvsoverpass_neus' added; 
                function 'timeseries_map_hourlyvsoverpass_neus' edited to 
                output mean O3 and bias plots as a two-panel figure
    21022019 -- function 'scatter_inventorynoo3' added to replace, in effect, 
                'scatter_inventorynonox_noxo3' (simplified version)
    22022019 -- function 'timeseries_transportchemistryo3_atpoints' added
    25022019 -- function 'scatter_inventorynoo3' changed to focus on O3 vs. NOx
    22042019 -- function 'map_allgmio3_do3dt_byprecip' added to start dealing
                with reviewers' comments...
"""
# # # # # # # # # # # # #
# change font
import matplotlib as mpl
import scipy as sp
import matplotlib.pyplot as plt    
prop = mpl.font_manager.FontProperties(fname = 
    '/Users/ghkerr/Library/Fonts/cmunbmr.ttf')
mpl.rcParams['font.family'] = prop.get_name()
prop = mpl.font_manager.FontProperties(fname = 
    '/Users/ghkerr/Library/Fonts/cmunbbx.ttf')
mpl.rcParams['mathtext.bf'] = prop.get_name()
# for unicode minus/negative sign implementation
mpl.rcParams['axes.unicode_minus'] = False
# change width and thickness of ticks/spines
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.width'] = 1.5    
# # # # # # # # # # # # #
class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)
    def __call__(self, value, clip=None):
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return sp.ma.masked_array(sp.interp(value, x, y))
# # # # # # # # # # # # #
def fill_oceans(ax, m):
    """fill oceans with solid color on Basemap; using maskoceans in 
    matplotlib.basemap leads to a white ocean mask with a resolution much lower
    than the coastline. This method is faster and cleaner.
    
    Parameters
    ----------  
    ax : matplotlib.axes._subplots.AxesSubplot
        Axis containing Basemap class instance
    m : mpl_toolkits.basemap.Basemap
        Basemap with specified coordinates and map projection        

    Returns
    ----------     
    None
    """
    import numpy as np
    from matplotlib.patches import Path, PathPatch
    # get the limits of the map:
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    map_edges = np.array([[x0,y0],[x1,y0],[x1,y1],[x0,y1]])
    # get all polygons used to draw the coastlines of the map
    polys = [p.boundary for p in m.landpolygons]
    # combining with map edges
    polys = [map_edges]+polys[:]
    # creating a PathPatch
    codes = [[Path.MOVETO] + [Path.LINETO for p in p[1:]]
            for p in polys]
    polys_lin = [v for p in polys for v in p]
    codes_lin = [c for cs in codes for c in cs]
    path = Path(polys_lin, codes_lin)
    patch = PathPatch(path, facecolor = 'lightgrey', lw=0)
    # masking the data:
    ax.add_patch(patch)
    return 
# # # # # # # # # # # # #
def outline_region(ax, m, focus_region_states):
    """function finds a polygon defined by a grouping of states and plots 
    an outline of this polygon on map. 
    
    Parameters
    ----------  
    ax : matplotlib.axes._subplots.AxesSubplot
        Axis containing Basemap class instance
    m : mpl_toolkits.basemap.Basemap
        Basemap with specified coordinates and map projection
    focus_region_states : list
        State names (string format) in region 

    Returns
    ----------     
    None    
    """
    from matplotlib.patches import Polygon
    import numpy as np
    import shapely.geometry as sg
    import shapely.ops as so
    from descartes import PolygonPatch        
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants        
    # read in shapefile
    m.readshapefile(pollutants_constants.PATH_SHAPEFILES + 
        'cb_2015_us_state_20m', name = 'states', drawbounds = True)
    state_names = []
    for shape_dict in m.states_info:
        state_names.append(shape_dict['NAME'])
    # dict values are the AQS state codes
    state_fips_code_listing = {'Alaska' : 2, 'Alabama' : 1, 'Arkansas' : 5, 
        'Arizona' : 4, 'California' : 6, 'Colorado' : 8, 'Connecticut' : 9, 
        'District of Columbia' : 11, 'Delaware' : 10, 'Florida' : 12, 
        'Georgia' : 13, 'Hawaii' : 15, 'Iowa' : 19, 'Idaho' : 16, 
        'Illinois' : 17, 'Indiana' : 18, 'Kansas' : 20, 'Kentucky' : 21, 
        'Louisiana' : 22, 'Massachusetts' : 25, 'Maryland' : 24, 'Maine' : 23, 
        'Michigan' : 26, 'Minnesota' : 27, 'Missouri' : 29, 'Mississippi' : 28, 
        'Montana' : 30, 'North Carolina' : 37, 'North Dakota' : 38, 
        'Nebraska' : 31, 'New Hampshire' : 33, 'New Jersey' : 34, 
        'New Mexico' : 35, 'Nevada' : 32, 'New York' : 36, 'Ohio' : 39, 
        'Oklahoma' : 40, 'Oregon' : 41, 'Pennsylvania' : 42, 
        'Rhode Island' : 44, 'South Carolina' : 45, 'South Dakota' : 46, 
        'Tennessee' : 47, 'Texas' : 48, 'Utah' : 49, 'Virginia' : 51, 
        'Vermont' : 50, 'Washington' : 53, 'Wisconsin' : 55, 
        'West Virginia' : 54, 'Wyoming' : 56}
    # iterate through states, if states are in the Northeastern United States, 
    # append shapely.geometry.polygon.Polygon obejcts to list
    patches = [] 
    patches_state = []
    for key, value in state_fips_code_listing.items():
        if key in focus_region_states:
                for info, shape in zip(m.states_info, m.states):
                    if info['NAME'] == key:
                        patches.append(sg.Polygon(shape))
                        patches_state.append(Polygon(np.array(shape), True))  
    # cascaded union can work on a list of shapes, adapted from 
    # https://stackoverflow.com/questions/34475431/plot-unions-of-polygons-in-matplotlib
    neus = so.cascaded_union(patches) 
    ax.add_patch(PolygonPatch(neus, fc = 'None', ec = 'k', alpha = 1.0, 
                 zorder = 11., linewidth = 2.))      
    return 
# # # # # # # # # # # # #    
def find_grid_in_region(m, focus_region_states, gmi_lat, gmi_lon): 
    """given a region comprised of states in variable 'focus_region_states' 
    and a gridded coordinate surface, function finds which grid cells are in 
    the region. The returned array is boolean: 0s mean that the grid cell is 
    not in region and 1s mean that grid cell is. 
    
    Parameters
    ----------  
    m : mpl_toolkits.basemap.Basemap
        Basemap with specified coordinates and map projection
    focus_region_states : list
        State names (string format) in region      
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
        
    Returns
    ----------     
    where_region : numpy.ndarray
        CTM grid where grid cells within region have a value of 1 and grid 
        cells not in region have a value of NaN, [lat, lon]         
    """
    from matplotlib.patches import Polygon
    import numpy as np
    import shapely.geometry as sg
    import shapely       
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants        
    print('finding grid cells in region...')
    # read in shapefile
    m.readshapefile(pollutants_constants.PATH_SHAPEFILES + 
        'cb_2015_us_state_20m', name = 'states', drawbounds = True)
    state_names = []
    for shape_dict in m.states_info:
        state_names.append(shape_dict['NAME'])
    # dict values are the AQS state codes
    state_fips_code_listing = {'Alaska' : 2, 'Alabama' : 1, 'Arkansas' : 5, 
        'Arizona' : 4, 'California' : 6, 'Colorado' : 8, 'Connecticut' : 9, 
        'District of Columbia' : 11, 'Delaware' : 10, 'Florida' : 12, 
        'Georgia' : 13, 'Hawaii' : 15, 'Iowa' : 19, 'Idaho' : 16, 
        'Illinois' : 17, 'Indiana' : 18, 'Kansas' : 20, 'Kentucky' : 21, 
        'Louisiana' : 22, 'Massachusetts' : 25, 'Maryland' : 24, 'Maine' : 23, 
        'Michigan' : 26, 'Minnesota' : 27, 'Missouri' : 29, 'Mississippi' : 28, 
        'Montana' : 30, 'North Carolina' : 37, 'North Dakota' : 38, 
        'Nebraska' : 31, 'New Hampshire' : 33, 'New Jersey' : 34, 
        'New Mexico' : 35, 'Nevada' : 32, 'New York' : 36, 'Ohio' : 39, 
        'Oklahoma' : 40, 'Oregon' : 41, 'Pennsylvania' : 42, 
        'Rhode Island' : 44, 'South Carolina' : 45, 'South Dakota' : 46, 
        'Tennessee' : 47, 'Texas' : 48, 'Utah' : 49, 'Virginia' : 51, 
        'Vermont' : 50, 'Washington' : 53, 'Wisconsin' : 55, 
        'West Virginia' : 54, 'Wyoming' : 56}
    # iterate through states, if states are in region append 
    # shapely.geometry.polygon.Polygon object to list
    patches = [] 
    patches_state = []
    for key, value in state_fips_code_listing.items():
        if key in focus_region_states:
                for info, shape in zip(m.states_info, m.states):
                    if info['NAME'] == key:
                        patches.append(sg.Polygon(shape))
                        patches_state.append(Polygon(np.array(shape), True))  
    where_region = np.empty(shape = np.meshgrid(gmi_lon, gmi_lat)[0].shape)
    where_region[:] = np.nan
    # loop through patches; each patch is a shapely.geometry.polygon.Polygon
    # object
    for patch in patches:
        # for each patch, loop through model grid and check if grid cells
        # are in patch (SEMI-SLOW)
        for i, slat in enumerate(gmi_lat):
            for j, slon in enumerate(gmi_lon): 
                x1, y1 = m(slon, slat)
                point = shapely.geometry.Point(x1, y1)
                if patch.contains(point) == True:
                    where_region[i, j] = 1.
    return where_region
# # # # # # # # # # # # #     
def get_merged_csv(flist, **kwargs):
    """function reads CSV files in the list comprehension loop, this list of
    DataFrames will be passed to the pd.concat() function which will return 
    single concatenated DataFrame
    From https://stackoverflow.com/questions/35973782/reading-multiple-csv-files-concatenate-list-of-file-names-them-into-a-singe-dat
    """
    import pandas as pd
    return pd.concat([pd.read_csv(f, **kwargs) for f in flist], 
                      ignore_index = True)
# # # # # # # # # # # # #    
def o3_givent2mpercentile(o3, t2m_overpass, ptile, gtlt):
    """for a given tempeature percentile function identifies the days which 
    exceed or are less than that threshold and the average O3 on those days 
    in each CTM grid cell. 

    Parameters
    ----------  
    o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, units of volume mixing 
        ratio, [time, lat, lon]
    t2m_overpass : numpy.ndarray
        MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]
    ptile : int
        The percentile of interest
    gtlt : string
        Either 'less' or 'greater.' If 'less' ('greater') function finds days 
        with 2-meter temperatures less than (greater than) percentile threshold 
        
    Returns
    ----------     
    o3_p : numpy.ndarray
        The mean O3 concentrations in each GMI CTM grid cell averaged over days
        which are less than/greater than 2-meter temperature thresholds, 
        units of ppbv, [lat, lon]
    """
    import numpy as np
    o3 = o3 * 1e9
    # find 90th percentile at each grid cell
    p = np.percentile(t2m_overpass, ptile, axis = 0)
    o3_p = np.empty(shape = (p.shape))
    for i in np.arange(0, len(gmi_lat), 1):
        for j in np.arange(0, len(gmi_lon), 1):
            # find days which exceed 90th percentile at each grid cell 
            if gtlt == 'less':
                pdays = np.where(t2m_overpass[:, i, j] < p[i, j])[0]
            if gtlt == 'greater':
                pdays = np.where(t2m_overpass[:, i, j] > p[i, j])[0]            
            # find mean O3 on hot days
            o3_pdays = o3[pdays, i, j]
            o3_p[i, j] = o3_pdays.mean()
    return o3_p
# # # # # # # # # # # # #    
def find_conus_aqsmda8(region_sc):
    """ function opens yearly AQS observations of MDA8 O3 for the measuring 
    period JJA 2000 - 2014. A regional average is calculated for the states 
    whose EPA-designated state codes are contained in the list 'region_sc.' 
    Dataframe column 'Mean Including All Data' is the average (arithmetic mean) 
    value including all data reported in the 8-hour block and column 'Mean
    Excluding All Flagged Data' is the average (arithmetic mean) value 
    excluding all data flagged by the data submitter as having been affected 
    by an exceptional event. For more information, see 
    https://aqs.epa.gov/aqsweb/airdata/FileFormats.html
    
    Parameters
    ---------- 
    region_sc : list
        EPA state codes for states in focus region (from 
        https://aqs.epa.gov/aqsweb/documents/codetables/states_and_counties.csv)
    
    Returns
    ----------      
    ozone_mean : pandas.core.frame.DataFrame
        JJA EPA AQS Data Mart MDA8 O3 for the focus region, units of ppm, 
        [number of years x number of days in summer, 7]    
    ozone_nomean : pandas.core.frame.DataFrame
        JJA EPA AQS Data Mart MDA8 O3 raw data for summertime in the selected 
        region, units of ppm, [number of obs., 12]   
    ozone_mda8 : pandas.core.frame.DataFrame
        MDA8 O3 raw data over CONUS, units of ppm
    """
    import numpy as np
    import pandas as pd
    import os, glob
    # column names for files of MDA8 ozone (parameter code 44201)
    cols = ['State Code', 'County Code', 'Site Num', 
            'Latitude', 'Longitude', 'Date Local',
            'Time Local', 'Observations with Events',
            'Null Observations', 'Mean Including All Data',
            'Mean Excluding All Flagged Data', 
            'Mean Excluding Concurred Flags']
    # specify data types for columns to avoid DtypeWarning being raised
    dtype = {'State Code' : np.str, 'County Code' : np.str, 'Site Num' : np.str, 
             'Latitude' : np.float64, 'Longitude' : np.float64, 
             'Date Local' : np.str, 'Time Local' : np.str, 
             'Observations with Events' : np.int64, 
             'Null Observations' : np.int64, 'Mean Including All Data' : np.float64, 
             'Mean Excluding All Flagged Data' : np.float64, 
             'Mean Excluding Concurred Flags' : np.float64}    
    # wildcard for multiple years of observations
    import time
    # track time to load    
    tic = time.clock()
    dpath = '/Users/ghkerr/phd/stagnation/data/MDA8/'    
    ozone_mda8 = os.path.join(dpath, '8hour_44201_*_JJA.csv')
    ozone_mda8 = get_merged_csv(glob.glob(ozone_mda8), dtype = dtype, 
                                index_col = None, usecols = cols)
    # select only observations from the ozone season, 1 June - 31 Aug
    ozone_mda8['Date Local'] = pd.to_datetime(ozone_mda8['Date Local']) 
    ozone_mda8 = ozone_mda8.loc[ozone_mda8['Date Local'
        ].dt.month.isin([6.0, 7.0, 8.0])]
    ozone_mda8 = ozone_mda8.loc[ozone_mda8['Date Local'
        ].dt.year.isin([2008, 2009, 2010])]
    # select states in focus region
    region_sc = [str(int(x)).zfill(2) for x in region_sc]
    ozone_nomean = ozone_mda8.loc[ozone_mda8['State Code'].isin(region_sc)]
    # find maxima of rolling 8-hour averages of O3 calculated throughout day
    # by grouping by State Code, County Code, and Site Num (these three 
    # constitute a uninque AQS station location) and by date
    ozone_nomean = ozone_nomean.groupby(['State Code', 'County Code', 
        'Site Num', 'Date Local'], as_index = False).max()    
    # daily regional mean O3
    ozone_mean = ozone_nomean.groupby(['Date Local']).mean()
    toc = time.clock()        
    print('AQS loaded in %.1f minutes!' %((toc - tic)/60.))
    return ozone_mean, ozone_nomean, ozone_mda8
# # # # # # # # # # # # #
def afternoon_MERRA(field, times, time_interest, dailyavg):
    """function finds hourly MERRA-2 grids in hours of interest (UTC) and, if 
    desired, computes a daily average over hours of interest.
    
    Parameters
    ----------    
    field : numpy.ndarray
        MERRA-2 gridded field, [time,]
    times : numpy.ndarray
        datetime.datetime objects corresponding to each hourly timestep of 
        MERRA-2 reanalysis [time,]    
    time_interest : list
        Times (UTC) during which MERRA-2 grids are retrieved
    dailyavg : str
        'yes' or other

    Returns
    ----------        
    field_intersect : numpy.ndarray
        MERRA-2 gridden field from hours in 'time_interest', [time,]
    """
    import numpy as np
    # find hours which fall into times of interest
    intersect = np.in1d([x.hour for x in times], time_interest)
    intersect = np.where(intersect == True)[0]
    field_intersect = field[intersect]
    # find daily average over hours in 'time_interest'
    if dailyavg == 'yes':
        # for MERRA-2 field
        field_intersect = np.reshape(field_intersect, 
            [int(field_intersect.shape[0]/len(time_interest)), 
            len(time_interest)])
        field_intersect = np.nanmean(field_intersect, axis = 1)
    return field_intersect
# # # # # # # # # # # # #
def aqs_mda8_r_do3dt2m(ozone_nomean_mda8, lat, lon): 
    """functions uses maximum daily 8-hour average O3 from AQS over the 
    continential U.S. and finds (1) the Pearson correlation coefficient (r) and 
    (2) O3-T2m sensitivity in each MERRA-2 grid cells with AQS observations 
    for JJA 2008-2010. 
    
    Parameters
    ----------    
    ozone_nomean_mda8 : pandas.core.frame.DataFrame
        JJA EPA AQS Data Mart MDA8 O3 raw data for summertime in the selected 
        region, units of ppm, [number of obs., 12]       
    lat : numpy.ndarray
        MERRA-2 latitude coordinates, units of degrees north, [lat,]
    lon : numpy.ndarray
        MERRA-2 longitude coordinates, units of degrees east, [lon,]   

    Returns
    ----------     
    r : numpy.ndarray
        Pearson correlation coefficients between MERRA-2 2 meter temperatures
        and AQS ozone at each MERRA-2 grid cell containing AQS O3 observations, 
        [lat, lon]
    do3dt : numpy.ndarray
        The O3-T2m sensitivity (slope of linear regression) 
        at each MERRA-2 grid cell containing AQS O3 observations, [lat, lon]    
    """
    import numpy as np
    import datetime
    import pytz
    import timezonefinder
    # arrays will be filled with correlation coefficients and O3-T2m 
    # sensitivities
    r = np.empty((lat.shape[0], lon.shape[0]))
    do3dt = np.empty((lat.shape[0], lon.shape[0]))
    # local times to fetch MERRA-2 fields
    afternoon_times = [11, 12, 13, 14, 15, 16]
    # make Date Local DataFrame index
    ozone_nomean_mda8.index = pd.DatetimeIndex(ozone_nomean_mda8['Date Local'])
    # define date range
    date_idx = []
    for year in np.arange(2008, 2011, 1):
        date_idx.append(pd.date_range('06-01-%d' %year, '08-31-%d' %year))
    # aggregate over years
    date_idx = np.hstack(date_idx)
    # MERRA-2 resolution
    latres = np.diff(lat).mean()
    lonres = np.diff(lon).mean()
    # loop through MERRA-2 lat/lon
    for i, ilat in enumerate(lat):
        for j, ilon in enumerate(lon):
            ulat = ilat + latres/2.
            llat = ilat - latres/2.
            llon = ilon - lonres/2.
            rlon = ilon + lonres/2.
            # find AQS observations in grid cell
            ozone_incell = ozone_nomean_mda8[(ozone_nomean_mda8
                ['Latitude'].between(llat, ulat, inclusive = True)) & 
                (ozone_nomean_mda8['Longitude'].between(llon, rlon, inclusive = True))]
            if ozone_incell.empty == True:
                r[i, j] = np.nan
                do3dt[i, j] = np.nan
            else: 
                # if > 1 AQS stations exist, average over stations
                ozone_incell = ozone_incell.groupby([ozone_incell.index.date]).mean()
                # fill missing value
                ozone_incell = ozone_incell.reindex(date_idx, fill_value = np.nan)
                # find UTC offset in grid node, from 
                # stackoverflow.com/questions/15742045/getting-time-zone-from-lat-long-coordinates
                tf = timezonefinder.TimezoneFinder()
                timezone_str = tf.certain_timezone_at(lat = ilat, lng = ilon)
                # stackoverflow.com/questions/27697117/get-the-offset-for-a-specific-time-and-timezone-in-python
                tz = pytz.timezone(timezone_str)
                # dummy summer date (to pick up on DST)
                dt = datetime.datetime.strptime('2010-06-01', '%Y-%m-%d')
                offset = int(tz.utcoffset(dt).total_seconds() / 3600.0)
                # since we are interested in 11am - 4pm local time, find these times 
                # in UTC for a given time zone
                times_incell = [time - (offset) for time in afternoon_times]
#                # # # # to check right temperatures were selected, n.b. the x 
#                # positions (i.e., [t+(24*2) for t in times_incell]) should match 
#                # up with the true/false in variable intersect in function 
#                # afternoon_MERRA
#                plt.plot(t2m[24*0:24*4, 21, 11], '-k')
#                plt.plot([t+(24*0) for t in times_incell], t2m[24*1-6:24*1, i, j], '-ro')
#                plt.plot([t+(24*1) for t in times_incell], t2m[24*2-6:24*2, i, j], '-ro')
#                plt.plot([t+(24*2) for t in times_incell], t2m[24*3-6:24*3, i, j], '-ro')
#                plt.plot([t+(24*3) for t in times_incell], t2m[24*4-6:24*4, i, j], '-ro')
#                plt.show()        
                # find MERRA-2 2-meter temperatures in grid cell during the afternoon
                t2m_incell = t2m[:, i, j]        
                t2m_afternoon = afternoon_MERRA(t2m_incell, times_all, 
                                                times_incell, 'yes')
                # mask for missing values
                mask = ~np.isnan(t2m_afternoon) & \
                ~np.isnan(ozone_incell['Mean Including All Data'].values)
                # calculate Pearson correlation coefficient and slope of linear 
                # regression
                r_incell = np.corrcoef(t2m_afternoon[mask], 
                    ozone_incell['Mean Including All Data'].values[mask])[0, 1]
                r[i, j] = r_incell
                do3dt_incell = np.polyfit(t2m_afternoon[mask], 
                    ozone_incell['Mean Including All Data'].values[mask] * 1000., 1)[0]
                do3dt[i, j] = do3dt_incell
    return r, do3dt
# # # # # # # # # # # # #
def find_conus_castnet(years): 
    """function opens CASTNet observations for the specified measuring period 
    and selects summer (JJA) observations 
    
    Parameters
    ----------    
    years : list
        Years of interest       

    Returns
    ----------     
    castnet_all : pandas.core.frame.DataFrame
        Station-based hourly CASTNet O3 concentrations, [no. obs., 8]
    """    
    import time
    import numpy as np
    import pandas as pd
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/GMI/')
    import commensurability
    # track time to load    
    tic = time.clock()
    # open siting information (lat/lon) for CASTNet stations
    csites = commensurability.open_castnetsiting()
    castnet_all = pd.DataFrame([])
    for year in years:
        # open CASTNET O3 observations for year of interest
        castnet = pd.read_csv((pollutants_constants.PATH_CASTNET + 
                               '%s/ozone_%s.csv' %(year, year)), header = 0,
                               low_memory = False)
        # add columns names to CASTNET dataframe
        castnet_cols = ['SITE_ID',
                        'DATE_TIME',
                        'OZONE',
                        'OZONE_F',                    
                        'QA_CODE',
                        'UPDATE_DATE']
        castnet = pd.DataFrame(castnet, columns = castnet_cols)
        # convert 'DATE_TIME' column from pandas.core.series.Series to a DateTime 
        # object 
        castnet['DATE_TIME'] = pd.to_datetime(castnet['DATE_TIME'])
        # only select summer months (JJA)
        castnet = castnet.loc[castnet['DATE_TIME'].dt.month.isin([6, 7, 8])]
        # remove the last three letters of SITE_ID column for direct comparison 
        # with Cooper/GMI station names 
        castnet['SITE_ID'] = castnet['SITE_ID'].astype(str).str[:-3].astype(np.str)
        castnet['DATE_TIME'] = castnet['DATE_TIME'].astype(str).str[:-6]    
        # append to multi-year DataFrame 
        castnet_all = castnet_all.append(castnet)
    # loop through each CASTNet observation and find latitude/longitude 
    # of station    
    lats_all, lons_all = [], []
    for csite in castnet_all['SITE_ID']:
        csite_atsite = csites.loc[csites['SITE_ID'].isin([csite])]
        lons_all.append(csite_atsite['LONGITUDE'].values[0])
        lats_all.append(csite_atsite['LATITUDE'].values[0])
    # add latitudes and longitudes to CASTNet data
    castnet_all['Latitude'] = lats_all
    castnet_all['Longitude'] = lons_all
    toc = time.clock()
    print('CASTNet data loaded in %.1f minutes!' %((toc - tic)/60.))
    return castnet_all
# # # # # # # # # # # # #     
def castnet_r_do3d2t(castnet, t2m, lat, lon, times_all):
    """function groups CASTNet observations by station, finds the afternoon 
    (1100-1600 local time) average ozone concentration and then finds the 
    nearest MERRA-2 grid node. The Pearson correlation coefficient and 
    O3-T2m sensitivity at each CASTNet site. n.b. 11/11/2018 variable 
    'afternoon_times' changed to [13, 14] to reflect that the overpass2 time 
    is not exactly at 1300 hours local time but sometime between 1300 and 1400
    hours. 

    Parameters
    ----------           
    castnet : pandas.core.frame.DataFrame
        Station-based hourly CASTNet O3 concentrations, [no. obs., 8]
    t2m : numpy.ndarray
        MERRA-2 hourly 2-meter air temperature, units of K, [time, lat, lon]
    lat : numpy.ndarray
        MERRA-2 latitude coordinates, units of degrees north, [lat,]
    lon : numpy.ndarray
        MERRA-2 longitude coordinates, units of degrees east, [lon,]   
    times_all : numpy.ndarray
        datetime.datetime objects corresponding to each hourly timestep of the 
        reanalysis, [time,]        

    Returns
    ----------   
    r_all : list
        Pearson correlation coefficients between MERRA-2 2 meter temperatures
        and CASTNet ozone at each CASTNet site, [no. sites,]
    do3dt2m_all : list
        The O3-T2m sensitivity (slope of linear regression fit between O3 
        versus T2m) at each CASTNet site, units of ppbv K^-1, [no. sites,]
    lat_all : list
        Latitudes of CASTNet sites, units of degrees north, [no. sites,]
    lon_all : list
        Longitudes of CASTNet sites, units of degrees west, [no. sites,]
    t_all : list
        Daily 1300 hours local time MERRA-2 2-meter temperatures at each 
        CASTNet site, units of K, [no. sites,]
    o3_all : list
        Daily 1300 hours local time O3 at each CASTNet site, units of ppbv, 
        [no. sites,]
    sites_all : list
        CASTNet site names, [no. sites,]
    """
    import numpy as np
    import pandas as pd
    import datetime
    import pytz
    import timezonefinder
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    afternoon_times = [13, 14]
    castnet['DATE_TIME'] = pd.to_datetime(castnet['DATE_TIME'])
    # define date range
    date_idx = []
    for year in np.arange(2008, 2011, 1):
        date_idx.append(pd.date_range('06-01-%d' %year, '08-31-%d' %year))
    # aggregate over years
    date_idx = np.hstack(date_idx)
    # lists to be filled with data from CASTNet sites
    lat_all, lon_all, r_all, do3dt2m_all, sites_all = [], [], [], [], []
    t_all, o3_all = [], []
    for siteid in np.unique(castnet['SITE_ID'].values):
        castnet_atsite = castnet.loc[castnet['SITE_ID'] == siteid]
        # only proceed if there are observations at site 
        if castnet_atsite['OZONE'].isnull().all() != True:
            lat_atsite = np.nanmean(castnet_atsite['Latitude'].values)
            lon_atsite = np.nanmean(castnet_atsite['Longitude'].values)
            # find CASTNet observations in afternoon 
            castnet_atsite = castnet_atsite.loc[
                    castnet_atsite['DATE_TIME'].dt.hour.isin(afternoon_times)]
            # sort by time
            castnet_atsite = castnet_atsite.sort_values(by = 'DATE_TIME')    
            # find daily average
            castnet_atsite = castnet_atsite.groupby(castnet_atsite['DATE_TIME'].dt.date).mean()
            # add missing values             
            if castnet_atsite.shape[0] != 276:
                castnet_atsite = castnet_atsite.reindex(date_idx, 
                    fill_value = np.nan)
            # find closest MERRA-2 grid cell
            merra_lat_idx = geo_idx(lat_atsite, lat)
            merra_lon_idx = geo_idx(lon_atsite, lon)
            # only consider if MERRA-2 gridcell is nearby CASTNet site
            if not (merra_lat_idx or merra_lon_idx):
                pass
            else:
                # find UTC offset in grid node
                tf = timezonefinder.TimezoneFinder()
                timezone_str = tf.certain_timezone_at(lat = lat_atsite, 
                                                      lng = lon_atsite)
                tz = pytz.timezone(timezone_str)
                dt = datetime.datetime.strptime('2010-06-01', '%Y-%m-%d')
                offset = int(tz.utcoffset(dt).total_seconds() / 3600.0)
                times_incell = [time - (offset) for time in afternoon_times]    
                # find MERRA-2 2-meter temperatures in grid cell during the afternoon
                t2m_atsite = t2m[:, merra_lat_idx, merra_lon_idx]        
                t2m_afternoon = afternoon_MERRA(t2m_atsite, times_all, 
                                                times_incell, 'yes')
                # add daily values at each sites to lists
                t_all.append(t2m_afternoon)
                o3_all.append(castnet_atsite['OZONE'].values)
                sites_all.append(siteid)
                # mask for missing values
                mask = ~np.isnan(t2m_afternoon) & \
                ~np.isnan(castnet_atsite['OZONE'].values)
                # calculate Pearson correlation coefficient and slope of linear 
                # regression
                r_atsite = np.corrcoef(t2m_afternoon[mask], 
                    castnet_atsite['OZONE'].values[mask])[0, 1]
                do3dt_atsite = np.polyfit(t2m_afternoon[mask], 
                    castnet_atsite['OZONE'].values[mask], 1)[0]
                # for sites with reliable data, save off correlation 
                # coefficients, O3-T2m sensitivities, and station latitudes and 
                # longitudes
                lat_all.append(lat_atsite)
                lon_all.append(lon_atsite)
                r_all.append(r_atsite)
                do3dt2m_all.append(do3dt_atsite)
#    # OPTIONAL: plot location and names of CASTNet sites
#    fig = plt.figure(figsize = (9, 7))
#    m = Basemap(projection = 'merc', llcrnrlon = -126., 
#                llcrnrlat = 24., urcrnrlon = -66.3, 
#                urcrnrlat = 50., resolution = 'h', area_thresh = 1000)
#    ax = plt.subplot2grid((1, 1), (0, 0))
#    x, y = m(lon_castnet, lat_castnet)
#    m.scatter(x, y, color = 'k', s = 4)
#    m.drawstates(color = 'k', linewidth = 0.5)
#    m.drawcountries(color = 'k', linewidth = 1.0)
#    m.drawcoastlines(color = 'k', linewidth = 1.0)    
#    for i, txt in enumerate(sites_castnet):
#        ax.annotate(txt, (x[i] + 10000, y[i]), fontsize = 8, color = 'r')
#    plt.savefig('/Users/ghkerr/Desktop/CASTNETmap.eps', dpi = 300)                
    return r_all, do3dt2m_all, lat_all, lon_all, t_all, o3_all, sites_all
# # # # # # # # # # # # #    
def calculate_gmi_r_do3dt2m(merra_lat, merra_lon, gmi_lat, gmi_lon, t2m, o3, 
    ps):
    """using hourly MERRA-2 temperature fields and CTM O3 fields at overpass
    time, function calculates the ozone-temperature sensitivity (dO3/dT2m) and 
    the O3-T2m correlation coefficient. Function first determines the UTC 
    offset in each grid cell (inefficient) and then finds the MERRA-2 1pm (1300 
    hours) local time 2-meter temperature on each day in each grid cell given 
    the UTC offsets. The MERRA-2 temperature grid at overpass time is then 
    interpolated to the resolution of the CTM and the sensitivity is found by 
    computing a least-squares fit. The correlation is found with the Pearson 
    product-moment correlation coefficient.
    
    Parameters
    ----------  
    merra_lat : numpy.ndarray
        MERRA-2 latitude coordinates, units of degrees north, [lat,]
    merra_lon : numpy.ndarray
        MERRA-2 longitude coordinates, units of degrees east, [lon,]       
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
    t2m : numpy.ndarray
        MERRA-2 hourly 2-meter air temperature, units of K, [time, lat, lon]
    o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, units of volume mixing 
        ratio, [time, lat, lon]        
    ps : numpy.ndarray
        MERRA-2 hourly sea level pressure, units of Pa, [time, lat, lon]        

    Returns
    ----------     
    do3dt2m : numpy.ndarray    
        The O3-T2m sensitivity at each GMI grid cell, units of ppbv K^-1,
        [lat, lon]         
    tls : numpy.ndarray
        The total least squares (or "orthogonal distance regression") of 
        O3 vs T2m. This is the same as variable "do3dt2m" but assumes that
        errors are present in both the independent (T2m) and dependent (O3)
        variables, units of ppbv K^-1, [lat, lon]
    r : numpy.ndarray
        The Pearson product-moment correlation coefficient calculated between 
        MERRA-2 2-meter temperatures and ozone at each GMI grid cell [lat, lon]
    t2m_atoverpass_interp : numpy.ndarray
        MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]
    ps_atoverpass_interp : numpy.ndarray
        MERRA-2 sea level pressure at overpass2 time interpolated to the 
        resolution of the CTM, units of Pa, [time, lat, lon]
    p : numpy.ndarray
        The p-values for the correlation coefficient calculated between 
        MERRA-2 2-meter temperatures and ozone at each GMI grid cell; the 
        p-values roughly indicates the probability of an uncorrelated system 
        producing datasets that have a Pearson correlation at least as extreme 
        as the one computed from these datasets, [lat, lon]        
    """
    import numpy as np
    import scipy.interpolate
    import datetime
    import pytz
    import timezonefinder
    import scipy.odr    
    import scipy.stats
    # # # #
    def fit_func(p, t):
        """first order linear regression for calculating total least squares
        """
        return p[0] * t + p[1]    
    # # # #
    # convert O3 units 
    o3 = o3 * 1e9
    # MERRA-2 output is in UTC, GMI overpass2 files are in local time, so we 
    # need to find the UTC offset at each MERRA-2 grid cell 
    print('calculating overpass time in UTC...')
    merra_tz = np.empty(shape = t2m[0].shape)
    merra_tz[:] = np.nan
    for i, ilat in enumerate(merra_lat):
        for j, jlon in enumerate(merra_lon):        
            tf = timezonefinder.TimezoneFinder()
            timezone_str = tf.closest_timezone_at(lat = ilat, lng = jlon)
            if not timezone_str:
                timezone_str = tf.closest_timezone_at(lat = ilat, 
                                                      lng = jlon, 
                                                      delta_degree = 11)
            tz = pytz.timezone(timezone_str)
            dt = datetime.datetime.strptime('2010-06-01', '%Y-%m-%d')
            offset = int(tz.utcoffset(dt).total_seconds() / 3600.0)
            merra_tz[i, j] = offset
    # now find the time of overpass (for simplicity 1300 hours local time) in 
    # UTC in each MERRA-2 grid cell and pull all temperatures in this cell at 
    # the overpass time         
    print('finding fields at overpass time...')        
    t2m_atoverpass = np.empty(shape = (int(t2m.shape[0]/24), 
                                       t2m.shape[1], t2m.shape[2]))
    t2m_atoverpass[:] = np.nan
    ps_atoverpass = np.empty(shape = (int(ps.shape[0]/24),
                                      ps.shape[1], ps.shape[2]))
    ps_atoverpass[:] = np.nan
    for i, ilat in enumerate(merra_lat):
        for j, jlon in enumerate(merra_lon):
            overpassutc = np.abs(merra_tz[i, j]) + 13
            t2m_atoverpass[:, i, j] = t2m[np.arange(int(overpassutc), 
                          len(t2m), 24), i, j]  
            ps_atoverpass[:, i, j] = ps[np.arange(int(overpassutc), 
                          len(ps), 24), i, j]                  
    # interpolate from resolution of MERRA-2 to resolution of CTM; n.b. 
    # here X, Y are the old grid dimensions and XI, YI are the new grid 
    print('interpolating fields to CTM resolution...')
    X, Y = np.meshgrid(merra_lon, merra_lat)
    XI, YI = np.meshgrid(gmi_lon, gmi_lat)
    t2m_atoverpass_interp = np.empty(o3.shape)
    ps_atoverpass_interp = np.empty(o3.shape)
    for z in np.arange(0, len(t2m_atoverpass), 1):    
        t2m_atoverpass_interp[z] = scipy.interpolate.griddata(
                (X.flatten(), Y.flatten()), t2m_atoverpass[z].flatten(), 
                (XI, YI), method = 'cubic')
        ps_atoverpass_interp[z] = scipy.interpolate.griddata(
                (X.flatten(), Y.flatten()), ps_atoverpass[z].flatten(), 
                (XI, YI), method = 'cubic')        
    # calculate the O3-temperature sensitivity at each CTM grid cell 
    print('calculating dO3-dT2m...')
    do3dt2m = np.empty(shape = t2m_atoverpass_interp.shape[1:])
    do3dt2m[:] = np.nan
    for i, ilat in enumerate(gmi_lat):
        for j, jlon in enumerate(gmi_lon):
            do3dt2m[i, j] = np.polyfit(t2m_atoverpass_interp[:, i, j], 
                                       o3[:, i, j], 1)[0]
    print('calculating total least squares...')            
    # total least squares/orthogonal distance regression
    Model = scipy.odr.Model(fit_func)
    tls = np.empty(shape = t2m_atoverpass_interp.shape[1:])
    tls[:] = np.nan
    for i, ilat in enumerate(gmi_lat):
        for j, jlon in enumerate(gmi_lon):
            odr_ij = scipy.odr.RealData(t2m_atoverpass_interp[:, i, j], 
                                        o3[:, i, j])
            odr_ij = scipy.odr.ODR(odr_ij, Model, 
                [np.polyfit(t2m_atoverpass_interp[:, i, j], o3[:, i, j], 1)[1], 
                 np.polyfit(t2m_atoverpass_interp[:, i, j], o3[:, i, j], 1)[0]],
                 maxit = 10000)
            output_ij = odr_ij.run()            
            beta_ij = output_ij.beta            
            tls[i, j] = beta_ij[0]
    # calculate the O3-temperature sensitivity at each CTM grid cell 
    print('calculating correlation coefficients...')
    r = np.empty(shape = t2m_atoverpass_interp.shape[1:])
    r[:] = np.nan
    p = np.empty(shape = t2m_atoverpass_interp.shape[1:])
    p[:] = np.nan    
    for i, ilat in enumerate(gmi_lat):
        for j, jlon in enumerate(gmi_lon):
            r[i, j] = np.corrcoef(t2m_atoverpass_interp[:, i, j], 
                                       o3[:, i, j], 1)[0, 1]
            p[i, j] = np.corrcoef(t2m_atoverpass_interp[:, i, j], 
                                       o3[:, i, j], 1)[0, 1]
    return do3dt2m, tls, r, t2m_atoverpass_interp, ps_atoverpass_interp, p
# # # # # # # # # # # # #
def calculate_castnet_r_do3dt2m_regionmean(t_castnet, o3_castnet, 
    sites_castnet, region_castnet):
    """function finds daily regionally-averaged CASTNet O3 and co-located 
    2-meter temperatures within a given region specified by the CASTNet site 
    IDs in variable 'sites_region.' The Pearson correlation coefficient and the 
    O3-T sensitivity are calculated. On 11/11/2018 function amended to 
    calculate regionally-averaged statistics (i.e., mean, standard deviation, 
    90th/10th percentiles, etc); previously these had been calculated 
    incorrectly, for example, by averaging over all O3-climate penalties in a
    region and reporting this as the mean. The regionally-averaged penalty 
    should be calculated using regionally-averagd O3 and T

    Parameters
    ----------  
    t_castnet : list
        Daily 1300 hours local time MERRA-2 2-meter temperatures at each 
        CASTNet site, units of K, [no. sites,]
    o3_castnet : list
        Daily 1300 hours local time O3 at each CASTNet site, units of ppbv, 
        [no. sites,]
    sites_castnet : list
        CASTNet site names, [no. sites,]
    region_castnet : list
        Site IDs for CASTNet sites in region, [no. sites in region,]
     
    Returns
    ----------     
    do3dt2m : numpy.float64
        The O3 - 2-meter temperatures sensitivity calculated with 
        regionally-averaged O3 from CASTNet and co-located 2-meter temperatures
    r : numpy.float64 
        The Pearson correlation coefficient between regionally-averaged O3 from 
        CASTNet and co-located 2-meter temperatures
    t2m_region : numpy.ndarray
        Daily regionally-averaged 2-meter temperatures co-located with CASTNet 
        stations, [time,]
    o3_region : numpy.ndarray
        Daily regionally-averaged O3 from CASTNet, [time,]          
    """
    import scipy.odr
    import numpy as np
    # # # #
    def fit_func(p, t):
        """first order linear regression for calculating total least squares
        """
        return p[0] * t + p[1]    
    # # # #
    where_region = np.in1d(sites_castnet, region_castnet)
    where_region = np.where(where_region == True)[0]
    # O3, T2m at CASTNet sites in region
    o3_region = np.array(o3_castnet)[where_region]
    t2m_region = np.array(t_castnet)[where_region]
    # find daily regional average
    o3_region = np.nanmean(o3_region, axis = 0)
    t2m_region = np.mean(t2m_region, axis = 0)
    r = np.corrcoef(t2m_region, o3_region)[0, 1]
    do3dt2m = np.polyfit(t2m_region, o3_region, 1)[0]
    Model = scipy.odr.Model(fit_func)
    odr = scipy.odr.RealData(t2m_region, o3_region)
    odr = scipy.odr.ODR(odr, Model, 
                [np.polyfit(t2m_region, o3_region, 1)[1],
                 np.polyfit(t2m_region, o3_region, 1)[0]], maxit = 10000)
    output = odr.run()            
    beta = output.beta    
    # identify days with temperature extremes 
    p90 = np.percentile(t2m_region, 90, axis = 0)
    p90 = np.where(t2m_region > p90)[0]
    p10 = np.percentile(t2m_region, 10, axis = 0)    
    p10 = np.where(t2m_region < p10)[0]
    o3_med = np.nanmedian(o3_region)    
    # find O3 on days with temperature extremes 
    o3_p90 = np.mean(o3_region[p90])
    o3_p10 = np.mean(o3_region[p10])
    # print out the following values 
    print('For CASTNet...')
    print('Mean = %.4f ppbv' %(np.nanmean(o3_region)))
    print('Standard deviation = %.4f ppbv' 
          %(np.nanstd(o3_region)))
    print('P10 - median = %.4f ppbv' 
          %(np.percentile(o3_region, 10) - o3_med))  
    print('P90 - median = %.4f ppbv' 
          %(np.percentile(o3_region, 90) - o3_med))
    print('O3-climate penalty with OLS = %.4f ppbv/K' %do3dt2m)
    print('O3-climate penalty with total least squares = %.4f ppbv/K' %beta[0])        
    print('O3 difference on hot - median days = %.4f ppbv' %(o3_p90 - o3_med))
    print('O3 difference on cold - median days = %.4f ppbv' %(o3_p10 - o3_med))    
    return do3dt2m, r, t2m_region, o3_region
# # # # # # # # # # # # #    
def calculate_gmi_r_do3dt2m_regionmean(t2m_overpass, o3, region, case):
    """function averages daily 2-meter temperature and surface-level O3 fields 
    into average values within a region and thereafter calculates the Pearson
    correlation coefficient and O3-T2m sensitivity.
    
    Parameters
    ----------  
    t2m_atoverpass : numpy.ndarry
        The MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]    
    o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, units of volume mixing 
        ratio, [time, lat, lon]   
    region : numpy.ndarray
        CTM grid where grid cells within region have a value of 1 and grid 
        cells not in region have a value of NaN, [lat, lon]   
    case : str
        Model simulation name
        
    Returns
    ----------     
    do3dt2m : numpy.float64
        The O3-T2m sensitivity calculated with regionally-averaged O3 and 
        T2m
    r : numpy.float64 
        The Pearson correlation coefficient between O3 and T2m calculated with 
        regionally-averaged O3 and T2m
    t2m_alldays : numpy.ndarray
        Daily regionally-averaged T2m, [time,]
    o3_alldays : numpy.ndarray
        Daily regionally-averaged O3, [time,]    
    """
    import scipy.odr
    import numpy as np
    # # # #
    def fit_func(p, t):
        """first order linear regression for calculating total least squares
        """
        return p[0] * t + p[1]    
    # # # #
    o3 = o3 * 1e9
    t2m_alldays, o3_alldays = [], []    
    for t2m_day, o3_day in zip(t2m_overpass, o3): 
        t2m_day = region * t2m_day
        o3_day = region * o3_day
        # find regionally-averaged temperatures, O3 on day in region 
        t2m_alldays.append(np.nanmean(t2m_day))
        o3_alldays.append(np.nanmean(o3_day))
    # calculate correlation coefficient
    r = np.corrcoef(t2m_alldays, o3_alldays)[0, 1]
    # calculate climate penalty with OLS 
    do3dt2m = np.polyfit(t2m_alldays, o3_alldays, 1)[0]
    # calculate climate penalty with total least squares (ODR)
    Model = scipy.odr.Model(fit_func)
    odr = scipy.odr.RealData(t2m_alldays, o3_alldays)
    odr = scipy.odr.ODR(odr, Model, 
                [np.polyfit(t2m_alldays, o3_alldays, 1)[1],
                 np.polyfit(t2m_alldays, o3_alldays, 1)[0]], maxit = 10000)
    output = odr.run()            
    beta = output.beta
    # identify days with temperature extremes 
    t_p90 = np.percentile(t2m_alldays, 90, axis = 0)
    t_p90 = np.where(t2m_alldays > t_p90)[0]
    t_p10 = np.percentile(t2m_alldays, 10, axis = 0)    
    t_p10 = np.where(t2m_alldays < t_p10)[0]
    o3_med = np.nanmedian(o3_alldays)
    # find O3 on days with temperature extremes 
    o3_p90 = np.mean(np.array(o3_alldays)[t_p90])
    o3_p10 = np.mean(np.array(o3_alldays)[t_p10])
    # print out the following values 
    print('For %s simulation...' %case)
    print('Mean = %.4f ppbv' 
          %(np.nanmean(o3 * region)))
    print('Standard deviation = %.4f ppbv' 
          %(np.std(np.nanmean(o3 * region, axis = tuple((1, 2))))))
    print('P10 - median = %.4f ppbv' 
          %(np.percentile(np.nanmean((o3 * region), axis = tuple((1, 2))), 10) - 
          o3_med))
    print('P90 - median = %.4f ppbv' 
          %(np.percentile(np.nanmean((o3 * region), axis = tuple((1, 2))), 90) - 
          o3_med))
    print('O3-climate penalty with OLS = %.4f ppbv/K' %do3dt2m)
    print('O3-climate penalty with total least squares = %.4f ppbv/K' %beta[0])    
    print('O3 difference on hot - median days = %.4f ppbv' %(o3_p90 - o3_med))
    print('O3 difference on cold - median days = %.4f ppbv' %(o3_p10 - o3_med))        
    return do3dt2m, r, np.array(t2m_alldays), np.array(o3_alldays)
# # # # # # # # # # # # #    
def map_ro3t2m_do3dt2m_conus(lat_castnet, lon_castnet, r_castnet, 
    do3dt2m_castnet, lat, lon, r, do3dt2m):
    """function plots maps of (1) the correlation of MDA8 AQS O3 and 1100-1600 
    localtime mean CASTNet O3 and MERRA-2 2 meter temperature and (2) the
    ozone-temperature sensitivity for the same datasets. 

    Parameters
    ----------    
    lat_castnet : list
        CASTNet latitude coordinates, units of degrees north, [no. CASTNet
        stations,]
    lon_castnet : list
        CASTNet longitude coordinates, units of degrees north, [no. CASTNet
        stations,]    
    r_castnet : list
        Pearson correlation coefficients between MERRA-2 2 meter temperatures
        and CASTNet ozone at each CASTNet site, [no. CASTNet stations,]
    do3dt2m_castnet : list
        The O3-T2m sensitivity (slope of linear regression) at each CASTNet
        site, [no. CASTNet stations,]        
    lat : numpy.ndarray
        MERRA-2 latitude coordinates, units of degrees north, [lat,]
    lon : numpy.ndarray
        MERRA-2 longitude coordinates, units of degrees east, [lon,]   
    r : numpy.ndarray
        Pearson correlation coefficients between MERRA-2 2 meter temperatures
        and AQS ozone at each MERRA-2 grid cell containing AQS O3 observations, 
        [lat, lon]
    do3dt2m : numpy.ndarray
        The O3-T2m sensitivity (slope of linear regression) 
        at each MERRA-2 grid cell containing AQS O3 observations, [lat, lon]            

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
    # continental U.S. focus region map 
    llcrnrlon = -130.
    llcrnrlat = 24.8
    urcrnrlon = -66.3
    urcrnrlat = 50.
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
        llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat, 
        resolution = 'h', area_thresh = 1000)
    # map for O3-T2m correlation
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    vmin = 0.0; vmax = 0.75
    cmap = plt.get_cmap('gist_earth', 6)
    norm = Normalize(vmin = vmin, vmax = vmax)
    # color mapper to covert values to colors
    mapper = ScalarMappable(norm = norm, cmap = cmap)
    # add CASTNet data as scatterpoints
    x, y = m(lon_castnet, lat_castnet)
    m.scatter(x, y, c = np.array(r_castnet), s = 30, cmap = cmap, 
              vmin = vmin, vmax = vmax, zorder = 30, edgecolor = 'k')
    # MERRA-2 resolution 
    rlon = np.diff(lon).mean()
    rlat = np.diff(lat).mean()
    # loop through GMI latitude and longitude coordinatesion
    for i, lat_gb in enumerate(lat):
        for j, lon_gb in enumerate(lon):
            lon_left = lon_gb - (rlon/2.)
            lon_right = lon_gb + (rlon/2.)
            lat_top = lat_gb + (rlat/2.)
            lat_bottom = lat_gb - (rlat/2.) 
            x1, y1 = m(lon_right, lat_bottom)
            x2, y2 = m(lon_right, lat_top)
            x3, y3 = m(lon_left, lat_top)
            x4, y4 = m(lon_left, lat_bottom)
            if np.isnan(r[i, j]) != True:
                color = rgb2hex(mapper.to_rgba(r[i, j]))
                p = Polygon([(x1, y1), (x2, y2), (x3, y3), (x4, y4)], 
                             facecolor = color,
                             edgecolor = 'none', linewidth = 0.0, zorder = 5, 
                             alpha = 1.0) 
                ax.add_patch(p) 
    m.drawmapboundary(fill_color = '#EBF4FA')
    m.fillcontinents(color = '#CCCCCC', lake_color = '#EBF4FA')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    # colorbar 
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_ticks(np.linspace(vmin, vmax, 7))
    cb.set_label(label = '$r_{\mathregular{O_{3}-T_{2m}}}$', size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/map_ro3t2m_conus.eps', dpi = 300)
    # map for O3-T2m sensitivity
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
                llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, 
                urcrnrlat = urcrnrlat, resolution = 'h', area_thresh = 1000)    
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    vmin = 0.0; vmax = 3.
    cmap = plt.get_cmap('gist_earth', 6)
    norm = Normalize(vmin = vmin, vmax = vmax)
    mapper = ScalarMappable(norm = norm, cmap = cmap)
    m.scatter(x, y, c = np.array(do3dt2m_castnet), s = 30, cmap = cmap, 
              vmin = vmin, vmax = vmax, zorder = 30, edgecolor = 'k')
    for i, lat_gb in enumerate(lat):
        for j, lon_gb in enumerate(lon):
            lon_left = lon_gb - (rlon/2.)
            lon_right = lon_gb + (rlon/2.)
            lat_top = lat_gb + (rlat/2.)
            lat_bottom = lat_gb - (rlat/2.) 
            x1, y1 = m(lon_right, lat_bottom)
            x2, y2 = m(lon_right, lat_top)
            x3, y3 = m(lon_left, lat_top)
            x4, y4 = m(lon_left, lat_bottom)
            if np.isnan(do3dt2m[i, j]) != True:
                color = rgb2hex(mapper.to_rgba(do3dt2m[i, j]))
                p = Polygon([(x1, y1), (x2, y2), (x3, y3), (x4, y4)], 
                             facecolor = color,
                             edgecolor = 'none', linewidth = 0.0, zorder = 5, 
                             alpha = 1.0) 
                ax.add_patch(p) 
    m.drawmapboundary(fill_color = '#EBF4FA')
    m.fillcontinents(color = '#CCCCCC', lake_color = '#EBF4FA')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_ticks(np.linspace(vmin, vmax, 7))
    cb.set_label(label = '$\partial$O$_{3}$ ' + 
                 '$\partial$T$_{\mathregular{2 m}}^{\mathregular{-1}}$', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/map_do3dt2m_conus.eps', dpi = 300)
    return
# # # # # # # # # # # # # 
def map_ro3t2m_do3dt2m_conus_gmi(gmi_lat, gmi_lon, do3dt2m, r, lat_castnet, 
    lon_castnet, r_castnet, do3dt2m_castnet, case):
    """function plots a map of the ozone-temperature sensitivity (dO3/dT2m)
    over the continential U.S. 

    Parameters
    ----------  
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
    do3dt2m : numpy.ndarray    
        The O3-T2m sensitivity at each GMI grid cell [lat, lon]   
    r : numpy.ndarray
        Pearson correlation coefficients between MERRA-2 2 meter temperatures
        and GMI ozone at each grid cell [lat, lon]
    lat_castnet : list
        CASTNet latitude coordinates, units of degrees north, [no. CASTNet
        stations,]
    lon_castnet : list
        CASTNet longitude coordinates, units of degrees north, [no. CASTNet
        stations,]    
    r_castnet : list
        Pearson correlation coefficients between MERRA-2 2 meter temperatures
        and CASTNet ozone at each CASTNet site, [no. CASTNet stations,]
    do3dt2m_castnet : list
        The O3-T2m sensitivity (slope of linear regression) at each CASTNet
        site, [no. CASTNet stations,]        
    case : str
        Model simulation name for output file

    Returns
    ----------     
    None
    """
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.colorbar import ColorbarBase
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    llcrnrlon = -126.
    llcrnrlat = 24.0
    urcrnrlon = -66.3
    urcrnrlat = 50.
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
        llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat, 
        resolution = 'h', area_thresh = 1000)
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)    
    fig = plt.figure(figsize=(5, 7.5))
    # O3-T2m correlation map
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
                llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, 
                urcrnrlat = urcrnrlat, resolution = 'h', area_thresh = 1000)
    ax = plt.subplot2grid((2, 1), (0, 0))
    ax.set_title('(a)', fontsize = 16, x = 0.1, y = 1.03)
    vmin = 0; vmax = 0.75
    cmap = plt.get_cmap('PuBu', 10)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m.contourf(x, y, r, clevs, cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    # superimpose O3-T2m sensitivity from observations
    x_castnet, y_castnet = m(lon_castnet, lat_castnet)
    m.scatter(x_castnet, y_castnet, c = r_castnet, s = 18, cmap = cmap, 
              vmin = vmin, vmax = vmax, zorder = 30, linewidth = 1., 
              edgecolor = 'orange')    
    fill_oceans(ax, m)    
    outline_region(ax, m, pollutants_constants.NORTHEAST_STATES)    
    caxt = fig.add_axes([0.21, 0.56, 0.6, 0.03]) 
    cb = ColorbarBase(caxt, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = r'$r$(T, O$_{\mathregular{3}}$)', size = 16, y = 0.25)    
    cb.ax.tick_params(labelsize = 12)
    cb.set_ticks(np.linspace(vmin, vmax, 6))    
    # O3-T2m sensitivity map         
    ax2 = plt.subplot2grid((2, 1), (1, 0))
    ax2.set_title('(b)', fontsize = 16, x = 0.1, y = 1.03)
    # set up colormap
    vmin = 0.0; vmax = 2.0
    cmap = plt.get_cmap('PuBu', 10)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m.contourf(x, y, do3dt2m, clevs, cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    # superimpose O3-T2m sensitivity from observations
    x_castnet, y_castnet = m(lon_castnet, lat_castnet)
    m.scatter(x_castnet, y_castnet, c = do3dt2m_castnet, s = 18, cmap = cmap, 
              vmin = vmin, vmax = vmax, zorder = 30, linewidth = 1., 
              edgecolor = 'orange')
    # mask oceans and outline NEUS
    fill_oceans(ax2, m)    
    outline_region(ax2, m, pollutants_constants.NORTHEAST_STATES)
    # add colorbar
    caxb = fig.add_axes([0.21, 0.13, 0.6, 0.03])
    cb = ColorbarBase(caxb, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_ticks(np.linspace(vmin, vmax, 6))
    cb.set_label(label = '$\mathregular{d}$O$_{\mathregular{3}}$/' + 
                 '$\mathregular{d}$T' +
                 '[ppbv K$^{\mathregular{-1}}$]', 
                 size = 16, y = 0.2)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.18, hspace=0.6)                 
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'map_ro3dt2m_do3dt2m_conus_gmi_%s.eps' %case, dpi = 350)    
    return
# # # # # # # # # # # # #
def map_sensitivityratio_conus(gmi_lat, gmi_lon, dat_r, dat_sens, mr2_sens):
    """function plot a map of the 2-meter temperature and O3 from the Transport
    simulation (a) correlation, (b) slope, and (c) the ratio of the Transport 
    O3-temperature sensitivity to the + Chemistry ozone-temperature 
    sensitivity (i.e. dO3/dT2m). A value of 1 would mean that the sensitivity 
    didn't change between simulations whereas a value of 0.5 would mean that 
    the sensitivity is halved between simulations.
    
    Parameters
    ----------  
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
    dat_r : numpy.ndarray
        The Pearson product-moment correlation coefficient calculated between 
        MERRA-2 2-meter temperatures and ozone at each GMI grid cell from the 
        Transport simulation, [lat, lon]        
    dat_sens : numpy.ndarray    
        The O3-T2m sensitivity at each GMI grid cell in the Transport 
        simulation [lat, lon]   
    mr2_sens : numpy.ndarray    
        The O3-T2m sensitivity at each GMI grid cell in the + Chemistry 
        simulation [lat, lon]           

    Returns
    ----------     
    None    
    """
    import numpy as np
    from matplotlib.colors import Normalize
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # continental U.S. focus region map 
    m = Basemap(projection = 'merc', llcrnrlon = -126., llcrnrlat = 24.0, 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'h', 
                area_thresh = 1000)
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)
    fig = plt.figure(figsize = (6, 8))
    # left subplot, correlation coefficients in Transport simulatiton
    ax1 = plt.subplot2grid((3, 1), (0, 0))
    ax1.set_title('(a)', fontsize = 16, x = 0.1, y = 1.03)
    vmin = 0.0; vmax = 0.75
    cmap = plt.get_cmap('PuBu', 10)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m1 = m.contourf(x, y, dat_r, clevs, cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax1, m)    
    outline_region(ax1, m, pollutants_constants.NORTHEAST_STATES)    
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size = '5%', pad = 0.25)
    cb = fig.colorbar(m1, cax = cax, norm = norm, orientation = 'vertical')
    cb.set_label(label = r'$r$(O$_{\mathregular{3}}$' + 
                 ', T$_{\mathregular{2\:m}}$)', labelpad = 8., size = 16)
    cb.ax.tick_params(labelsize = 12)
    # middle subplot, O3-climate penalty in Transport simulation
    ax2 = plt.subplot2grid((3, 1), (1, 0))
    ax2.set_title('(b)', fontsize = 16, x = 0.1, y = 1.03)
    vmin = 0.0; vmax = 2.0
    cmap = plt.get_cmap('PuBu', 10)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m2 = m.contourf(x, y, dat_sens, clevs, cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax2, m)    
    outline_region(ax2, m, pollutants_constants.NORTHEAST_STATES)    
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size = '5%', pad = 0.25)
    cb = fig.colorbar(m2, cax = cax, norm = norm, orientation = 'vertical')
    cb.set_label(label = '$\mathregular{\partial}$O$_{\mathregular{3}}$ ' + 
                 '$\mathregular{\partial}$T$_{\mathregular{2 m}}^' +
                 '{\mathregular{-1}}$ [ppbv K$^{\mathregular{-1}}$]', 
                 labelpad = 8., size = 16)
    cb.ax.tick_params(labelsize = 12)
    # right subplot, ratio of slope in Transport/+ Chemistry simulations
    ax3 = plt.subplot2grid((3, 1), (2, 0))
    ax3.set_title('(c)', fontsize = 16, x = 0.1, y = 1.03)
    vmin = 0.0; vmax = 1.0
    cmap = plt.get_cmap('gist_heat', 10)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m3 = m.contourf(x, y, dat_sens/mr2_sens, clevs, cmap = cmap, 
                    extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax3, m)    
    outline_region(ax3, m, pollutants_constants.NORTHEAST_STATES)    
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size = '5%', pad = 0.25)
    cb = fig.colorbar(m3, cax = cax,  norm = norm, orientation = 'vertical')
    cb.set_label(label = 'Ratio $\mathregular{\partial}$O$_{\mathregular{3}}$ ' + 
                 '$\mathregular{\partial}$T$_{\mathregular{2 m}}^' +
                 '{\mathregular{-1}}$', labelpad = 8., size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(hspace = 0.35)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_sensitivityratio_conus.eps', dpi = 300)  
    return 
# # # # # # # # # # # # #    
def map_std_90ptile(t2m_overpass, o3, dat_o3, gmi_lat, gmi_lon):
    """function plots standard deviations of MERRA-2 2-meter temperature 
    fields at overpass2 time, O3 from Transport and + Chemistry simulations, 
    the difference of standard deviation between O3 from Transport and 
    + Chemistry simulations, and the ratio of these simulations' standard 
    deviations. 9 Oct 2018: in addition to the aforementioned plots, function 
    also plots the same but using the 90th percentile of O3 as the quantity 
    plotted.
    
    Parameters
    ----------  
    t2m_atoverpass_interp : numpy.ndarry
        The MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]
    o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from MR2/+ Chemistry 
        simulation, units of volume mixing ratio, [time, lat, lon]        
    dat_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from Diurnal-AvgT/
        Transport simulation, units of volume mixing ratio, [time, lat, lon]      
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  

    Returns
    ----------     
    None      
    """
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.colorbar import ColorbarBase
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from matplotlib.cm import ScalarMappable   
    from MidpointNormalize import MidpointNormalize    
    m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 25., 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'h', 
                area_thresh = 1000)
    # map of 2-meter temperature variability
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)
    vmin = 2.0; vmax = 6.0
    cmap = plt.get_cmap('terrain_r', 8)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 9, endpoint = True)
    m.contourf(x, y, np.std(t2m_overpass, axis = 0), clevs, cmap = cmap, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    fill_oceans(ax, m)
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = '$\mathregular{\sigma}_{\mathregular{T_{\mathregular{2 m}}}}$ [K]', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/map_stdt2m_overpass.eps', 
                dpi = 300)
    # map of O3 variability in + Chemistry simulation (standard deviation)
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    vmin = 4.0; vmax = 12.0
    cmap = plt.get_cmap('terrain_r', 8)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 9, endpoint = True)
    m.contourf(x, y, np.std(o3 * 1e9, axis = 0), clevs, cmap = cmap, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    fill_oceans(ax, m)
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = '$\mathregular{\sigma}_{\mathregular{O_{' +
                 '\mathregular{3, + Chemistry}}}}$ [ppbv]', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/map_stdo3_conus_gmi_MR2.eps', 
                dpi = 300)
    # (90th percentile)
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    vmin = 40.; vmax = 65.
    cmap = plt.get_cmap('terrain_r', 10)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m.contourf(x, y, np.percentile(o3 * 1e9, 90, axis = 0), clevs, cmap = cmap, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    fill_oceans(ax, m)
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = 'P$_{\mathregular{90}}$' + 
                 ' O$_{\mathregular{3, + Chemistry}}$ [ppbv]', size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/map_90ptileo3_conus_gmi_MR2.eps', 
                dpi = 300)
    # map of O3 variability in Transport simulation (standard deviation)
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    vmin = 4.0; vmax = 12.0
    cmap = plt.get_cmap('terrain_r', 8)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 9, endpoint = True)
    m.contourf(x, y, np.std(dat_o3 * 1e9, axis = 0), clevs, cmap = cmap, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    fill_oceans(ax, m)
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = '$\mathregular{\sigma}_{\mathregular{O_{' +
                 '\mathregular{3, Transport}}}}$ [ppbv]', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'map_stdo3_conus_gmi_Diurnal-AvgT.eps', dpi = 300)
    # (90th percentile)
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    vmin = 40.; vmax = 65.
    cmap = plt.get_cmap('terrain_r', 10)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m.contourf(x, y, np.percentile(dat_o3 * 1e9, 90, axis = 0), clevs, cmap = cmap, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    fill_oceans(ax, m)
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = 'P$_{\mathregular{90}}$' + 
                 ' O$_{\mathregular{3, Transport}}$ [ppbv]', size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/map_90ptileo3_conus_gmi_Diurnal-AvgT.eps', 
                dpi = 300)
    # map of the O3 variability difference between + Chemistry and Transport 
    # simulations
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    vmin = -2.; vmax = 2.
    cmap = plt.get_cmap('bwr', 8)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 9, endpoint = True)
    m.contourf(x, y, np.std(o3 * 1e9, axis = 0) - 
               np.std(dat_o3 * 1e9, axis = 0), clevs, cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    fill_oceans(ax, m)
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = '$\mathregular{\Delta\sigma}_{\mathregular{O_' +
                 '{\mathregular{3, + Chemistry - Transport}}}}$ [ppbv]', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'map_diffstdo3_conus_gmi_MR2_Diurnal-AvgT.eps', dpi = 300)
    # difference between 90th percentile 
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    vmin = -2.; vmax = 2.
    cmap = plt.get_cmap('bwr', 8)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 9, endpoint = True)
    m.contourf(x, y, np.percentile(o3 * 1e9, 90, axis = 0) - 
               np.percentile(dat_o3 * 1e9, 90, axis = 0), clevs, 
               cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    fill_oceans(ax, m)
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = 'P$_{\mathregular{90}}$ O$_{\mathregular{3, + Chemistry}}$' + 
                 '$-$ P$_{\mathregular{90}}$ O$_{\mathregular{3, Transport}}$', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'map_diff90ptileo3_conus_gmi_MR2_Diurnal-AvgT.eps', dpi = 300)
    # map of the ratio of O3 variability between Transport and + Chemistry 
    # simulations    
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    vmin = 0.78; vmax = 1.06
    cmap = plt.get_cmap('Spectral_r', 14)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 15, endpoint = True)
    m.contourf(x, y, np.std(dat_o3 * 1e9, axis = 0)/
               np.std(o3 * 1e9, axis = 0), clevs, cmap = cmap, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)
    fill_oceans(ax, m)    
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = '$\mathregular{\sigma}_{\mathregular{O_' +
                 '{\mathregular{3, Transport}}}}:\mathregular{\sigma}_' +
                 '{\mathregular{O_{\mathregular{3, +Chemistry}}}}$', 
                 size = 16)
    cb.set_ticks(np.linspace(vmin, vmax, 8))
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'map_ratiostdo3_conus_gmi_MR2_Diurnal-AvgT.eps', dpi = 300)
    # ratio of 90th percentiles
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    vmin = 0.96; vmax = 1.04
    cmap = plt.get_cmap('terrain_r', 8)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 9, endpoint = True)
    m.contourf(x, y, np.percentile(o3 * 1e9, 90, axis = 0)/
               np.percentile(dat_o3 * 1e9, 90, axis = 0), clevs, cmap = cmap, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    fill_oceans(ax, m)    
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = 'P$_{\mathregular{90}}$ O$_{\mathregular{3, + ' + 
                 'Chemistry}}$ : P$_{\mathregular{90}}$ O$_{\mathregular{3, Transport}}$', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'map_ratio90ptileo3_conus_gmi_MR2_Diurnal-AvgT.eps', dpi = 300)
    return
# # # # # # # # # # # # #    
def timeseriesscatter_castnetdo3dt2m_cases(do3dt2m_castnet, t_castnet, 
    o3_castnet): 
    """function identifies the highest, median, and lowest O3-T2m sensitivities
    from CASTNet/MERRA-2 data in the continental U.S. and plots timeseries 
    and scatterplots of O3, T2m for these three cases.
    
    Parameters
    ----------  
    do3dt2m_castnet : list
        The O3-T2m sensitivity (slope of linear regression fit between O3 
        versus T2m) at each CASTNet site, units of ppbv K^-1, [no. sites,]
    t_castnet : list
        Daily 1300 hours local time MERRA-2 2-meter temperatures at each 
        CASTNet site, units of K, [no. sites,]
    o3_castnet : list
        Daily 1300 hours local time O3 at each CASTNet site, units of ppbv, 
        [no. sites,]

    Returns
    ----------     
    None     
    """
    import numpy as np
    import matplotlib.pyplot as plt
    # initialize figure, axes
    fig = plt.figure()
    axt = plt.subplot2grid((3, 4), (0, 0), colspan = 3)
    axm = plt.subplot2grid((3, 4), (1, 0), colspan = 3)
    axb = plt.subplot2grid((3, 4), (2, 0), colspan = 3)
    axtr = plt.subplot2grid((3, 4), (0, 3))
    axmr = plt.subplot2grid((3, 4), (1, 3))
    axbr = plt.subplot2grid((3, 4), (2, 3))
    # CASTNet sites with highest, mid-range, and lowest O3-T2m sensitivity
    high_idx = np.where(np.array(do3dt2m_castnet) == 
                        np.array(do3dt2m_castnet).max())[0][0]
    medium_idx = do3dt2m_castnet.index(np.percentile(do3dt2m_castnet, 50, 
                                                     interpolation = 'nearest'))
    low_idx = np.where(np.array(do3dt2m_castnet) == 
                       np.array(do3dt2m_castnet).min())[0][0]
    # mask for missing O3 observations
    high_mask = ~np.isnan(t_castnet[high_idx]) & ~np.isnan(o3_castnet[high_idx])
    medium_mask = ~np.isnan(t_castnet[medium_idx]) & ~np.isnan(o3_castnet[medium_idx])
    low_mask = ~np.isnan(t_castnet[low_idx]) & ~np.isnan(o3_castnet[low_idx])
    ylim1 = [20, 85]
    ylim2 = [288, 315]
    # max. sensitivity
    axt.plot(o3_castnet[high_idx])
    axt.set_ylim(ylim1)
    axt2 = axt.twinx()
    axt2.plot(t_castnet[high_idx], '-r')
    axt2.set_ylim(ylim2)
    axt2.set_title('Max. O$_{3}$-T$_{2m}$ sensitivity (= %.2f ppbv K$^{-1}$)' 
                   %do3dt2m_castnet[high_idx])
    axtr.scatter(t_castnet[high_idx][high_mask], 
                 o3_castnet[high_idx][high_mask], c = 'k', s = 5)
    axtr.set_xlim(ylim2)
    axtr.set_ylim(ylim1)
    # median sensitivity 
    axm.plot(o3_castnet[medium_idx])
    axm.set_ylim(ylim1)
    axm2 = axm.twinx()
    axm2.plot(t_castnet[medium_idx], '-r')
    axm2.set_ylim(ylim2)
    axm2.set_title('Median O$_{3}$-T$_{2m}$ sensitivity (= %.2f ppbv K$^{-1}$)' 
                   %do3dt2m_castnet[medium_idx])
    axmr.scatter(t_castnet[medium_idx][medium_mask], 
                 o3_castnet[medium_idx][medium_mask], c = 'k', s = 5)
    axmr.set_xlim(ylim2)
    axmr.set_ylim(ylim1)
    #  min. sensitivity
    axb.plot(o3_castnet[low_idx], label = 'O$_{3}$')
    axb.set_ylim(ylim1)
    axb.legend(loc = 8, frameon = False)
    axb2 = axb.twinx()
    axb2.plot(t_castnet[low_idx], '-r', label = 'T$_{2m}$')
    axb2.set_ylim(ylim2)
    axb2.legend(loc = 4, frameon = False)
    axb2.set_title('Min. O$_{3}$-T$_{2m}$ sensitivity (= %.2f ppbv K$^{-1}$)' 
                   %do3dt2m_castnet[low_idx])
    axbr.scatter(t_castnet[low_idx][low_mask], 
                 o3_castnet[low_idx][low_mask], c = 'k', s = 5)
    axbr.set_xlim(ylim2)
    axbr.set_ylim(ylim1)
    plt.tight_layout()
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'timeseriesscatter_castnetdo3dt2m_cases.eps', dpi = 300)    
    return 
# # # # # # # # # # # # #    
def map_rdato3mr2o3_conus(o3, dat_o3, gmi_lat, gmi_lon):
    """function calculates Pearson correlation coefficient between Transport
    and + Chemistry simulations and plots map. 
    
    Parameters
    ----------  
    o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from MR2/+ Chemistry 
        simulation, units of volume mixing ratio, [time, lat, lon]        
    dat_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from Diurnal-AvgT/
        Transport simulation, units of volume mixing ratio, [time, lat, lon]      
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  

    Returns
    ----------     
    None     
    """
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.colorbar import ColorbarBase
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    o3 = o3 * 1e9
    dat_o3 = dat_o3 * 1e9
    # calculate Pearson correlation coefficient between Transport and 
    # + Chemistry simulations
    mr2_dat_r = np.empty(shape = o3.shape[1:])
    mr2_dat_r[:] = np.nan
    for i, ilat in enumerate(gmi_lat):
        for j, jlon in enumerate(gmi_lon):
            mr2_dat_r[i, j] = np.corrcoef(o3[:, i, j], 
                                          dat_o3[:, i, j], 1)[0, 1]
    m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 24.8, 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'h', 
                area_thresh = 1000)
    # map of the correlation between Transport and + Chemistry simulations
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)
    vmin = 0.85; vmax = 1.0
    cmap = plt.get_cmap('terrain_r', 15)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 16, endpoint = True)
    m.contourf(x, y, mr2_dat_r, clevs, cmap = cmap, extend = 'min')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'min')
    cb.set_label(label = r'$\rho$(O$_{\mathregular{3, Transport}}$' + 
                 ', O$_{\mathregular{3, + Chemistry}}$)', size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/map_rdato3mr2o3_conus.eps', 
                dpi = 300)
    return
# # # # # # # # # # # # #    
def map_rmr2o3gridboxheight_conus(o3, gridboxheight, gmi_lat, gmi_lon): 
    """function plots the Pearson correlation coefficient between ozone 
    from the MR2/+ Chemistry simulation and the CTM grid box height for the 
    first model level.

    Parameters
    ----------  
    o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from MR2/+ Chemistry 
        simulation, units of volume mixing ratio, [time, lat, lon]        
    gridboxheight : numpy.ndarray
        GMI CTM grid box height for model levels, units of m, [time, level, 
        lat, lon]      
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  

    Returns
    ----------     
    None     
    """
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.colorbar import ColorbarBase
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap    
    r = np.empty(shape = gridboxheight.shape[2:])
    r[:] = np.nan
    for i, ilat in enumerate(gmi_lat):
        for j, jlon in enumerate(gmi_lon):
            r[i, j] = np.corrcoef(gridboxheight[:, 0, i, j], 
                                  o3[:, i, j], 1)[0, 1]
    m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 24.8, 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'h', 
                area_thresh = 1000)
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)
    vmin = -0.4; vmax = 0.8
    cmap = plt.get_cmap('PuBu', 6)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 7, endpoint = True)
    m.contourf(x, y, r, clevs, cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)  
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = r'$\rho$(O$_{\mathregular{3, + Chemistry}}$' + 
                 ', grid box height$_{\mathregular{surface}}$)', size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/map_rmr2o3gridboxheight_conus.eps', 
                dpi = 300)
    return 
# # # # # # # # # # # # #    
def map_meanmr2o3meancastneto3_conus(o3, o3_castnet, gmi_lat, gmi_lon,
    lat_castnet, lon_castnet, region_castnet, sites_castnet):
    """function plots mean O3 fields from + Chemistry/MR2 simulation and 
    mean CASTNet observations for the measuring period 2008-2010. To illustrate
    the model resolution, parallels and meridians at the CTM's resolution 
    are also plotted. A second figure showing the bias at each CASTNet site 
    is plotted and the mean CASTNet-GMI bias and the bias in a specified region
    are output. Here bias = GMI O3 - CASTNet O3.
    
    Parameters
    ----------  
    o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from MR2/+ Chemistry 
        simulation, units of volume mixing ratio, [time, lat, lon]        
    o3_castnet : list
        Daily 1300 hours local time O3 at each CASTNet site, units of ppbv, 
        [no. sites,]        
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]    
    lat_castnet : list
        Latitudes of CASTNet sites, units of degrees north, [no. sites,]
    lon_castnet : list
        Longitudes of CASTNet sites, units of degrees west, [no. sites,]
    region_castnet : list
        Site IDs for CASTNet sites in region, [no. sites in region,]
    sites_castnet : list
        CASTNet site names, [no. sites,]        
        
    Returns
    ----------     
    None 
    """
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.colorbar import ColorbarBase
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    from geo_idx import geo_idx    
    # # # # plot mean JJA O3 as simulated by + Emissions/HindcastMR2 simulation
    o3_mean = np.mean(o3, axis = 0)
    # from CASTNet
    o3_castnet_mean = np.nanmean(np.array(o3_castnet), axis = 1)
    m = Basemap(projection = 'merc', llcrnrlon = -126., llcrnrlat = 24., 
            urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'h', 
            area_thresh = 1000)
    fig = plt.figure(figsize=(5, 7.5))
    ax = plt.subplot2grid((2, 1), (0, 0))
    ax.set_title('(a)', fontsize = 16, x = 0.03, y = 1.03, ha = 'left')    
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)
    vmin = 26; vmax = 56
    cmap = plt.get_cmap('gist_earth', 15)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 16, endpoint = True)
    m.contourf(x, y, o3_mean * 1e9, clevs, cmap = cmap, extend = 'both')
    fill_oceans(ax, m)
    # draw parallels and meridians to indicate the resolution of the CTM
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(gmi_lat[0], gmi_lat[-1], np.diff(gmi_lat).mean())
    # labels = [left, right, top, bottom]
    m.drawparallels(parallels, color = 'lightgrey', linewidth = 0.3,
                    dashes = [1, 0], labels = [False, False, False, False])
    meridians = np.arange(gmi_lon[0], gmi_lon[-1], np.diff(gmi_lon).mean())
    m.drawmeridians(meridians, color = 'lightgrey', linewidth = 0.3, 
                    dashes = [1, 0], labels = [False, False, False, False])
    # add mean CASTNet data
    x_castnet, y_castnet = m(lon_castnet, lat_castnet)
    m.scatter(x_castnet, y_castnet, c = o3_castnet_mean, s = 20, cmap = cmap, 
              vmin = vmin, vmax = vmax, zorder = 30, linewidth = 1., 
              edgecolor = 'k')    
    # add countries, states, etc here so they sit atop the gridlines
    m.drawstates(color = 'k', linewidth = 0.5, zorder = 10)
    m.drawcountries(color = 'k', linewidth = 1.0, zorder = 10)
    m.drawcoastlines(color = 'k', linewidth = 1.0, zorder = 10)  
    outline_region(ax, m, pollutants_constants.NORTHEAST_STATES)    
    caxt = fig.add_axes([0.21, 0.56, 0.6, 0.03])    
    cb = ColorbarBase(caxt, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_ticks(np.linspace(vmin, vmax, 6))
    cb.set_label(label = '$\mathregular{O}_3}$ [ppbv]', size=16)
    cb.ax.tick_params(labelsize = 12)
    # # # # plot bias at each CASTNet site
    i = 0
    bias_all = []
    for ilat, ilon in zip(lat_castnet, lon_castnet): 
        # loop through CASTNet sites and find nearest GMI grid cell
        lat_gmi_near = geo_idx(ilat, gmi_lat)
        lon_gmi_near = geo_idx(ilon, gmi_lon)
        # O3 at GMI grid cell 
        gmi_near_castnet = o3[:, lat_gmi_near, lon_gmi_near] * 1e9
        castnet_atsite = o3_castnet[i]
        bias_atsite = np.nanmean(gmi_near_castnet)  - np.nanmean(castnet_atsite)
        bias_all.append(bias_atsite)
        i = i + 1
    print('Max bias of %.3f ppbv at %s' %(np.max(bias_all),
          sites_castnet[np.where(np.array(bias_all) == np.max(bias_all))[0][0]]))
    print('Min bias of %.3f ppbv at %s' %(np.min(bias_all),
          sites_castnet[np.where(np.array(bias_all) == np.min(bias_all))[0][0]]))    
    ax2 = plt.subplot2grid((2, 1), (1, 0))
    ax2.set_title('(b)', fontsize = 16, x = 0.03, y = 1.03, ha = 'left')
    vmin = -12; vmax = 12
    cmap = plt.get_cmap('RdBu_r', 12)
    clevs = np.linspace(vmin, vmax, 13, endpoint = True)
    sc = m.scatter(x_castnet, y_castnet, c = np.array(bias_all), s = 20, 
              cmap = cmap, vmin = vmin, vmax = vmax,
              zorder = 30, linewidth = 1., edgecolor = 'k')    
    fill_oceans(ax2, m)    
    # add countries, states, etc here so they sit atop the gridlines
    m.drawstates(color = 'k', linewidth = 0.5, zorder = 10)
    m.drawcountries(color = 'k', linewidth = 1.0, zorder = 10)
    m.drawcoastlines(color = 'k', linewidth = 1.0, zorder = 10)  
    outline_region(ax2, m, pollutants_constants.NORTHEAST_STATES)    
    cax = fig.add_axes([0.21, 0.13, 0.6, 0.03])
    cb = fig.colorbar(sc, cax=cax, orientation='horizontal', extend = 'both')
    cb.set_ticks(np.linspace(vmin, vmax, 7, endpoint = True))
    cb.set_label(label = 'Bias [ppbv]', size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.18, hspace=0.6)
    print('mean nationwide bias (GMI - CASTNet) = %.4f ppbv' 
          %np.mean(bias_all))
    # find O3 bias in region 
    where_region = np.in1d(sites_castnet, region_castnet)
    where_region = np.where(where_region == True)[0]
    bias_region = np.array(bias_all)[where_region]
    print('mean bias in region (GMI - CASTNet) = %.4f ppbv'
          %np.mean(bias_region))
    plt.savefig('/Users/ghkerr/phd/GMI/figs/'
                'map_meanemisso3meancastneto3_withbias_conus.eps', dpi = 300)
    return
# # # # # # # # # # # # #  
def scatter_castnett2mo3_gmit2mo3_slopes(castnet_o3_region, castnet_t2m_region, 
    mr2_o3_region, mr2_t2m_region, region):
    """function plots regionally-averaged O3 vs. 2-meter temperatures from
    CASTNet (left subplot) and GMI (right subplot) and finds (1) least squares
    regression, (2) Theil-Sen regression, and (3) quadratic regression and 
    plots these best fit lines for both observations and modeled O3. 
    
    Parameters
    ----------  
    castnet_o3_region : numpy.ndarray
        Daily regionally-averaged O3 from CASTNet, [time,]         
    castnet_t2m_region : numpy.ndarray
        Daily regionally-averaged 2-meter temperatures co-located with CASTNet 
        stations, [time,]        
    mr2_o3_region : numpy.ndarray
        Daily regionally-averaged O3, [time,]      
    mr2_t2m_region : numpy.ndarray
        Daily regionally-averaged T2m, [time,]
    region : string
        Name or abbreviation of region for output file
     
    Returns
    ----------     
    None
    """
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize = (9, 4))
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))    
    # regression using Theil-Sen with 95% confidence intervals; note: syntax is
    # stats.theilslopes(y, x, CI)
    res_castnet = stats.theilslopes(castnet_o3_region, castnet_t2m_region, 0.95)
    res_mr2 = stats.theilslopes(mr2_o3_region, mr2_t2m_region, 0.95)
    # least squares regression
    lsq_res_castnet = stats.linregress(castnet_t2m_region, castnet_o3_region)
    lsq_res_mr2 = stats.linregress(mr2_t2m_region, mr2_o3_region)
    # Jing & Steiner [2017] found that a second-order regression was most 
    # suitable for the O3-T2m relationship; form is O3 = b0 + (b1 * T2m) + 
    # (b2 * t2m^2)
    so_res_mr2 = np.polyfit(mr2_t2m_region, mr2_o3_region, 2) # returns 
    # coefficients with the highest power first
    # using polyfit with the CASTNet data leads to a poorly conditioned fit
    # i.e., np.polyfit(castnet_t2m_neus, castnet_o3_neus, 2)
    # instead, use least squares solution 
    X = np.vstack((np.ones(castnet_t2m_region.shape), 
                  castnet_t2m_region, 
                  castnet_t2m_region**2))
    X = X.T
    b = np.matmul(np.matmul(np.linalg.inv(np.matmul(X.T, X)), X.T), 
                  castnet_o3_region)
    # total least squares/orthogonal distance regression
    # adapted from https://www.tutorialspoint.com/scipy/scipy_odr.htm and
    # https://stackoverflow.com/questions/24804397/numpy-polyfit-versus-scipy-odr
    # define function for scipy.odr
    import scipy.odr        
    def fit_func(p, t):
        return p[0] * t + p[1]
    # fit the data using scipy.odr for observations and model 
    Model = scipy.odr.Model(fit_func)
    Data = scipy.odr.RealData(castnet_t2m_region, castnet_o3_region)
    Odr = scipy.odr.ODR(Data, Model, [0., 0.], maxit = 10000)
    output = Odr.run()
    beta = output.beta
    #betastd = output.sd_beta    
    Data_gmi = scipy.odr.RealData(mr2_t2m_region, mr2_o3_region)
    Odr_gmi = scipy.odr.ODR(Data_gmi, Model, [0., 0.], maxit = 10000)
    output_gmi = Odr_gmi.run()
    beta_gmi = output_gmi.beta
    #betastd = output_gmi.sd_beta        
    # first subplot, CASTNet O3
    ax1.plot(castnet_t2m_region, castnet_o3_region, 'o', markersize = 3, 
             color = 'darkgrey', zorder = 1)
#    # total least squares
#    ax1.plot(castnet_t2m_region, fit_func(beta, castnet_t2m_region), 
#             '#a6cee3', lw = 2, zorder = 2)    
    castnet_t2m_region = np.sort(castnet_t2m_region)
    # least squares
    ax1.plot(castnet_t2m_region, lsq_res_castnet[1] + lsq_res_castnet[0] * 
             castnet_t2m_region, '-', color = '#1f78b4', lw = 2., zorder = 5)
    # Theil-Sen
    ax1.plot(castnet_t2m_region, res_castnet[1] + res_castnet[0] * 
             castnet_t2m_region, '--',  color = '#fb9a99', lw = 2., zorder = 6)
    # second-order
    ax1.plot(castnet_t2m_region, (castnet_t2m_region**2 * b[2]) + 
             (castnet_t2m_region * b[1]) + b[0], '-', zorder = 3,
             color = '#33a02c', lw = 2.)
    xlim = ax1.get_xlim()
    ylim = ax1.get_ylim()
    for t in ax1.get_xticklabels():
        t.set_fontsize(12)
    for t in ax1.get_yticklabels():
        t.set_fontsize(12)    
    ax1.set_title('(a) CASTNet', fontsize = 16, x = 0.03, y = 1.03, ha = 'left')
    ax1.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    ax1.get_yaxis().set_label_coords(-0.15, 0.53)
    ax1.set_xlabel('T$_{}$ [K]', fontsize = 16, x = 1.1)
    ax1.tick_params(right = True, left = True)   
    print('CASTNet Total Least Squares slope = %.3f ppbv K^-1' %beta[0])
    print('CASTNet Least Squares slope = %.3f ppbv K^-1' %lsq_res_castnet.slope)
    print('CASTNet Theil-Sen slope = %.3f ppbv K^-1' %res_castnet[0])
    print('CASTNet second-order slope, b1, b2, b3 = %.5f, %.5f, %.5f'
          %(b[0], b[1], b[2]))
    # second subplot, GMI O3
    ax2.plot(mr2_t2m_region, mr2_o3_region, 'o', markersize = 3, 
             color = 'darkgrey', zorder = 1)
#    # total least squares
#    ax2.plot(mr2_t2m_region, fit_func(beta_gmi, mr2_t2m_region), 
#             '#a6cee3', lw = 2, label = 'TLS', 
#             zorder = 2)    
    mr2_t2m_region = np.sort(mr2_t2m_region)
    # least squares
    ax2.plot(mr2_t2m_region, lsq_res_mr2[1] + lsq_res_mr2[0] * mr2_t2m_region, 
             '-', color = '#1f78b4', lw = 2., label = 'OLS', zorder = 5)
    # Theil-Sen
    ax2.plot(mr2_t2m_region, res_mr2[1] + res_mr2[0] * mr2_t2m_region, '--',  
             color = '#fb9a99', lw = 2., label = 'Theil-Sen', zorder = 6)
    # second-order
    ax2.plot(mr2_t2m_region, (mr2_t2m_region**2 * so_res_mr2[0]) + 
             (mr2_t2m_region * so_res_mr2[1]) + so_res_mr2[2], '-', 
             color = '#33a02c', lw = 2., label = 'Second order', zorder = 3)
    # set x/y limits of second subplot to be consistent with first 
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)
    for t in ax2.get_xticklabels():
        t.set_fontsize(12)
    ax2.set_yticklabels([''])
    ax2.yaxis.tick_right()
    ax2.set_title('(b) CTM', fontsize = 16, x = 0.03, y = 1.03, ha = 'left')
    plt.subplots_adjust(bottom = 0.28)
    ax2.legend(loc = 3, ncol = 4, frameon = False, fontsize = 16, 
               bbox_to_anchor = (-0.9, -0.48))
    # uncomment and replace if TLS is included
#               bbox_to_anchor = (-1.2, -0.48))    
    print('GMI Total Least Squares slope = %.3f ppbv K^-1' %beta_gmi[0])    
    print('GMI Least Squares slope = %.3f ppbv K^-1' %lsq_res_mr2.slope)
    print('GMI Theil-Sen slope = %.3f ppbv K^-1' %res_mr2[0])
    print('GMI second-order slope, b1, b2, b3 = %.5f, %.5f, %.5f'
          %(so_res_mr2[2], so_res_mr2[1], so_res_mr2[0]))
    plt.savefig('/Users/ghkerr/phd/GMI/figs/'
                'scatter_castnett2mo3_gmit2mo3_slopeswithoutTLS_%s.eps' 
                %region, dpi = 300)    
    return 
# # # # # # # # # # # # #    
def timeseries_castneto3allgmio3(transport, chemistry, emissions, obs, region,
    year, years):
    """function plots time series of regionally-averaged O3 from GMI 
    Transport, + Chemistry, and + Emissions simulations and regionally-averaged
    CASTNet observations. Correlation coefficients for the year of interest and 
    for the multi-year period are calculated. 
    
    Parameters
    ----------       
    transport : numpy.ndarray
        GMI CTM daily regionally-averaged O3 concentrations for model case 
        HindcastMR2-DiurnalAvgT, units of ppbv, [time,]
    chemistry : numpy.ndarray
        GMI CTM daily regionally-averaged O3 concentrations for model case 
        HindcastMR2, units of ppbv, [time,]
    emissions : numpy.ndarray
        GMI CTM daily regionally-averaged O3 concentrations for model case 
        GHKerr-DailyEmiss, units of ppbv, [time,]
    obs : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [time,]    
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
    year : int
        Year of interest
    years : list
        Years in measuring period
        
    Returns
    ----------      
    None         
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # for single year
    yearpos = np.where(np.array(years) == year)[0][0]    
    fig = plt.figure(figsize = (10, 4))
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.plot(obs[yearpos*92:(yearpos+1)*92], '-', color = 'darkgrey', 
             zorder = 1, lw = 2, label = 'CASTNet')
    ax.plot(emissions[yearpos*92:(yearpos+1)*92], '-', 
            color = pollutants_constants.COLOR_EMISSIONS, 
            label = '$+$AEmissions')
    ax.plot(chemistry[yearpos*92:(yearpos+1)*92],
            '-', color = pollutants_constants.COLOR_CHEMISTRY, 
            label = '$+$Chemistry')
    ax.plot(transport[yearpos*92:(yearpos+1)*92], '-', 
            color = pollutants_constants.COLOR_TRANSPORT, label = 'Transport')
    ax.set_xlim([0, 91])
    ax.set_ylim([30, 65])
    # axis labels
    ax.set_xticks([0, 14, 30, 44, 61, 75])
    ax.set_xticklabels(['1 June %d' %year, '', '1 July', '', '1 Aug', ''], 
                       ha = 'center', fontsize = 12)
    ax.set_ylabel('Ozone [ppbv]', fontsize = 16)
    for tl in ax.get_yticklabels():
        tl.set_fontsize(12)
    ax.tick_params(right = True, left = True, top = True)
    # reorder and place legend
    handles, labels = ax.get_legend_handles_labels()
    handles = [handles[0], handles[3], handles[2], handles[1]]
    labels = [labels[0], labels[3], labels[2], labels[1]]
    ax.legend(handles,labels, bbox_to_anchor = (0.49, -0.18), loc = 'center', 
              ncol = 4, frameon = False, fontsize = 16)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/'
                'timeseries_castneto3allgmio3_%d_%s.eps' %(year, region), 
                dpi = 300)      
    # print correlation coefficients for year of interest and for 2008-2010
    print('CASTNet/Transport correlation = %.3f for %d' %(
            np.corrcoef(obs[yearpos*92:(yearpos+1)*92], 
            transport[yearpos*92:(yearpos+1)*92])[0, 1], year))
    print('CASTNet/Transport correlation = %.3f for all years' %(
            np.corrcoef(obs, transport)[0, 1]))
    print('CASTNet/+Chemistry correlation = %.3f for %d' %(
            np.corrcoef(obs[yearpos*92:(yearpos+1)*92], 
            chemistry[yearpos*92:(yearpos+1)*92])[0, 1], year))
    print('CASTNet/+Chemistry correlation = %.3f for all years' %(
            np.corrcoef(obs, chemistry)[0, 1]))
    print('CASTNet/+Emissions correlation = %.3f for %d' %(
            np.corrcoef(obs[yearpos*92:(yearpos+1)*92], 
            emissions[yearpos*92:(yearpos+1)*92])[0, 1], year)) 
    print('CASTNet/+Emissions correlation = %.3f for all years' %(
            np.corrcoef(obs, emissions)[0, 1]))
    print('Transport/+Emissions correlation = %.3f for %d' %(
            np.corrcoef(transport[yearpos*92:(yearpos+1)*92], 
            emissions[yearpos*92:(yearpos+1)*92])[0, 1], year)) 
    print('Transport/+Emissions correlation = %.3f for all years ' %(
            np.corrcoef(transport, emissions)[0, 1]))
    print('Transport/+Chemistry correlation = %.3f for %d' %(
            np.corrcoef(transport[yearpos*92:(yearpos+1)*92], 
            chemistry[yearpos*92:(yearpos+1)*92])[0, 1], year)) 
    print('Transport/+Chemistry correlation = %.3f for all years' %(
            np.corrcoef(transport, chemistry)[0, 1])) 
    print('+Chemistry/+Emissions correlation = %.3f for %d' %(
            np.corrcoef(chemistry[yearpos*92:(yearpos+1)*92], 
            emissions[yearpos*92:(yearpos+1)*92])[0, 1], year)) 
    print('+Chemistry/+Emissions correlation = %.3f for all years' %(
            np.corrcoef(chemistry, emissions)[0, 1])) 
    return 
# # # # # # # # # # # # #    
def scatterhist_castneto3allgmio3(transport, chemistry, emissions, obs, t2m, 
    castnet_t2m, region):
    """function plots a scatterplot of regionally-averaged O3 versus MERRA-2
    2-meter temperatures for the transport, chemistry, and emissions 
    simulations and observations. KDEs of distributions are plotted along the 
    spines of the scatterplot. 

    Parameters
    ----------       
    transport : numpy.ndarray
        GMI CTM daily regionally-averaged O3 concentrations for model case 
        HindcastMR2-DiurnalAvgT, units of ppbv, [time,]
    chemistry : numpy.ndarray
        GMI CTM daily regionally-averaged O3 concentrations for model case 
        HindcastMR2, units of ppbv, [time,]
    emissions : numpy.ndarray
        GMI CTM daily regionally-averaged O3 concentrations for model case 
        GHKerr-DailyEmiss, units of ppbv, [time,]
    obs : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [time,] 
    t2m : numpy.ndarray
        MERRA-2 2-meter temperatures interpolated to resolution of CTM and 
        averaged over all grid cells in region, [time,]
    castnet_t2m : 
        MERRA-2 2-meter temperatures interpolated to the resolution of CTM and 
        sampled at grid cells containing CASTNet sites in a particular region 
        and averaged over grid cells in this region, [time, ]
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
    year : int
        Year of interest
    years : list
        Years in measuring period
        
    Returns
    ----------      
    None         
    """
    import numpy as np
    import scipy.stats as stats
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    nullfmt = NullFormatter()
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    # start with a rectangular figure
    fig = plt.figure(1, figsize = (8, 8))
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    # remove labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    # scatter plot, change absolute O3 concentrations versus 
    # absolute temperature 
    axScatter.scatter(castnet_t2m, obs, label = 'CASTNet', s = 20, 
                      color = 'darkgrey', zorder = 2)    
    axScatter.scatter(t2m, transport, label = 'Transport', s = 20, 
                      color = pollutants_constants.COLOR_TRANSPORT, zorder = 9)
    axScatter.scatter(t2m, chemistry, label = '$\mathregular{+}$Chemistry', s = 20, 
                      color = pollutants_constants.COLOR_CHEMISTRY, zorder = 6)
    axScatter.scatter(t2m, emissions, label = '$\mathregular{+}$AEmissions', s = 20, 
                      color = pollutants_constants.COLOR_EMISSIONS, zorder = 3)
    # add in lines of best fit 
    obs_m = np.poly1d(np.polyfit(castnet_t2m, obs, 1))[1]
    transport_m = np.poly1d(np.polyfit(t2m, transport, 1))[1]
    chemistry_m = np.poly1d(np.polyfit(t2m, chemistry, 1))[1]
    emissions_m = np.poly1d(np.polyfit(t2m, emissions, 1))[1]
    axScatter.plot(np.unique(t2m), np.poly1d(np.polyfit(t2m, obs, 
        1))(np.unique(t2m)), zorder = 30, color = 'darkgrey', lw = 2.5)
    axScatter.plot(np.unique(t2m), np.poly1d(np.polyfit(t2m, transport, 
        1))(np.unique(t2m)), zorder = 30,
        color = pollutants_constants.COLOR_TRANSPORT, lw = 2.5)
    axScatter.plot(np.unique(t2m), np.poly1d(np.polyfit(t2m, chemistry, 
        1))(np.unique(t2m)), zorder = 20,
        color = pollutants_constants.COLOR_CHEMISTRY, lw = 2.5)
    axScatter.plot(np.unique(t2m), np.poly1d(np.polyfit(t2m, emissions, 
        1))(np.unique(t2m)), zorder = 10, 
        color = pollutants_constants.COLOR_EMISSIONS, lw = 2.5)    
    axScatter.set_xlabel('T [K]', fontsize = 16)
    axScatter.set_ylabel('Ozone [ppbv]', fontsize = 16)
    for tl in axScatter.get_xticklabels():
        tl.set_fontsize(12)
    for tl in axScatter.get_yticklabels():
        tl.set_fontsize(12)    
    axScatter.get_xaxis().set_label_coords(0.5, -0.07)        
    axScatter.get_yaxis().set_label_coords(-0.09, 0.5) 
    axScatter.text(0.04, 0.92, '(a)',
                 transform = axScatter.transAxes, fontsize = 20)  
    # manually set xticks and remove final tick to avoid clutter 
    axScatter.set_xlim([290, 308])
    axScatter.set_xticks([290, 292, 294, 296, 298, 300, 302, 304, 306, 308])
    plt.setp(axScatter.get_xticklabels()[0], visible=False)
    plt.setp(axScatter.get_xticklabels()[-1], visible=False)
    # histograms
    nbins = 18
    n, x, _  = axHistx.hist(t2m, bins = nbins, histtype = u'step', 
                            density = True, color = 'k', lw = 0.)
    density = stats.gaussian_kde(t2m)
    axHistx.plot(x, density(x), lw = 2., color = 'k')
    axHistx.set_ylabel('Density', fontsize = 16)
    #axHistx.text(0.25, 0.7, 'T$_{\mathregular{2\:m}}$', 
    #             transform = axHistx.transAxes, fontsize = 20)   
    axHistx.text(0.04, 0.92, '(b)',
                 transform = axHistx.transAxes, fontsize = 20)    
    axHistx.get_yaxis().set_label_coords(-0.11, 0.5)        
    # observations
    n, x, _  = axHisty.hist(obs, bins = nbins, histtype = u'step', 
                            orientation = 'horizontal', density = True, lw = 0.)
    density = stats.gaussian_kde(obs)
    axHisty.plot(density(x), x, zorder = 2, lw = 2.,color = 'darkgrey')
    # transport
    n, x, _  = axHisty.hist(transport, bins = nbins, histtype = u'step', 
                            orientation = 'horizontal', density = True, lw = 0.)
    density = stats.gaussian_kde(transport)
    axHisty.plot(density(x), x, zorder = 9, lw = 2.,
                 color = pollutants_constants.COLOR_TRANSPORT)
    # chemistry 
    n, x, _  = axHisty.hist(chemistry, bins = nbins, histtype = u'step', 
                            orientation = 'horizontal', density = True, lw = 0.)
    density = stats.gaussian_kde(chemistry)
    axHisty.plot(density(x), x, zorder = 6, lw = 2.,
                 color = pollutants_constants.COLOR_CHEMISTRY)
    # emissions
    n, x, _  = axHisty.hist(emissions, bins = nbins, histtype = u'step', 
                            orientation = 'horizontal', density = True, lw = 0.)
    density = stats.gaussian_kde(emissions)
    axHisty.plot(density(x), x, zorder = 3, lw = 2.,
                 color = pollutants_constants.COLOR_EMISSIONS)
    axHisty.set_xlabel('Density', fontsize = 16)
    #axHisty.text(0.25, 0.8, 'Ozone', transform = axHisty.transAxes, 
    #             fontsize = 20)
    axHisty.text(0.25, 0.92, '(c)', transform = axHisty.transAxes, 
                 fontsize = 20)
    axHisty.get_xaxis().set_label_coords(0.5, -0.07)        
    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_xlim([0, 0.1])
    axHistx.set_ylim([0, 0.15])
    axHistx.spines['top'].set_visible(False)
    axHistx.spines['right'].set_visible(False)
    axHisty.set_ylim(axScatter.get_ylim())
    axHisty.spines['top'].set_visible(False)
    axHisty.spines['right'].set_visible(False)
    axScatter.legend(ncol = 4, frameon = False, fontsize = 16, 
                     bbox_to_anchor = (1.45, -0.12))
    for ax in [axScatter, axHistx, axHisty]:
        for tl in ax.get_xticklabels():
            tl.set_fontsize(12)
        for tl in ax.get_yticklabels():
            tl.set_fontsize(12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatterhist_castneto3allgmio3_%s.eps' %(region), dpi = 300, 
                bbox_inches = 'tight')     
    return 
# # # # # # # # # # # # #   
def scatter_rt2memisso3_demisso3t2m_latlon_conus(emiss_sens, emiss_r, gmi_lat, 
    gmi_lon):
    """function plots the O3-climate penalty versus the correlation coefficient
    calculated between 2-meter temperatures and O3 at every land-based grid 
    cell in CTM. The colorscale corresponds to the latitude (left subplot) and
    the longitude (right subplot) of the individual grid cell.
    
    Parameters
    ----------  
    emiss_sens : numpy.ndarray    
        The O3-climate penalty at each GMI grid cell from the + Emissions 
        simulation, units of ppbv K^-1, [lat, lon]         
    emiss_r : numpy.ndarray
        The Pearson product-moment correlation coefficient calculated between 
        MERRA-2 2-meter temperatures and ozone at each GMI grid cell from the 
        + Emissions simulation, [lat, lon]
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  

    Returns
    ----------     
    None
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # find grid cells in CONUS to create land-sea mask 
    m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 24.0, 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'c', 
                area_thresh = 1000)
    states = pollutants_constants.NORTHEAST_STATES + \
        pollutants_constants.SOUTH_STATES + pollutants_constants.MIDWEST_STATES + \
        pollutants_constants.IMW_STATES + pollutants_constants.WEST_STATES
    lsmask = find_grid_in_region(m, states, gmi_lat, gmi_lon)
    # O3-climate penalty, O3-T2m correlation on land from simulations
    emiss_sens_land = emiss_sens * lsmask
    emiss_r_land = emiss_r * lsmask
    # latitude, longitude coordinates on land
    lon_gmi_mesh, lat_gmi_mesh = np.meshgrid(gmi_lon, gmi_lat)
    lon_gmi_mesh = lon_gmi_mesh * lsmask
    lat_gmi_mesh = lat_gmi_mesh * lsmask
    # stack into 1D array
    emiss_sens_land = np.hstack(emiss_sens_land)
    emiss_r_land = np.hstack(emiss_r_land)
    lat_gmi_mesh = np.hstack(lat_gmi_mesh)
    lon_gmi_mesh = np.hstack(lon_gmi_mesh)    
    fig = plt.figure(figsize = (8, 4))
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    # O3-climate penalty versus O3-T2m correlation colored by latitude
    mp1 = ax1.scatter(emiss_r_land, emiss_sens_land, c = np.hstack(lat_gmi_mesh), 
                      s = 7)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size = '5%', pad = 0.05)
    fig.colorbar(mp1, cax = cax, orientation = 'vertical', 
                 label = 'Latitude [$^{\circ}$N]')
    ax1.set_xlabel('r(T$_{\mathregular{2 m}}$, O$_{\mathregular{3}}$)', 
                   size = 16)
    ax1.set_ylabel('$\mathregular{\partial}$O$_{\mathregular{3}}$ ' + 
                   '$\mathregular{\partial}$T$_{\mathregular{2 m}}^' +
                   '{\mathregular{-1}}$ [ppbv K$^{\mathregular{-1}}$]', 
                   size = 16)
    # O3-climate penalty versus O3-T2m correlation colored by longitude
    mp2 = ax2.scatter(emiss_r_land, emiss_sens_land, 
                      c = np.hstack(lon_gmi_mesh), s = 7)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size = '5%', pad = 0.05)
    fig.colorbar(mp2, cax = cax, orientation = 'vertical',
                 label = 'Longitude [$^{\circ}$E]')             
    ax2.set_xlabel('r(T$_{\mathregular{2 m}}$, O$_{\mathregular{3}}$)', 
                   size = 16)
    for ax in [ax1, ax2]:
        for tl in ax.get_xticklabels():
            tl.set_fontsize(12)
        for tl in ax.get_yticklabels():
            tl.set_fontsize(12)
    plt.subplots_adjust(wspace = 0.3, bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_rt2memisso3_demisso3t2m_latlon_conus.eps', dpi = 300)       
    return
# # # # # # # # # # # # #    
def scatter_dmr2o3dt2m_ddato3dt2m_rt2memisso3_conus(mr2_sens, dat_sens, r, 
    gmi_lat, gmi_lon): 
    """function plots the O3-climate penalty from the Transport simulation 
    versus the penalty from the + Chemistry simulation for every GMI grid
    cell over land in the CONUS. Scatterpoints' colors are the correlation 
    coefficient calculated between 2-meter temperatures and O3 from the 
    + Emissions simulation. 

    Parameters
    ----------  
    mr2_sens : numpy.ndarray    
        The O3-climate penalty at each GMI grid cell from the + Chemistry
        simulation, units of ppbv K^-1, [lat, lon]         
    dat_sens : numpy.ndarray    
        The O3-climate penalty at each GMI grid cell from the Transport
        simulation, units of ppbv K^-1, [lat, lon]        
    r : numpy.ndarray
        The Pearson product-moment correlation coefficient calculated between 
        MERRA-2 2-meter temperatures and ozone at each GMI grid cell from the 
        + Emissions simulation, [lat, lon]
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  

    Returns
    ----------     
    None    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # find grid cells in CONUS to create land-sea mask 
    m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 24.0, 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'c', 
                area_thresh = 1000)
    states = pollutants_constants.NORTHEAST_STATES + \
        pollutants_constants.SOUTH_STATES + pollutants_constants.MIDWEST_STATES + \
        pollutants_constants.IMW_STATES + pollutants_constants.WEST_STATES
    lsmask = find_grid_in_region(m, states, gmi_lat, gmi_lon)
    r_land = r * lsmask
    mr2_sens_land = mr2_sens * lsmask
    dat_sens_land = dat_sens * lsmask
    # latitude, longitude coordinates on land
    lon_gmi_mesh, lat_gmi_mesh = np.meshgrid(gmi_lon, gmi_lat)
    lon_gmi_mesh = lon_gmi_mesh * lsmask
    lat_gmi_mesh = lat_gmi_mesh * lsmask
    # stack into 1D array
    mr2_sens_land = np.hstack(mr2_sens_land)
    dat_sens_land = np.hstack(dat_sens_land)
    r_land = np.hstack(r_land)
    lat_gmi_mesh = np.hstack(lat_gmi_mesh)
    lon_gmi_mesh = np.hstack(lon_gmi_mesh)
    # initialize figure, axis
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    cmap = plt.get_cmap('PuBu', 10)
    mp = ax.scatter(mr2_sens_land, dat_sens_land, c = r_land, cmap = cmap, 
                    s = 6)
    # plot 1:1, 1:2, 1:3 lines 
    ax.plot(np.linspace(-0.5, 3.5, 50), np.linspace(-0.5, 3.5, 50), '--k')
    ax.plot(np.linspace(-0.5, 3.5, 50), 
            0.5 * np.linspace(-0.5, 3.5, 50), '--k')
    ax.plot(np.linspace(-0.5, 3.5, 50), 
            (1/3.) * np.linspace(-0.5, 3.5, 50), '--k')
    ax.plot(np.linspace(-0.5, 3.5, 50), 
            (1/4.) * np.linspace(-0.5, 3.5, 50), '--k')
    ax.set_xlim([0, 2.5])
    ax.set_ylim([-0.5, 1.5])
    ax.set_xlabel('$\mathregular{\partial}$O$_{\mathregular{3,+\:Chemistry}}$ ' + 
                  '$\mathregular{\partial}$T$_{\mathregular{2\:m}}^' +
                  '{\mathregular{-1}}$ [ppbv K$^{\mathregular{-1}}$]', 
                   size = 16)    
    ax.set_ylabel('$\mathregular{\partial}$O$_{\mathregular{3, Transport}}$ ' + 
                  '$\mathregular{\partial}$T$_{\mathregular{2\:m}}^' +
                  '{\mathregular{-1}}$ [ppbv K$^{\mathregular{-1}}$]', 
                   size = 16)    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '5%', pad = 0.05)
    fig.colorbar(mp, cax = cax, orientation = 'vertical',
                 label = 'r(T$_{\mathregular{2 m}}$, O$_{\mathregular{3}}$)')
    plt.subplots_adjust(left = 0.15, bottom = 0.15)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_dmr2o3dt2m_ddato3dt2m_rt2memisso3_conus.eps', 
                dpi = 300)      
    return 
# # # # # # # # # # # # #    
def scatter_lat_rt2memisso3_conus(emiss_r, mr2_sens, dat_sens, emiss_sens, 
    lat_gmi, lon_gmi, lat_castnet, r_castnet, do3dt2m_castnet):
    """function plots two figures: (1) the correlation coefficients calculated 
    between O3 from the + Emissions simulation and 2-meter temperatures versus
    grid cell latitude for every latitude over land. Colors in this figure are 
    the O3-climate penalty from the + Emissions simulation. (2) the same as 
    (1) but the colors correspond to the percent change in the O3-climate 
    penalty between the Transport and + Chemistry simulations.
    
    Parameters
    ----------  
    emiss_r : numpy.ndarray
        The Pearson product-moment correlation coefficient calculated between 
        MERRA-2 2-meter temperatures and ozone at each GMI grid cell from the 
        + Emissions simulation, [lat, lon]
    mr2_sens : numpy.ndarray    
        The O3-climate penalty at each GMI grid cell from the + Chemistry
        simulation, units of ppbv K^-1, [lat, lon]         
    dat_sens : numpy.ndarray    
        The O3-climate penalty at each GMI grid cell from the Transport
        simulation, units of ppbv K^-1, [lat, lon]  
    emiss_sens : numpy.ndarray    
        The O3-climate penalty at each GMI grid cell from the + Emissions
        simulation, units of ppbv K^-1, [lat, lon]          
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
    lat_castnet : list
        CASTNet latitude coordinates, units of degrees north, [no. CASTNet
        stations,]
    r_castnet : list
        Pearson correlation coefficients between MERRA-2 2 meter temperatures
        and CASTNet ozone at each CASTNet site, [no. CASTNet stations,]
    do3dt2m_castnet : list
        The O3-T2m sensitivity (slope of linear regression) at each CASTNet
        site, [no. CASTNet stations,]         
        
    Returns
    ----------     
    None        
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # find grid cells in CONUS to create land-sea mask 
    m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 24.0, 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'c', 
                area_thresh = 1000)
    states = pollutants_constants.NORTHEAST_STATES + \
        pollutants_constants.SOUTH_STATES + pollutants_constants.MIDWEST_STATES + \
        pollutants_constants.IMW_STATES + pollutants_constants.WEST_STATES
    lsmask = find_grid_in_region(m, states, gmi_lat, gmi_lon)  
    r_land = emiss_r * lsmask
    mr2_sens_land = mr2_sens * lsmask
    dat_sens_land = dat_sens * lsmask
    emiss_sens_land = emiss_sens * lsmask
    # latitude, longitude coordinates on land
    lon_gmi_mesh, lat_gmi_mesh = np.meshgrid(gmi_lon, gmi_lat)
    lon_gmi_mesh = lon_gmi_mesh * lsmask
    lat_gmi_mesh = lat_gmi_mesh * lsmask
    # stack into 1D array
    mr2_sens_land = np.hstack(mr2_sens_land)
    dat_sens_land = np.hstack(dat_sens_land)
    emiss_sens_land = np.hstack(emiss_sens_land)    
    r_land = np.hstack(r_land)
    lat_gmi_mesh = np.hstack(lat_gmi_mesh)
    lon_gmi_mesh = np.hstack(lon_gmi_mesh)
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    mp = ax.scatter(lat_gmi_mesh, r_land, c = (emiss_sens_land), s = 12, 
                    vmin = 0., vmax = 1.5, zorder = 3)
    ax.set_xlabel('Latitude [$^{\circ}$N]')
    ax.set_ylabel('r(T$_{\mathregular{2 m}}$, O$_{\mathregular{3}}$)')
    fig.colorbar(mp, extend = 'both', 
                 label = '$\mathregular{\partial}$O$_{\mathregular{3}}$ ' + 
                 '$\mathregular{\partial}$T$_{\mathregular{2 m}}^' +
                 '{\mathregular{-1}}$ [ppbv K$^{\mathregular{-1}}$]')
    ax.scatter(lat_castnet, r_castnet,  c = (do3dt2m_castnet), s = 40, 
               vmin = 0., vmax = 1.5, edgecolor = 'k', zorder = 5)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_lat_rt2memisso3_demisso3t2m_conus.eps', 
                dpi = 300)
    plt.show()
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    pc = (mr2_sens_land - dat_sens_land)/(mr2_sens_land) * 100.
    mp = ax.scatter(lat_gmi_mesh, r_land, c = (pc), s = 12, vmin = 0, 
                    vmax = 100)
    fig.colorbar(mp, extend = 'both', label = 'Percent Change [%]')
    ax.set_xlabel('Latitude [$^{\circ}$N]')
    ax.set_ylabel('r(T$_{\mathregular{2 m}}$, O$_{\mathregular{3}}$)')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_lat_rt2memisso3_percentchange_conus.eps', 
                dpi = 300)    
    plt.show()
    return
# # # # # # # # # # # # #
def map_allgmio3_hotcold(dat_o3, mr2_o3, emiss_o3, t2m_overpass, gmi_lat, 
    gmi_lon, neus):
    """for each GMI simulation function finds the difference in O3 between hot
    (i.e., 2-meter temperatures > 90th percentile) and median days and cold
    (i.e., 2-meter temperatures < 10th percentile) and median days. These 
    differences for every simulation are cast in terms of their contribution 
    to the total difference in the + Emissions simulation. For example for the 
    "Transport" subplot, this is the hot-median O3 average from the Transport 
    simulation divided by the hot-median O3 average from the + Emissions 
    simulations. For the "+ Emissions" subplot it is the difference in
    the hot-median O3 average between + Emissions and + Chemistry simulations
    divided by the difference in the + Emissions simulation. Function also 
    prints out Northeast-averaged values. 
    
    Parameters
    ----------       
    dat_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, for model case 
        HindcastMR2-DiurnalAvgT, units of volume mixing ratio, [time, lat, lon] 
    mr2_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, for model case 
        HindcastMR2, units of volume mixing ratio, [time, lat, lon] 
    emiss_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, for model case 
        GHKerr-DailyEmiss, units of volume mixing ratio, [time, lat, lon] 
    t2m_overpass : numpy.ndarray
        MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
    neus : numpy.ndarray
        CTM grid where grid cells within NEUS have a value of 1 and grid 
        cells not in region have a value of NaN, [lat, lon]         
        
    Returns
    ----------      
    None             
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from matplotlib.colors import Normalize
    from matplotlib.colors import LinearSegmentedColormap    
    from matplotlib.colorbar import ColorbarBase   
    # from https://stackoverflow.com/questions/3279560/reverse-colormap-in-matplotlib
    def reverse_colourmap(cmap, name = 'my_cmap_r'):
        """
        In: 
        cmap, name 
        Out:
        my_cmap_r
    
        Explanation:
        t[0] goes from 0 to 1
        row i:   x  y0  y1 -> t[0] t[1] t[2]
                       /
                      /
        row i+1: x  y0  y1 -> t[n] t[1] t[2]
    
        so the inverse should do the same:
        row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                       /
                      /
        row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
        """        
        reverse = []
        k = []   
        for key in cmap._segmentdata:    
            k.append(key)
            channel = cmap._segmentdata[key]
            data = []
            for t in channel:                    
                data.append((1-t[0],t[2],t[1]))            
            reverse.append(sorted(data))    
        LinearL = dict(zip(k,reverse))
        my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL, N=10) 
        return my_cmap_r    
    # find mean O3 in each simulation on days where 2-meter temperatures exceed
    # the 90th percentile
    dat_o3_p90 = o3_givent2mpercentile(dat_o3, t2m_overpass, 90, 'greater')
    mr2_o3_p90 = o3_givent2mpercentile(mr2_o3, t2m_overpass, 90, 'greater')
    # are less than the 10th percentile 
    dat_o3_p10 = o3_givent2mpercentile(dat_o3, t2m_overpass, 10, 'less')
    mr2_o3_p10 = o3_givent2mpercentile(mr2_o3, t2m_overpass, 10, 'less')
    # find median O3 values in each simulation and convert from volume mixing
    # ratio to ppbv
    emiss_o3_med = np.median(emiss_o3, axis = 0) * 1e9
    # O3 on hot days versus median days from simulation
    delta_dat_hot = dat_o3_p90 - emiss_o3_med
    delta_mr2_hot = mr2_o3_p90 - emiss_o3_med
    # O3 on cold days versus median days 
    delta_dat_cold = dat_o3_p10 - emiss_o3_med
    delta_mr2_cold = mr2_o3_p10 - emiss_o3_med
    # map and colorscheme
    # for continental U.S. 
    m = Basemap(projection = 'merc', llcrnrlon = -126., llcrnrlat = 24., 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'h', 
                area_thresh = 1000)
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)    
    vmin = 0.; vmax = 12.
    cmap = plt.get_cmap('PuBu', 12)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 13, endpoint = True)
    # # # # plot for hot days 
    # contributions from Transport simulation
    fig = plt.figure(figsize = (4, 7))    
    ax1 = plt.subplot2grid((3, 1), (0, 0))
    ax1.set_title('(a) Transport', ha = 'left', fontsize = 16, x = 0.03, 
                  y = 1.03)
    m.contourf(x, y, delta_dat_hot, clevs, cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax1, m)    
    # contributions from + Chemistry simulation
    ax2 = plt.subplot2grid((3, 1), (1, 0))
    ax2.set_title('(b) $\mathregular{+}$Chemistry', ha = 'left', fontsize = 16, x = 0.03, 
                  y = 1.03)
    m.contourf(x, y, delta_mr2_hot, clevs, cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax2, m)    
    # colorbar for first two plots 
    plt.subplots_adjust(right = 0.7)
    cbar_ax = fig.add_axes([0.75, 0.43, 0.05, 0.4])
    cb = ColorbarBase(cbar_ax, cmap = cmap, norm = norm, extend = 'both',
                      orientation = 'vertical')
    cb.set_ticks(np.linspace(vmin, vmax, 7))
    cb.set_label(label = 'O$_{\mathregular{3}}$(T$_{\mathregular{P}' +
                 '_{\mathregular{90}}})$ $-$ O$_{\mathregular{3}}$' +
                 '(T$_{\mathregular{P}_{\mathregular{50}}})$ [ppbv]', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)
    # ratio of simulations
    ax3 = plt.subplot2grid((3, 1), (2, 0))
    ax3.set_title('(c) Transport/$\mathregular{+}$Chemistry', ha = 'left', fontsize = 16, 
                  x = 0.03, y = 1.03)
    vmin = 0.; vmax = 100.
    # only use half of RdBu colormap since there aren't many places exceeding
    # 100% by evaluating an existing colormap from 0.5 (midpoint) to 1 
    # (upper end)
    cmap = plt.get_cmap('RdBu', 10)
    colors = cmap(np.linspace(0.5, 1, cmap.N // 2))
    # create a new colormap from those colors
    cmap = LinearSegmentedColormap.from_list('Upper Half', colors, 10)
    #cmap = reverse_colourmap(cmap)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m.contourf(x, y, (delta_dat_hot/delta_mr2_hot) * 100., clevs, cmap = cmap,
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)    
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax3, m)    
    m.contour(x, y, emiss_r, levels = [0.4], colors = ('orange',),
              linestyles = ('--',), linewidths = (1,), zorder = 1)     
    cbar_ax = fig.add_axes([0.75, 0.16, 0.05, 0.16])
    cb = ColorbarBase(cbar_ax, cmap = cmap, norm = norm, extend = 'both',
                      orientation = 'vertical')
    cb.set_ticks(np.linspace(vmin, vmax, 6))
    cb.set_label(label = '[%]', size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_allgmio3_p90mp50_hot.eps', dpi = 300)        
    plt.show()
    # # # #
    vmin = -12.; vmax = 0.
    cmap = plt.get_cmap('PuBu_r', 12)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 13, endpoint = True)
    # contributions from Transport simulation
    fig = plt.figure(figsize = (4, 7))    
    ax1 = plt.subplot2grid((3, 1), (0, 0))
    ax1.set_title('(a) Transport', ha = 'left', fontsize = 16, x = 0.03, 
                  y = 1.03)
    m.contourf(x, y, delta_dat_cold, clevs, cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax1, m)    
    # contributions from + Chemistry simulation
    ax2 = plt.subplot2grid((3, 1), (1, 0))
    ax2.set_title('(b) $\mathregular{+}$Chemistry', ha = 'left', fontsize = 16, x = 0.03, 
                  y = 1.03)
    m.contourf(x, y, delta_mr2_cold, clevs, cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax2, m)    
    # colorbar for first two plots 
    plt.subplots_adjust(right = 0.7)
    cbar_ax = fig.add_axes([0.75, 0.43, 0.05, 0.4])
    cb = ColorbarBase(cbar_ax, cmap = cmap, norm = norm, extend = 'both',
                      orientation = 'vertical')
    cb.set_ticks(np.linspace(vmin, vmax, 7))
    cb.set_label(label = 'O$_{\mathregular{3}}$(T$_{\mathregular{P}' +
                 '_{\mathregular{10}}})$ $-$ O$_{\mathregular{3}}$' +
                 '(T$_{\mathregular{P}_{\mathregular{50}}})$ [ppbv]', size = 16)
    cb.ax.tick_params(labelsize = 12)
    # ratio from simulations
    ax3 = plt.subplot2grid((3, 1), (2, 0))
    ax3.set_title('(c) Transport/$\mathregular{+}$Chemistry', ha = 'left', fontsize = 16, 
                  x = 0.03, y = 1.03)
    vmin = 0.; vmax = 100.
    cmap = plt.get_cmap('RdBu', 10)
    colors = cmap(np.linspace(0.5, 1, cmap.N // 2))
    # create a new colormap from those colors
    cmap = LinearSegmentedColormap.from_list('Upper Half', colors, 10)
    #cmap = reverse_colourmap(cmap)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m.contourf(x, y, (delta_dat_cold/delta_mr2_cold) * 100., clevs, 
               cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)    
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax3, m)    
    m.contour(x, y, emiss_r, levels = [0.4], colors = ('orange',),
              linestyles = ('--',), linewidths = (1,), zorder = 1)     
    cbar_ax = fig.add_axes([0.75, 0.16, 0.05, 0.16])
    cb = ColorbarBase(cbar_ax, cmap = cmap, norm = norm, extend = 'both',
                      orientation = 'vertical')
    cb.set_ticks(np.linspace(vmin, vmax, 6))
    cb.set_label(label = '[%]', size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_allgmio3_p10mp50_cold.eps', dpi = 300)     
    plt.show()
    return 
# # # # # # # # # # # # #    
def map_allgmio3_do3dt(dat_sens, mr2_sens, emiss_sens, emiss_r, gmi_lat, 
    gmi_lon, fstr):
    """plot maps of the O3-climate penalty from Transport, + Chemistry, and
    + Emissions simulatios. Outline regions where the O3-T correlation (from
    + Emissions simulation) falls below r = 0.3.

    Parameters
    ----------       
    dat_sens : numpy.ndarray    
        The O3-climate penalty at each GMI grid cell from the + Emissions 
        simulation, units of ppbv K^-1, [lat, lon]         
    mr2_sens : numpy.ndarray    
        The O3-climate penalty at each GMI grid cell from the + Emissions 
        simulation, units of ppbv K^-1, [lat, lon]         
    emiss_sens : numpy.ndarray    
        The O3-climate penalty at each GMI grid cell from the + Emissions 
        simulation, units of ppbv K^-1, [lat, lon]         
    emiss_r : numpy.ndarray
        The Pearson product-moment correlation coefficient calculated between 
        MERRA-2 2-meter temperatures and ozone at each GMI grid cell from the 
        + Emissions simulation, [lat, lon]
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
    fstr : str
        Output filename suffix
        
    Returns
    ----------      
    None                
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from matplotlib.colors import Normalize
    from matplotlib.colorbar import ColorbarBase    
    from matplotlib.colors import LinearSegmentedColormap    
    # from https://stackoverflow.com/questions/3279560/reverse-colormap-in-matplotlib
    def reverse_colourmap(cmap, name = 'my_cmap_r'):
        """
        In: 
        cmap, name 
        Out:
        my_cmap_r
    
        Explanation:
        t[0] goes from 0 to 1
        row i:   x  y0  y1 -> t[0] t[1] t[2]
                       /
                      /
        row i+1: x  y0  y1 -> t[n] t[1] t[2]
    
        so the inverse should do the same:
        row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                       /
                      /
        row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
        """        
        reverse = []
        k = []   
        for key in cmap._segmentdata:    
            k.append(key)
            channel = cmap._segmentdata[key]
            data = []
            for t in channel:                    
                data.append((1-t[0],t[2],t[1]))            
            reverse.append(sorted(data))    
        LinearL = dict(zip(k,reverse))
        my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL, N=10) 
        return my_cmap_r 
    m = Basemap(projection = 'merc', llcrnrlon = -126., llcrnrlat = 24., 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'h', 
                area_thresh = 1000)
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)    
    vmin = 0.0; vmax = 2.0
    cmap = plt.get_cmap('PuBu', 8)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 9, endpoint = True)
    # contributions from Transport simulation
    fig = plt.figure(figsize = (4, 7))    
    ax1 = plt.subplot2grid((3, 1), (0, 0))
    ax1.set_title('(a) Transport', ha = 'left', fontsize = 16, x = 0.03, 
                  y = 1.03)
    m.contourf(x, y, dat_sens, clevs, cmap = cmap, norm = norm, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax1, m)    
    # + Chemistry simulation
    ax2 = plt.subplot2grid((3, 1), (1, 0))
    ax2.set_title('(b) $\mathregular{+}$Chemistry', ha = 'left', fontsize = 16, x = 0.03, 
                  y = 1.03)
    m.contourf(x, y, mr2_sens, clevs, cmap = cmap, norm = norm, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax2, m)    
    plt.subplots_adjust(right = 0.7)
    cbar_ax = fig.add_axes([0.75, 0.43, 0.05, 0.4])
    cb = ColorbarBase(cbar_ax, cmap = cmap, norm = norm, extend = 'both',
                      orientation = 'vertical')
    cb.set_ticks(np.linspace(vmin, vmax, 5))
    cb.set_label(label = '$\mathregular{d}$O$_{\mathregular{3}}$/' + 
                 '$\mathregular{d}$T ' +
                 '[ppbv K$^{\mathregular{-1}}$]', size = 16)    
    cb.ax.tick_params(labelsize = 12)    
    # ratio from simulations
    ax3 = plt.subplot2grid((3, 1), (2, 0))
    ax3.set_title('(c) Transport/$\mathregular{+}$Chemistry', ha = 'left', fontsize = 16, 
                  x = 0.03, y = 1.03)
    vmin = 0.; vmax = 100.
    cmap = plt.get_cmap('RdBu', 10)
    colors = cmap(np.linspace(0.5, 1, cmap.N // 2))
    # create a new colormap from those colors
    cmap = LinearSegmentedColormap.from_list('Upper Half', colors, 10)
    #cmap = reverse_colourmap(cmap)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m.contourf(x, y, (dat_sens/mr2_sens) * 100., clevs, cmap = cmap, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)    
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax3, m)    
    m.contour(x, y, emiss_r, levels = [0.4], colors = ('orange',),
              linestyles = ('--',), linewidths = (1,), zorder = 1)     
    cbar_ax = fig.add_axes([0.75, 0.16, 0.05, 0.16])
    cb = ColorbarBase(cbar_ax, cmap = cmap, norm = norm, extend = 'both',
                      orientation = 'vertical')
    cb.set_ticks(np.linspace(vmin, vmax, 5))
    cb.set_label(label = '[%]', size = 16)
    cb.ax.tick_params(labelsize = 12)    
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_allgmio3_do3dt%s.eps' %fstr, dpi = 300)         
    return        
# # # # # # # # # # # # #  
def timeseries_t2m_castneto3_cemsnox(castnet_o3, castnet_t2m, dat_o3_neus, 
    mr2_o3_neus, emiss_o3_neus, year, years, region):
    """for a given summer (JJA) function plots regionally-averaged 2-meter
    temperatures from MERRA-2, regionally-averaged O3 from CASTNet, and 
    regionally-summed NOx emissions from CEMS and determines the correlation 
    coefficients between these variables. A scatterplot of O3 versus 
    temperature is also plotted. 

    Parameters
    ----------       
    castnet_o3 : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period, days in months in 'sampling_months']   
    castnet_t2m : numpy.ndarray
        MERRA-2 2-meter temperatures co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of K, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']  
    year : int
        Year of interest
    years : list
        Years in measuring period
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
        
    Returns
    ----------      
    None             
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/emissions/')
    import AQSCEMSobs      
    sampling_months = [6, 7, 8]
    # load CEMS NOx emissions in NEUS
    states_ab = ['CT', 'DC', 'DE', 'MA', 'MD', 'ME', 'NH', 'NJ', 'NY', 'PA', 
                 'RI', 'VA', 'VT', 'WV']
    nox_state, nox_lat, nox_lon = AQSCEMSobs.cems_specifystates_dailymean(
            '/Volumes/GAIGEKERR/emissions/CEMS/', states_ab, sampling_months)
    nox = nox_state['%d-0%d-01'%(year, sampling_months[0]):
                    '%d-0%d-31'%(year, sampling_months[-1])].values
    # select year of interest
    yearpos = np.where(np.array(years) == year)[0][0]
    # data in year of interest
    castnet_t2m_ty = castnet_t2m[(yearpos*92):((yearpos+1)*92)]
    castnet_o3_ty = castnet_o3[(yearpos*92):((yearpos+1)*92)]
    # initialize figure, axes
    fig = plt.figure(figsize = (7, 5.5))
    ax1 = plt.subplot2grid((3, 4), (0, 0), colspan = 4)
    ax2 = plt.subplot2grid((3, 4), (2, 0), colspan = 4)
    ax3 = plt.subplot2grid((3, 4), (1, 0), colspan = 4)
    # MERRA-2 2-meter temperature plot
    ax1.set_title('(a)', fontsize = 16, x = 0.07, y = 1.03)
    ax1.plot(castnet_t2m_ty, lw = 2., color = '#ff7f00')
    ax1.set_xlim([0, len(castnet_t2m_ty) - 1])
    ax1.set_xticks([0, 30, 61])
    ax1.set_xticklabels([''])
    for t in ax1.get_yticklabels():
        t.set_fontsize(12)        
    ax1.get_yaxis().set_label_coords(-0.1, 0.5)
    ax1.set_ylabel('T [K]', fontsize = 16)
    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_ticks_position('both')
    # CASTNet/GMI O3
    ax2.set_title('(c)', fontsize = 16, x = 0.07, y = 1.03)
    ax2.plot(castnet_o3_ty, lw = 2., color = '#999999', label = 'CASTNet')
    ax2.plot(emiss_o3_neus[yearpos*92:(yearpos+1)*92], '-', 
            color = pollutants_constants.COLOR_EMISSIONS, 
            label = '$\mathregular{+}$AEmissions')
    ax2.plot(mr2_o3_neus[yearpos*92:(yearpos+1)*92],
            '-', color = pollutants_constants.COLOR_CHEMISTRY, 
            label = '$\mathregular{+}$Chemistry')
    ax2.plot(dat_o3_neus[yearpos*92:(yearpos+1)*92], '-', 
            color = pollutants_constants.COLOR_TRANSPORT, label = 'Transport')
    ax2.set_xlim([0, len(castnet_o3_ty) - 1])
    ax2.set_xticks([0, 30, 61])    
    ax2.set_xticklabels(['1 June %d' %year, '1 July', '1 Aug'], fontsize = 12)
    ax2.set_yticks([30, 45, 60])
    for t in ax2.get_yticklabels():
        t.set_fontsize(12)
    ax2.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    ax2.get_yaxis().set_label_coords(-0.1, 0.5)
    ax2.xaxis.set_ticks_position('both')
    ax2.yaxis.set_ticks_position('both')   
    # CEMS NOx plot
    ax3.set_title('(b)', fontsize = 16, x = 0.07, y = 1.03)
    ax3.plot(nox, lw = 2., color = '#377eb8')
    ax3.set_xlim([0, len(nox) - 1])
    ax3.set_xticks([0, 30, 61])
    ax3.set_xticklabels([''])    
    ax3.set_ylabel('NO$_{x}$ [tons]', fontsize = 16)
    for t in ax3.get_yticklabels():
        t.set_fontsize(12)
    ax3.get_yaxis().set_label_coords(-0.1, 0.5)
    pc = (nox - np.mean(nox))/np.mean(nox) * 100.
    ax3t = ax3.twinx()
    ax3t.plot(pc, lw = 1.0)
    for t in ax3t.get_yticklabels():
        t.set_fontsize(12)    
    ax3t.set_ylabel('Change [%]', rotation = 270, fontsize = 16)   
    ax3t.get_yaxis().set_label_coords(1.12, 0.48)
    ax3.xaxis.set_ticks_position('both')
    plt.subplots_adjust(bottom=0.2, hspace=0.5)
    # reorder and place legend
    handles, labels = ax2.get_legend_handles_labels()
    handles = [handles[0], handles[3], handles[2], handles[1]]
    labels = [labels[0], labels[3], labels[2], labels[1]]
    ax2.legend(handles,labels, bbox_to_anchor = (0.5, -.8), loc = 'center', 
              ncol = 2, frameon = False, fontsize = 16)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'timeseries_t2m_castneto3_cemsnox_%d_%s.eps' %(year, region), 
                dpi = 300)    
    print('For year %d...' %year)
    print('T2m-O3 correlation = %.3f' %(np.corrcoef(castnet_t2m_ty, castnet_o3_ty)[0, 1]))
    print('T2m-NOx correlation = %.3f' %(np.corrcoef(castnet_t2m_ty, nox)[0, 1]))
    print('O3-NOx correlation = %.3f' %(np.corrcoef(castnet_o3_ty, nox)[0, 1])) 
    print('Observed O3-climate penalty = %.3f ppbv K-1' %(np.polyfit(castnet_t2m_ty, castnet_o3_ty, 1)[0])) 
    return 
# # # # # # # # # # # # #  
def map_r2o3t2m_conus(lat_castnet, lon_castnet, r_castnet, gmi_lat, gmi_lon, 
    r, case): 
    """function plots maps of the coefficient of determination ("r-squared")
    calculated between O3 and MERRA2 2-meter temperatures from GMI and CASTNet.
    
    Parameters
    ----------    
    lat_castnet : list
        CASTNet latitude coordinates, units of degrees north, [no. CASTNet
        stations,]
    lon_castnet : list
        CASTNet longitude coordinates, units of degrees north, [no. CASTNet
        stations,]    
    r_castnet : list
        Pearson correlation coefficients between MERRA-2 2 meter temperatures
        and CASTNet ozone at each CASTNet site, [no. CASTNet stations,]
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
    do3dt2m : numpy.ndarray    
        The O3-T2m sensitivity at each GMI grid cell [lat, lon]   
    r : numpy.ndarray
        Pearson correlation coefficients between MERRA-2 2 meter temperatures
        and GMI ozone at each grid cell [lat, lon]
    case : str
        Model simulation name for output file

    Returns
    ----------     
    None
    """
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.colorbar import ColorbarBase
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    llcrnrlon = -126.
    llcrnrlat = 24.0
    urcrnrlon = -66.3
    urcrnrlat = 50.
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
        llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat, 
        resolution = 'h', area_thresh = 1000)
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)    
    fig = plt.figure(figsize = (5.5, 5))
    ax = plt.subplot2grid((1, 1), (0, 0))
    vmin = 0; vmax = 0.7
    cmap = plt.get_cmap('PuBu', 7)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 8, endpoint = True)
    m.contourf(x, y, r**2, clevs, cmap = cmap, extend = 'max')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    # superimpose O3-T2m sensitivity from observations
    x_castnet, y_castnet = m(lon_castnet, lat_castnet)
    m.scatter(x_castnet, y_castnet, c = np.array(r_castnet)**2, s = 18, 
              cmap = cmap, vmin = vmin, vmax = vmax, zorder = 30, 
              linewidth = 1., edgecolor = 'orange')    
    fill_oceans(ax, m)    
    outline_region(ax, m, pollutants_constants.NORTHEAST_STATES)    
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'max')
    cb.set_label(label = r'$r^{\mathregular{2}}$(O$_{\mathregular{3}}$' + 
                 ', T$_{}$)', size = 16)    
    cb.ax.tick_params(labelsize = 12)
    cb.set_ticks(np.linspace(vmin, vmax, 8)) 
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'map_r2o3t2m_conus_%s.eps' %case, dpi = 300)
    return 
# # # # # # # # # # # # #    
def timeseries_mr2o3dato3t2m_atpoint(t2m_overpass, o3, dat_o3, gmi_lat, 
    gmi_lon):
    """since there are differences in the O3-temperature (T) correlations 
    across the country, function plots timeseries of O3 from Transport and 
    + Chemistry simulations and co-located T. Even if correlations between O3 
    and T and the O3-climate sensitivity are low; however, the difference
    between the Transport and + Chemistry simulations still can give 
    information about the role of transport and deposition on O3 variability.
    Function also calculates relevant correlation coefficients.
    
    Parameters
    ----------  
    t2m_atoverpass : numpy.ndarry
        The MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]    
    o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from MR2/+ Chemistry 
        simulation, units of volume mixing ratio, [time, lat, lon]        
    dat_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from Diurnal-AvgT/
        Transport simulation, units of volume mixing ratio, [time, lat, lon]      
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
        
    Returns
    ----------     
    None         
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # indices of grid cell in South
    slon_idx_south = 33
    slat_idx_south = 8 + 2
    print('(%.3f, %.3f) is chosen coordinate point in South...' 
          %(gmi_lat[slat_idx_south], gmi_lon[slon_idx_south]))   
    # select values at point
    o3_atpoint = o3[184:, slat_idx_south, slat_idx_south] * 1e9
    dato3_atpoint = dat_o3[184:, slat_idx_south, slat_idx_south] * 1e9
    t2m_atpoint = t2m_overpass[184:, slat_idx_south, slat_idx_south]
    # calculate correlation coefficients/sensitivity
    print('r(Transport O3, + Chemistry O3) = %.3f' 
          %(np.corrcoef(o3_atpoint, dato3_atpoint)[0, 1]))
    print('r(+ Chemistry O3, T2m) = %.3f' 
          %(np.corrcoef(o3_atpoint, t2m_atpoint)[0, 1]))
    print('d(+ Chemistry O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_atpoint, o3_atpoint, 1)[0]))
    print('d(Transport O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_atpoint, dato3_atpoint, 1)[0]))
    print('Mean/Std Transport O3 = %.3f, %.3f' 
          %(np.mean(dato3_atpoint), np.std(dato3_atpoint)))
    print('Mean + Chemistry O3 = %.3f, %.3f' 
          %(np.mean(o3_atpoint), np.std(o3_atpoint)))
    slon_idx_north = 39
    slat_idx_north = 18
    print('(%.3f, %.3f) is chosen coordinate point in North...' 
          %(gmi_lat[slat_idx_north], gmi_lon[slon_idx_north]))   
    # select values at point
    o3_atpoint_north = o3[184:, slat_idx_north, slon_idx_north] * 1e9
    dato3_atpoint_north = dat_o3[184:, slat_idx_north, slon_idx_north] * 1e9
    t2m_atpoint_north = t2m_overpass[184:, slat_idx_north, slon_idx_north]
    # calculate correlation coefficients/sensitivity
    print('r(Transport O3, + Chemistry O3) = %.3f' 
          %(np.corrcoef(o3_atpoint_north, dato3_atpoint_north)[0, 1]))
    print('r(+ Chemistry O3, T2m) = %.3f' 
          %(np.corrcoef(o3_atpoint_north, t2m_atpoint_north)[0, 1]))
    print('d(+ Chemistry O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_atpoint_north, o3_atpoint_north, 1)[0]))
    print('d(Transport O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_atpoint_north, dato3_atpoint_north, 1)[0]))
    print('Mean/Std Transport O3 = %.3f, %.3f' 
          %(np.mean(dato3_atpoint_north), np.std(dato3_atpoint_north)))
    print('Mean + Chemistry O3 = %.3f, %.3f' 
          %(np.mean(o3_atpoint_north), np.std(o3_atpoint_north)))
    # initialize figure, axis 
    fig = plt.figure()
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan = 2)
    ax2 = plt.subplot2grid((2, 2), (1, 0), colspan = 2)
    # first subplot, O3 and T timeseries for a point that has high O3-T 
    # correlations
    ax1.set_title('(a)', ha = 'left', fontsize = 16, x = 0.02, y = 0.76)    
    lns1 = ax1.plot(o3_atpoint_north, lw = 2., 
                   c = pollutants_constants.COLOR_CHEMISTRY, 
                   label = 'Transport', zorder = 10)
    lns2 = ax1.plot(dato3_atpoint_north, lw = 2., 
                   c = pollutants_constants.COLOR_TRANSPORT,
                   label = '+ Chemistry', zorder = 5)
    ax1.set_xlim([0, len(o3_atpoint_north) - 1])
    for t in ax1.get_yticklabels():
        t.set_fontsize(12)   
    ax1.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    ax1b = ax1.twinx()
    lns3 = ax1b.plot(t2m_atpoint_north, lw = 2., color = '#ff7f00',
             label = 'T', zorder = 1)
    ax1b.set_xlim([0, len(t2m_atpoint) - 1])
    ax1b.set_xticks([0, 14, 29, 45, 61, 76])
    ax1b.set_xticklabels([])
    ax1.tick_params(right = False, left = True, top = True, bottom = False,
                    labelbottom = False, labeltop = True)
    ax1b.set_ylabel('T [K]', rotation = 270, fontsize = 16, color = '#ff7f00')
    ax1b.get_yaxis().set_label_coords(1.15, 0.53)
    for t in ax1b.get_yticklabels():
        t.set_fontsize(12)    
        t.set_color('#ff7f00')
    # second subplot, O3 and T timeseries for a point that has low O3-T 
    # correlations
    ax2.set_title('(b)', ha = 'left', fontsize = 16, x = 0.02, y = 0.76)    
    lns1 = ax2.plot(o3_atpoint, lw = 2., 
                   c = pollutants_constants.COLOR_CHEMISTRY, zorder = 10)
    lns2 = ax2.plot(dato3_atpoint, lw = 2., 
                   c = pollutants_constants.COLOR_TRANSPORT, zorder = 5)
    ax2.set_xlim([0, len(o3_atpoint) - 1])
    for t in ax2.get_yticklabels():
        t.set_fontsize(12)   
    ax2.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    ax2b = ax2.twinx()
    lns3 = ax2b.plot(t2m_atpoint, lw = 2., color = '#ff7f00', zorder = 1)
    ax2b.set_xlim([0, len(t2m_atpoint) - 1])
    ax2b.set_xticks([0, 14, 29, 45, 61, 76])
    ax2b.set_xticklabels(['1 June 2010', '', '1 July', '', '1 August', ''], 
                         fontsize = 16)
    ax2b.set_ylabel('T [K]', rotation = 270, fontsize = 16, color = '#ff7f00')
    ax2b.get_yaxis().set_label_coords(1.15, 0.53)
    for t in ax2.get_xticklabels():
        t.set_fontsize(12)
    for t in ax2b.get_yticklabels():
        t.set_fontsize(12)    
        t.set_color('#ff7f00')
    # add inset map showing location of cells
    left, bottom, width, height = [0.6, 0.4, 0.2, 0.22]
    axi = fig.add_axes([left, bottom, width, height])        
    m_small = Basemap(projection = 'merc', llcrnrlon = -93.5, 
                      llcrnrlat = 24., urcrnrlon = -65.3, 
                      urcrnrlat = 50., resolution = 'h', area_thresh = 1000)
    m_small.drawstates(color = 'k', linewidth = 0.5)
    m_small.drawcountries(color = 'k', linewidth = 1.0)
    m_small.drawcoastlines(color = 'k', linewidth = 1.0)   
    fill_oceans(axi, m_small)   
    x, y = m_small(gmi_lon[slon_idx_south], gmi_lat[slat_idx_south])
    m_small.scatter(x, y, 30, color = 'r', marker = '.', edgecolor = 'r', 
              zorder = 20)
    x, y = m_small(gmi_lon[slon_idx_north], gmi_lat[slat_idx_north])
    m_small.scatter(x, y, 30, color = 'r', marker = '.', edgecolor = 'r', 
              zorder = 20)
    # add legend
    leg = ax1.legend(loc = 9, bbox_to_anchor = (0.5, 1.35), fontsize = 16,
                    ncol = 2)
    leg.get_frame().set_linewidth(0.0)
    plt.subplots_adjust(right = 0.88)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'timeseries_mr2o3dato3t_atpoint_2010.eps', dpi = 300)
    return
# # # # # # # # # # # # #    
def map_simulationschematic(gmi_lat, gmi_lon): 
    """plot location of CEMS facilities in the Eastern U.S. (east of 
    Mississippi River). The size and colors of the markers denoting their 
    locations is based on dividing their cumulative emissions over the
    measuring period 2000 - 2014 into percentiles. The shaded regions 
    corresponds to the CTM domain that was perturbed in the + Emissions 
    simulation. 
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Polygon
    sys.path.append('/Users/ghkerr/phd/')
    from geo_idx import geo_idx
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/emissions/')
    import cems_nei
    # # # #
    def get_marker_color(amount_array):
        """find scatterplot marker color and size based on value with colormap
        from http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=5
    
        Parameters
        ----------
        amount_array : numpy.array
            the average yearly sum of NOx emissions (tons) from each CEMS
            emitting facility in the focus region
            
        Returns
        ----------
        sizes : list
            marker sizes
        colors : list
            hex color code for markers
        """
        sizes = []
        colors = []        
        for value in amount_array: 
            if (value < 25.618500000000004):
                sizes.append(18)
                colors.append('#a6cee3')
            if (value >= 25.618500000000004) and (value < 134.91200000000012):
                sizes.append(26)
                colors.append('#1f78b4')
            if (value >= 134.91200000000012) and (value < 1065.5449999999996):
                sizes.append(34)
                colors.append('#b2df8a')                          
            if (value >= 1065.5449999999996) and (value < 3282.3063999999977):
                sizes.append(42)
                colors.append('#33a02c')                                  
            if (value >= 3282.3063999999977):
                sizes.append(50)
                colors.append('#e31a1c')                           
        return sizes, colors
    # # # # 
    # initialize figure, axis
    fig = plt.figure(figsize = (8, 10))
    ax = plt.subplot2grid((1, 1), (0, 0))
    # focus region map 
    m = Basemap(projection = 'merc', llcrnrlon = -93.5, 
                llcrnrlat = 24., urcrnrlon = -65.3, 
                urcrnrlat = 50., resolution = 'h', area_thresh = 1000)
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
#    m.fillcontinents(color = '#f9f6d8', lake_color = '#dcf0fa')
    fill_oceans(ax, m)    
    outline_region(ax, m, pollutants_constants.NORTHEAST_STATES)    
    # three bounding boxes: 
    # 1. 47˚N, -80˚E, 28˚N, -92˚E
    bb1t = geo_idx(47., gmi_lat)
    bb1r = geo_idx(-80., gmi_lon)
    bb1b = geo_idx(28., gmi_lat)
    bb1l = geo_idx(-92., gmi_lon)
    # 2. 48˚N, -66˚E, 41˚N, -80˚N
    bb2t = geo_idx(48., gmi_lat)
    bb2r = geo_idx(-66., gmi_lon)
    bb2b = geo_idx(41., gmi_lat)
    bb2l = geo_idx(-80., gmi_lon)
    # 3. 41˚N, -74˚E, 28˚N, -80˚E
    bb3t = geo_idx(41., gmi_lat)
    bb3r = geo_idx(-74., gmi_lon)
    bb3b = geo_idx(28., gmi_lat)
    bb3l = geo_idx(-80., gmi_lon)
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)
    # bounding box #1
    x1, y1 = m(gmi_lon[bb1l], gmi_lat[bb1b])
    x2, y2 = m(gmi_lon[bb1l], gmi_lat[bb1t])
    x3, y3 = m(gmi_lon[bb1r], gmi_lat[bb1t])
    x4, y4 = m(gmi_lon[bb1r], gmi_lat[bb1b])
    poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)], facecolor = '#ff7f00', 
                    alpha = 0.6, linewidth = 0)
    plt.gca().add_patch(poly)
    # bounding box #2
    x1, y1 = m(gmi_lon[bb2l], gmi_lat[bb2b])
    x2, y2 = m(gmi_lon[bb2l], gmi_lat[bb2t])
    x3, y3 = m(gmi_lon[bb2r], gmi_lat[bb2t])
    x4, y4 = m(gmi_lon[bb2r], gmi_lat[bb2b])
    poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)], facecolor = '#ff7f00', 
                    alpha = 0.6, linewidth = 0)
    plt.gca().add_patch(poly)
    # bounding box #3 
    x1, y1 = m(gmi_lon[bb3l], gmi_lat[bb3b])
    x2, y2 = m(gmi_lon[bb3l], gmi_lat[bb3t])
    x3, y3 = m(gmi_lon[bb3r], gmi_lat[bb3t])
    x4, y4 = m(gmi_lon[bb3r], gmi_lat[bb3b])
    poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)], facecolor = '#ff7f00', 
                    alpha = 0.6, linewidth = 0)
    plt.gca().add_patch(poly)
    # import emissions information
    facilityinfo = cems_nei.facility_attributes('/Volumes/GAIGEKERR/emissions/CEMS/')
    emissions_nomean = cems_nei.read_emissions_nomean_easternus('/Volumes/GAIGEKERR/emissions/CEMS/')
    # given JJA CEMS measurements for measuring period, find unique stack IDs
    unique_orispl = np.unique(emissions_nomean['Facility ID (ORISPL)'].values)
    # look up unique stack IDs in facility attribute files
    facilityinfo_reduced = facilityinfo.loc[facilityinfo['Facility ID (ORISPL)'
                                                         ].isin(unique_orispl)]
    # group by Facility ID and find unique latitude and longtiudes at inidividual
    # facilities; data is screwy and, for instance, some years don't have 
    # values for facility latitude or longitude (nan) or unreasonable values.
    # Filter away these measurements
    facilitylon = facilityinfo_reduced.groupby('Facility ID (ORISPL)')['Facility Longitude'].unique()
    facilitylat = facilityinfo_reduced.groupby('Facility ID (ORISPL)')['Facility Latitude'].unique()
    # group measurements by stack ID, find sum of total NOx emissions over the 
    # entire measuring period
    summed_nox = emissions_nomean.groupby('Facility ID (ORISPL)')['NOx (tons)'].sum()
    # total NOx emissions over the measuring period will be the size of the 
    # scatterpoints on map; ensure that the stack IDs of the stacks' locations
    # (lon/lat) in the same order as the IDs of total NOx emissions
    if summed_nox.index.difference(facilitylat.index).shape != (0,):
        print('mismatch between NOx sum and stack location!')
    # find marker sizes, colors for plot
    sizes, colors = get_marker_color(summed_nox.values)
    xcems, ycems = m(np.concatenate(facilitylon.values), 
             np.concatenate(facilitylat.values))      
    m.scatter(xcems, ycems, s = sizes, color = colors, marker = 'o', 
              edgecolor = colors, zorder = 30)
    # legend for a scatter plot using a proxy artists 
    # see http://matplotlib.sourceforge.net/users/legend_guide.html#using-proxy-artist
    circle1 = Line2D(range(1), range(1), color = 'w', marker = 'o', 
        markersize = np.log(18)**1.5, zorder = 50, markerfacecolor = '#a6cee3', 
        markeredgecolor = '#a6cee3', label = '$\mathregular{\Sigma}$$\,' + 
        '$NO$_x$$\,$<$\,$25$^{\mathregular{th}}$')
    circle2 = Line2D(range(1), range(1), color = 'w', marker = 'o', 
        markersize = np.log(26)**1.5, markerfacecolor = '#1f78b4', 
        markeredgecolor = '#1f78b4', zorder = 50, label = '25$^{\mathregular'+
        '{th}}$$\,$$\mathregular{\leq}$$\,$$\mathregular{\Sigma}$$\,$NO$_x$'+
        '$\,$<$\,$50$^{\mathregular{th}}$')
    circle3 = Line2D(range(1), range(1), color = 'w', marker = 'o', 
        markersize = np.log(34)**1.5, zorder = 50, markerfacecolor = '#b2df8a', 
        markeredgecolor = '#b2df8a', label = '50$^{\mathregular{th}}$$\,$'+
        '$\mathregular{\leq}$$\,$$\mathregular{\Sigma}$$\,$NO$_x$$\,$<$\,$'+
        '75$^{\mathregular{th}}$')
    circle4 = Line2D(range(1), range(1), color = 'w', marker = 'o', 
        markersize = np.log(42)**1.5, zorder = 50, markerfacecolor = '#33a02c', 
        markeredgecolor = '#33a02c', label = '75$^{\mathregular{th}}$$\,$'+
        '$\mathregular{\leq}$$\,$$\mathregular{\Sigma}$ NO$_x$$\,$<$\,$90'+
        '$^{\mathregular{th}}$')
    circle5 = Line2D(range(1), range(1), color = 'w', marker = 'o', 
        markersize = np.log(50)**1.5, zorder = 50, markerfacecolor = '#e31a1c', 
        markeredgecolor = '#e31a1c', label = '$\mathregular{\Sigma}$$\,'+
        '$NO$_x$$\,$$\mathregular{\geq}$$\,$90$^{\mathregular{th}}$')
    # add legend 
    leg = ax.legend(loc = 9, bbox_to_anchor = (0.5, 1.2),
                    handles = [circle1, circle2, circle3, circle4, circle5], 
                    ncol = 2, fontsize = 16, numpoints = 1, facecolor = 'w')
    leg.get_frame().set_linewidth(0.0)
    plt.subplots_adjust(right = 0.85)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 'map_simulationschematic.eps', 
                dpi = 300)
    return
# # # # # # # # # # # # #    
def compare_regional_averaging(neus, t_castnet, o3_castnet, castnet_t2m_neus, 
    castnet_o3_neus, sites_castnet, neus_castnet, t2m_overpass, emiss_o3, 
    emiss_t2m_neus, emiss_o3_neus, emiss_sens, do3dt2m_castnet):
    """examine two different ways to calculate regional averages. The first
    method ("regional quantities") uses regionally-averaged primary fields 
    (i.e. O3, T) to calculate secondary quantities like the O3-climate penalty. 
    The second method ("regional average") calculates the secondary quantity at
    every grid cell and thereafter averages the secondary quantity over the 
    region. This is done for both CASTNet and GMI for the O3-climate penalty,
    the standard deviation, and the O3 increase on hot days. 
    
    Parameters
    ----------  
    neus : numpy.ndarray
        CTM grid where grid cells within NEUS have a value of 1 and grid 
        cells not in region have a value of NaN, [lat, lon]
    t_castnet : list
        Daily 1300 hours local time MERRA-2 2-meter temperatures at each 
        CASTNet site, units of K, [no. sites,]        
    o3_castnet : list
        Daily 1300 hours local time O3 at each CASTNet site, units of ppbv, 
        [no. sites,]
    castnet_t2m_neus : numpy.ndarray
        Daily regionally-averaged 2-meter temperatures co-located with CASTNet 
        stations, [time,]      
    castnet_o3_neus : numpy.ndarray
        Daily regionally-averaged O3 from CASTNet, [time,]         
    sites_castnet : list
        CASTNet site names, [no. sites,]
    neus_castnet : list
        Site IDs for CASTNet sites in region, [no. sites in region,]    
    t2m_overpass : numpy.ndarray
        MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]
    emiss_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, for model case 
        GHKerr-DailyEmiss, units of volume mixing ratio, [time, lat, lon] 
    emiss_t2m_neus : numpy.ndarray
        Daily regionally-averaged T2m from + Emissions simulation, [time,]    
    emiss_o3_neus : numpy.ndarray
        Daily regionally-averaged O3 from + Emissions simulation, [time,]          
    emiss_sens : numpy.ndarray    
        The O3-climate penalty at each GMI grid cell from the + Emissions 
        simulation, units of ppbv K^-1, [lat, lon]  
    do3dt2m_castnet : list
        The O3-T2m sensitivity (slope of linear regression) at each CASTNet
        site, [no. CASTNet stations,]     
        
    Returns
    ----------     
    None     
    """
    import numpy as np
    # rq stands for "regional quantities," i.e. quantity has been calculated 
    # using primary quantities that have been regionally-averaged. For example, 
    # the O3-climate penalty calculated with regional quantities means that it 
    # was calculated with regionally-averaged T and O3. ra stands for "regional
    # averaged" meaning that the final product has been averaged over the 
    # region. For the O3-climate penalty example, this would mean that the O3-
    # climate penalty was calculated at local areas/points/sites within a 
    # region using local quantities and then the penalty itself was averaged. 
    # # # # for O3-climate penalty
    # GMI
    emiss_m_rq = np.polyfit(emiss_t2m_neus, emiss_o3_neus, 1)[0]
    emiss_m_ra = np.nanmean(neus*emiss_sens)
    print('Regional quantity: GMI penalty = %.4f ppbv K-1' %emiss_m_rq)
    print('Regional average: GMI penalty = %.4f ppbv K-1' %emiss_m_ra)    
    where_neus = np.in1d(sites_castnet, neus_castnet)
    where_neus = np.where(where_neus == True)[0]
    # CASTNet
    castnet_m_rq = np.polyfit(castnet_t2m_neus, castnet_o3_neus, 1)[0]
    castnet_m_ra = np.array(do3dt2m_castnet)[where_neus].mean()
    print('Regional quantity: CASTNet penalty = %.4f ppbv K-1' %castnet_m_rq)
    print('Regional average: CASTNet penalty = %.4f ppbv K-1' %castnet_m_ra)      
    # # # # for variability
    # GMI
    emiss_std_rq = np.std(emiss_o3_neus)
    emiss_std_ra = neus*(np.std(emiss_o3, axis = 0) * 1e9)
    emiss_std_ra = np.nanmean(emiss_std_ra)
    print('Regional quantity: GMI standard deviation = %.4f ppbv'
          %emiss_std_rq)
    print('Regional average: GMI standard deviation = %.4f ppbv'
          %emiss_std_ra)      
    # CASTNet 
    castnet_std_rq = np.std(castnet_o3_neus)
    castnet_std_ra = np.nanstd(np.array(o3_castnet)[where_neus], axis = 1).mean()
    print('Regional quantity: CASTNet standard deviation = %.4f ppbv' 
          %castnet_std_rq)
    print('Regional average: CASTNet standard deviation = %.4f ppbv' 
          %castnet_std_ra)       
    # # # # for hot/cold days 
    # GMI
    # find mean O3 in each simulation on days where 2-meter temperatures exceed
    # the 90th percentile
    emiss_o3_p90 = o3_givent2mpercentile(emiss_o3, t2m_overpass, 90, 'greater')
    # find median O3 values in each simulation and convert from volume mixing
    # ratio to ppbv
    emiss_o3_med = np.median(emiss_o3, axis = 0) * 1e9
    delta_emiss_hot = emiss_o3_p90 - emiss_o3_med
    emiss_o3_p90_rq = emiss_o3_neus[np.where(emiss_t2m_neus >
        np.percentile(emiss_t2m_neus, 90))[0]]
    emiss_o3_med_rq = np.median(emiss_o3_neus)
    emiss_hot_rq = emiss_o3_p90_rq.mean() - emiss_o3_med_rq
    emiss_hot_ra = np.nanmean(delta_emiss_hot*neus)
    print('Regional quantity: GMI hot day increase = %.4f ppbv' 
          %emiss_hot_rq)
    print('Regional average: GMI hot day increase = %.4f ppbv' 
          %emiss_hot_ra)    
    # CASTNet
    castnet_o3_med = np.median(castnet_o3_neus)
    castnet_o3_p90 = castnet_o3_neus[np.where(castnet_t2m_neus > 
                              np.percentile(castnet_t2m_neus, 90))[0]]
    castnet_hot_ra = []
    for o3_as, t2m_as in zip(np.array(o3_castnet)[where_neus], 
                             np.array(t_castnet)[where_neus]):
        o3_as_p90 = np.nanmean(o3_as[np.where(t2m_as > np.percentile(t2m_as, 90))[0]])
        o3_as_med = np.nanmedian(o3_as)
        castnet_hot_ra.append(o3_as_p90 - o3_as_med)
    castnet_hot_rq = castnet_o3_p90.mean() - castnet_o3_med
    castnet_hot_ra = np.mean(castnet_hot_ra)
    print('Regional quantity: CASTNet hot day increase = %.4f ppbv' 
          %castnet_hot_rq)
    print('Regional average: CASTNet hot day increase = %.4f ppbv' 
          %castnet_hot_ra)           
    return 
# # # # # # # # # # # # #    
def NO_inventory_atpoint(t2m_overpass, ilat, ilon, year, gmi_lat, gmi_lon): 
    """given a latitude, longitude coordinate and year, function finds the 
    closest grid cell in the emission inventory of the GMI CTM and plots the 
    (1) monthly mean values from the control run and (2) the scaled (daily-
    varying) values from the emissions-temperature sensitivity simulations. 
    Here "values" are NO flux per grid box from the fossil fuel sector for 
    summer months. Also plotted are MERRA2 2-meter temperatures from the grid 
    cell nearest the given latitude, longitude coordinate. A small inset map 
    showing the location of the grid cell is also plotted. n.b. if the 
    inventory/temperature over Baltimore is desired, ilat = 39.2904 and 
    ilon = 283.387. 
    
    Parameters
    ----------  
    t2m_overpass : numpy.ndarray
        MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]
    ilat : float
        Latitude of interest
    ilon : float
        Longitude (0 - 360˚) of interest
    year : int
        Year of interest
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
        
    Returns
    ----------
    None      
    """
    import numpy as np
    from netCDF4 import Dataset
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap    
    from matplotlib.lines import Line2D
    from matplotlib.patches import Polygon
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    # load perturbed emissions inventory 
    month = 'JUN'
    dataset = Dataset(pollutants_constants.PATH_GMIEMISSIONS + 
                     'hindcast_IAVanthQFED.daily.1x1.25_%d_%s_NOscaling.nc' 
                     %(year, month), 'r')
    species_raw = dataset['species'][:]
    species = []
    for row in np.arange(0, np.shape(species_raw)[0]):
        temp = []
        for element in species_raw[row]:
            temp.append(element.decode('UTF-8'))
        si = ''.join(temp)
        si = si.rstrip()
        species.append(si[:])      
    # fetch perturbed fields and dimensions
    lat = dataset['latitude_dim'][:]
    lon = dataset['longitude_dim'][:]
    NO_ff_prime = dataset['emiss_2d'][:, np.where(np.array(species) == 'NO_ff')[0][0], :, :]
    # total surface NO
    NO_total_prime_jun = NO_ff_prime
    bmore_lat = geo_idx(ilat, lat)
    bmore_lon = geo_idx(ilon, lon)
    NO_total_prime_jun = NO_total_prime_jun[:, bmore_lat, bmore_lon]
    # # # #
    month = 'JUL'
    dataset = Dataset(pollutants_constants.PATH_GMIEMISSIONS + 
                     'hindcast_IAVanthQFED.daily.1x1.25_%d_%s_NOscaling.nc' 
                     %(year, month), 'r')    
    NO_ff_prime = dataset['emiss_2d'][:, np.where(np.array(species) == 'NO_ff')[0][0], :, :]
    NO_total_prime_jul = NO_ff_prime
    NO_total_prime_jul = NO_total_prime_jul[:, bmore_lat, bmore_lon]
    # # # # August of summer of interest
    month = 'AUG'
    dataset = Dataset(pollutants_constants.PATH_GMIEMISSIONS + 
                     'hindcast_IAVanthQFED.daily.1x1.25_%d_%s_NOscaling.nc' 
                     %(year, month), 'r')     
    NO_ff_prime = dataset['emiss_2d'][:, np.where(np.array(species) == 'NO_ff')[0][0], :, :]
    NO_total_prime_aug = NO_ff_prime
    NO_total_prime_aug = NO_total_prime_aug[:, bmore_lat, bmore_lon]
    del bmore_lat, bmore_lon
    # find 2-meter temperatures
    year_where = np.where(np.array([2008, 2009, 2010]) == year)[0][0]
    bmore_lat = geo_idx(ilat, gmi_lat)
    bmore_lon = geo_idx(ilon - 360., gmi_lon)
    t2m_bmore = t2m_overpass[92*year_where:(92*(year_where + 1) + 1), 
                             bmore_lat, bmore_lon]
    # initialize figure, axis
    fig = plt.figure()
    ax = plt.subplot2grid((1, 3), (0, 0), rowspan = 1, colspan = 3)
    # twin axis, add temperature 
    ax2 = ax.twinx()
    ax2.plot(t2m_bmore, linewidth = 2., color = '#ff7f00', clip_on = True, 
             zorder = 2)
    ax.plot(np.concatenate([NO_total_prime_jun, NO_total_prime_jul, NO_total_prime_aug]), 
            linewidth = 1.5, clip_on = True, label = '$\mathregular{' +
            r'\partial}$Emissions $\mathregular{\partial}$T$^{\mathregular{-1}}\neq$0', 
            color = '#377eb8', linestyle = ':')
    ax.plot(np.arange(0, 30, 1), np.repeat(np.mean(NO_total_prime_jun), 30), 
            linewidth = 2., clip_on = True,  color = '#377eb8', zorder = 6)
    ax.plot(np.arange(30, 61, 1), np.repeat(np.mean(NO_total_prime_jul), 31), 
            linewidth = 2., clip_on = True, color = '#377eb8', zorder = 6)
    ax.plot(np.arange(61, 92, 1), np.repeat(np.mean(NO_total_prime_aug), 31), 
            linewidth = 2., clip_on = True, label = '$\mathregular{' +
            r'\partial}$Emissions $\mathregular{\partial}$T$^{\mathregular{-1}}=$0', 
            color = '#377eb8')
    ax.set_xlim([0, 91])
    ax.set_xticks([0, 14, 30, 44, 61, 75])
    ax.set_xticklabels(['1 June %d' %year, '', '1 July', '', '1 Aug', ''], 
                       ha = 'center', fontsize = 12)
    # add legend and title
    leg = ax.legend(loc = 9, bbox_to_anchor = (0.5, 1.18), fontsize = 16,
                    ncol = 2)
    leg.get_frame().set_linewidth(0.0)
    ax.set_ylabel('Emissions [kg s$^{\mathregular{-1}}$'+ 
                        ' grid cell$^{\mathregular{-1}}$]', color = '#377eb8', 
                        fontsize = 16)
    ax2.set_ylabel('T [K]', color = '#ff7f00', fontsize = 16, rotation = 270)
    ax2.get_yaxis().set_label_coords(1.18, 0.50)
    for label in ax.get_yticklabels():
        label.set_fontsize(12)   
        label.set_color('#377eb8')
    for label in ax2.get_yticklabels():
        label.set_fontsize(12)    
        label.set_color('#ff7f00')
    # add inset axis
    left, bottom, width, height = [0.46, 0.16, 0.28, 0.28]
    axi = fig.add_axes([left, bottom, width, height])        
    m = Basemap(projection = 'merc', llcrnrlon = -93.5, 
                llcrnrlat = 24., urcrnrlon = -65.3, 
                urcrnrlat = 50., resolution = 'h', area_thresh = 1000)
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(axi, m)    
    # three bounding boxes: 
    # 1. 47˚N, -80˚E, 28˚N, -92˚E
    bb1t = geo_idx(47., gmi_lat)
    bb1r = geo_idx(-80., gmi_lon)
    bb1b = geo_idx(28., gmi_lat)
    bb1l = geo_idx(-92., gmi_lon)
    # 2. 48˚N, -66˚E, 41˚N, -80˚N
    bb2t = geo_idx(48., gmi_lat)
    bb2r = geo_idx(-66., gmi_lon)
    bb2b = geo_idx(41., gmi_lat)
    bb2l = geo_idx(-80., gmi_lon)
    # 3. 41˚N, -74˚E, 28˚N, -80˚E
    bb3t = geo_idx(41., gmi_lat)
    bb3r = geo_idx(-74., gmi_lon)
    bb3b = geo_idx(28., gmi_lat)
    bb3l = geo_idx(-80., gmi_lon)
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)
    # bounding box #1
    x1, y1 = m(gmi_lon[bb1l], gmi_lat[bb1b])
    x2, y2 = m(gmi_lon[bb1l], gmi_lat[bb1t])
    x3, y3 = m(gmi_lon[bb1r], gmi_lat[bb1t])
    x4, y4 = m(gmi_lon[bb1r], gmi_lat[bb1b])
    poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)], facecolor = '#ff7f00', 
                    alpha = 0.6, linewidth = 0)
    plt.gca().add_patch(poly)
    # bounding box #2
    x1, y1 = m(gmi_lon[bb2l], gmi_lat[bb2b])
    x2, y2 = m(gmi_lon[bb2l], gmi_lat[bb2t])
    x3, y3 = m(gmi_lon[bb2r], gmi_lat[bb2t])
    x4, y4 = m(gmi_lon[bb2r], gmi_lat[bb2b])
    poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)], facecolor = '#ff7f00', 
                    alpha = 0.6, linewidth = 0)
    plt.gca().add_patch(poly)
    # bounding box #3 
    x1, y1 = m(gmi_lon[bb3l], gmi_lat[bb3b])
    x2, y2 = m(gmi_lon[bb3l], gmi_lat[bb3t])
    x3, y3 = m(gmi_lon[bb3r], gmi_lat[bb3t])
    x4, y4 = m(gmi_lon[bb3r], gmi_lat[bb3b])
    poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)], facecolor = '#ff7f00', 
                    alpha = 0.6, linewidth = 0)
    plt.gca().add_patch(poly)
    # add point corresponding to grid cell 
    x, y = m(ilon - 360, ilat) 
    m.scatter(x, y, 20, color = 'w', marker = 's', edgecolor = 'w', zorder = 20)
    plt.subplots_adjust(right = 0.85)
    plt.savefig('/Users/ghkerr/phd/emissions/figs/' + 
                'NO_inventory_atpoint.eps', dpi = 300)
    return
# # # # # # # # # # # # #      
def scatter_inventorynonox_noxo3(t2m_overpass, mr2_no, mr2_no2, mr2_o3, 
    emiss_no, emiss_no2, emiss_o3, neus_states, gmi_lat, gmi_lon):
    """function opens output from Strode et al. [2015] EmFix and Std 
    simulations and the emissions inventories from these simulations along with 
    the emissions inventory from the + Emissions and + Chemistry simulations. 
    Thereafter we calculate regionally-averaged quantities in the
    Northeast. These results are compared with results from + Emissions and 
    + Chemistry simulations through two scatterplots and related regressions. 
    The first scatterplot is modeled NOx versus NO from the models' emissions
    inventories and thes second is modeled O3 versus modeled NOx. 

    Parameters
    ----------  
    t2m_overpass : numpy.ndarray
        MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]
    mr2_no : numpy.ndarray
        GMI CTM surface-level nitrogen oxide at overpass time from the 
        + Chemistry simulation, units of volume mixing ratio, [time, lat, lon]
    mr2_no2 : numpy.ndarray
        GMI CTM surface-level nitrogen dioxide at overpass time from the 
        + Chemistry simulation, units of volume mixing ratio, [time, lat, lon]            
    mr2_o3 : numpy.ndarray    
        GMI CTM surface-level ozone at overpass time from the + Chemistry 
        simulation, units of volume mixing ratio, [time, lat, lon]
    emiss_no : numpy.ndarray
        GMI CTM surface-level nitrogen oxide at overpass time from the 
        + Emissions simulation, units of volume mixing ratio, [time, lat, lon]
    emiss_no2 : numpy.ndarray
        GMI CTM surface-level nitrogen dioxide at overpass time from the 
        + Emissions simulation, units of volume mixing ratio, [time, lat, lon]        
    emiss_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from the + Emissions 
        simulation, units of volume mixing ratio, [time, lat, lon]
    neus_states : list
        State names (string format) in region      
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
                
    Returns
    ----------     
    std_inventory_daily : numpy.ndarray   
        Daily Northeast-averaged NO emissions from Strode et al. [2015] 
        emissions inventory Std simulation for the measuring period 2000 - 
        2010, units of kg grid cell-1 s-1, [time,]
    std_o3_neus : numpy.ndarray
        Daily Northeast-averaged O3 from Strode et al. [2015] 
        emissions inventory Std simulation for the measuring period 2000 - 
        2010, units of ppbv, [time,]
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/GMI/')
    import commensurability     
    # load Std
    (gmi_lat_c, gmi_lon_c, eta_c, times_c, std_co, std_no, std_no2, std_o3, 
     std_cloudfraction, std_gridboxheight) = \
    commensurability.open_overpass2('HindcastFFIgac2', [2000, 2001, 2002, 
        2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010])
    # load EmFix
    (gmi_lat_c, gmi_lon_c, eta_c, times_c, emfix_co, emfix_no, emfix_no2, emfix_o3, 
     std_cloudfraction, std_gridboxheight) = \
    commensurability.open_overpass2('Hindcast3Igac2', [2000, 2001, 2002, 
        2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010])
    gmi_lon_c = np.mod(gmi_lon_c - 180.0, 360.0) - 180.0
    # emissions inventory from + Emissions simulation
    emiss_no_inventory, lon_inventory, lat_inventory = \
     commensurability.open_perturbed_emissions()
    # emissions inventory from + Chemistry simulation 
    mr2_no_inventory, lon_inventory, lat_inventory = \
     commensurability.open_unperturbed_emissions('1x1.25_IAVanthGFED4')
    lon_inventory = np.mod(lon_inventory[0] - 180.0, 360.0) - 180.0
    # emissions inventory from Std/EmFix simulations
    std_no_inventory, lon_inventory_c, lat_inventory_c = \
     commensurability.open_unperturbed_emissions('2x2.5_IAVanthGFED3gcEF')
    lon_inventory_c = np.mod(lon_inventory_c[0] - 180.0, 360.0) - 180.0
    # find grid cells in Northeast for different simulations/inventories
    m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 24.0, 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'c', 
                area_thresh = 1000)
    neus_states = pollutants_constants.NORTHEAST_STATES
    neus = find_grid_in_region(m, neus_states, gmi_lat, gmi_lon)
    neus_c = find_grid_in_region(m, neus_states, gmi_lat_c, gmi_lon_c) 
    neus_inventory = find_grid_in_region(m, neus_states, lat_inventory[0], 
                                         lon_inventory)  
    neus_inventory_c = find_grid_in_region(m, neus_states, lat_inventory_c[0], 
                                           lon_inventory_c)       
     # # # find Northeast-averaged inventories, NOx, O3 
    # from Strode et al. [2015] simulations
    std_no_inventory_c_neus = np.array(std_no_inventory) * neus_inventory_c
    std_no_inventory_c_neus = np.nanmean(std_no_inventory_c_neus, axis = tuple((2, 3)))
    std_o3_neus = np.nanmean((std_o3 * neus_c), axis = tuple((1, 2))) * 1e9
    std_nox_neus = np.nanmean(((std_no + std_no2) * neus_c), 
                              axis = tuple((1, 2))) * 1e9
    emfix_o3_neus = np.nanmean((emfix_o3 * neus_c), axis = tuple((1, 2))) * 1e9
    emfix_nox_neus = np.nanmean(((emfix_no + emfix_no2) * neus_c), 
                              axis = tuple((1, 2))) * 1e9
    # from + Emissions/+ Chemistry simulations
    emiss_no_inventory_neus = np.vstack(emiss_no_inventory) * neus_inventory
    emiss_no_inventory_neus = np.nanmean(emiss_no_inventory_neus, axis = tuple((1, 2)))             
    mr2_no_inventory_neus = np.vstack(mr2_no_inventory) * neus_inventory
    mr2_no_inventory_neus = np.nanmean(mr2_no_inventory_neus, axis = tuple((1, 2)))    
    # since the inventory from MR2 only varies from month-to-month, repeat 
    # monthly values for every day in month
    mr2_inventory_daily = []
    for year in (mr2_no_inventory_neus[0:3], 
                 mr2_no_inventory_neus[3:6],
                 mr2_no_inventory_neus[6:]): 
        mr2_inventory_daily.append(np.repeat(year[0], 30))
        mr2_inventory_daily.append(np.repeat(year[1], 31))
        mr2_inventory_daily.append(np.repeat(year[2], 31))
    mr2_no_inventory_neus = np.hstack(mr2_inventory_daily)
    emiss_o3_neus = np.nanmean((emiss_o3 * neus), axis = tuple((1, 2))) * 1e9
    emiss_nox_neus = np.nanmean(((emiss_no + emiss_no2) * neus), 
                                axis = tuple((1, 2))) * 1e9
    mr2_o3_neus = np.nanmean((mr2_o3 * neus), axis = tuple((1, 2))) * 1e9
    mr2_nox_neus = np.nanmean(((mr2_no + mr2_no2) * neus), 
                              axis = tuple((1, 2))) * 1e9
    # find regionally-averaged temperatures 
    t2m_neus = np.nanmean((t2m_overpass * neus), axis = tuple((1, 2)))
    t2m_hot = np.where(t2m_neus > np.percentile(t2m_neus, 90))[0]
    t2m_cold = np.where(t2m_neus < np.percentile(t2m_neus, 10))[0]
    # separate into years and months
    std_o3_neus = np.reshape(std_o3_neus, (11, 92))
    std_nox_neus = np.reshape(std_nox_neus, (11, 92))
    emfix_o3_neus = np.reshape(emfix_o3_neus, (11, 92))
    emfix_nox_neus = np.reshape(emfix_nox_neus, (11, 92))
    emiss_o3_neus = np.reshape(emiss_o3_neus, (3, 92))
    emiss_nox_neus = np.reshape(emiss_nox_neus, (3, 92)) 
    mr2_o3_neus = np.reshape(mr2_o3_neus, (3, 92))
    mr2_nox_neus = np.reshape(mr2_nox_neus, (3, 92))
    # for emissions inventories, n.b. repeat monthly mean values for every day and 
    # save off 2000 values for EmFix
    std_inventory_daily, emfix_inventory_daily = [], []  
    emfix_inventory_mm = []  
    for year in std_no_inventory_c_neus: 
        std_inventory_daily.append(np.repeat(year[0], 30))
        std_inventory_daily.append(np.repeat(year[1], 31))
        std_inventory_daily.append(np.repeat(year[2], 31))
        emfix_inventory_daily.append(np.repeat(std_no_inventory_c_neus[0][0], 30))
        emfix_inventory_daily.append(np.repeat(std_no_inventory_c_neus[0][1], 31))
        emfix_inventory_daily.append(np.repeat(std_no_inventory_c_neus[0][2], 31))
        emfix_inventory_mm.append(std_no_inventory_c_neus[0][0])
        emfix_inventory_mm.append(std_no_inventory_c_neus[0][1])
        emfix_inventory_mm.append(std_no_inventory_c_neus[0][2])
    # # # # + Emissions/+ Chemistry simulations on hot and cold days
    # for inventories    
    emiss_inventory_neus_hot = emiss_no_inventory_neus[t2m_hot]
    emiss_inventory_neus_cold = emiss_no_inventory_neus[t2m_cold]
    delta_emiss_inventory = emiss_inventory_neus_hot - emiss_inventory_neus_cold
    mr2_inventory_neus_hot = mr2_no_inventory_neus[t2m_hot]
    mr2_inventory_neus_cold = mr2_no_inventory_neus[t2m_cold]
    delta_mr2_inventory = mr2_inventory_neus_hot - mr2_inventory_neus_cold
    # for NOx
    emiss_nox_hot = np.hstack(emiss_nox_neus)[t2m_hot]
    emiss_nox_cold = np.hstack(emiss_nox_neus)[t2m_cold]
    delta_emiss_nox = emiss_nox_hot - emiss_nox_cold
    mr2_nox_hot = np.hstack(mr2_nox_neus)[t2m_hot]
    mr2_nox_cold = np.hstack(mr2_nox_neus)[t2m_cold]
    delta_mr2_nox = mr2_nox_hot - mr2_nox_cold
    # for O3
    emiss_o3_hot = np.hstack(emiss_o3_neus)[t2m_hot]
    emiss_o3_cold = np.hstack(emiss_o3_neus)[t2m_cold]
    delta_emiss_o3 = emiss_o3_hot - emiss_o3_cold
    mr2_o3_hot = np.hstack(mr2_o3_neus)[t2m_hot]
    mr2_o3_cold = np.hstack(mr2_o3_neus)[t2m_cold]
    delta_mr2_o3 = mr2_o3_hot - mr2_o3_cold
    # stack arrays
    mr2_inventory_daily = np.hstack(mr2_inventory_daily)
    std_inventory_daily = np.hstack(std_inventory_daily)
    emfix_inventory_daily = np.hstack(emfix_inventory_daily)
    emiss_no_inventory_neus = np.hstack(emiss_no_inventory_neus)
    # initialize figure, axes
    fig = plt.figure()
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    # daily values (i.e., EmFix - Std and (+ Emissions) - (+ Chemistry))
    # of modeled NOx versus inventory
    ax1.plot((emfix_inventory_daily - std_inventory_daily), 
             (np.hstack(emfix_nox_neus) - np.hstack(std_nox_neus)), 'o', 
             color = '#666666', markersize = 3, alpha = 1.)
    ax1.plot((delta_emiss_inventory - delta_mr2_inventory), 
             (delta_emiss_nox - delta_mr2_nox), 'o', 
             color = '#d95f02', markersize = 3, alpha = 1.)
    ax1.set_xlabel('$\mathregular{\Delta}$ Emiss$_{\mathregular{NO}}$ '+
                   '[kg s$^{\mathregular{-1}}$ grid cell$^{\mathregular{-1}}$]', 
                   fontsize = 16)         
    ax1.set_ylabel(r'$\mathregular{\Delta}$ NO$_{x}$ [ppbv]', fontsize = 16)
    for t in ax1.get_xticklabels():
        t.set_fontsize(12)
    for t in ax1.get_yticklabels():
        t.set_fontsize(12)  
    ax1.set_title('(a)', fontsize = 16, x = 0.03, y = 1.03, ha = 'left')
    ax1.tick_params(right = True, left = True, top = True, 
                     labelright = False, labelleft = True)    
    # daily values of modeled O3 versus modeled NOx
    ax2.plot((np.hstack(emfix_nox_neus) - np.hstack(std_nox_neus)), 
             (np.hstack(emfix_o3_neus) - np.hstack(std_o3_neus)), 'o', 
             markersize = 3, color = '#666666', alpha = 1., 
             label = 'Strode et al. (2015)')
    ax2.plot((delta_emiss_nox - delta_mr2_nox), (delta_emiss_o3 - delta_mr2_o3), 'o', 
             color = '#d95f02', markersize = 3, alpha = 1., 
             label = 'this study')
    ax2.set_xlabel(r'$\mathregular{\Delta}$ NO$_{x}$ [ppbv]', fontsize = 16)
    ax2.get_xaxis().set_label_coords(0.5, -0.14)    
    ax2.set_ylabel(r'$\mathregular{\Delta}$ O$_{\mathregular{3}}$ [ppbv]', 
                   fontsize = 16)
    for t in ax2.get_xticklabels():
        t.set_fontsize(12)
    for t in ax2.get_yticklabels():
        t.set_fontsize(12)        
    ax2.set_title('(b)', fontsize = 16, x = 0.03, y = 1.03, ha = 'left')
    ax2.tick_params(right = True, left = True, top = True, 
                     labelright = False, labelleft = True)    
    plt.subplots_adjust(wspace = 0.5)
    # add legend 
    leg = ax2.legend(loc = 8, bbox_to_anchor = (-0.4, -0.52), ncol = 2, 
                    fontsize = 16, numpoints = 1)
    leg.get_frame().set_linewidth(0.0)
    plt.subplots_adjust(bottom = 0.3)    
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_inventorynonox_noxo3.eps', dpi = 300)
    # calculate regression between fields 
    print('m(Strode Emiss, Strode NOx) = %.4f ppbv s grid cell kg-1' %
          np.polyfit((emfix_inventory_daily - std_inventory_daily), 
          (np.hstack(emfix_nox_neus) - np.hstack(std_nox_neus)), 1)[0])
    print('m(Emiss, NOx) = %.4f ppbv s grid cell kg-1' %
          np.polyfit((delta_emiss_inventory - delta_mr2_inventory), 
          (delta_emiss_nox - delta_mr2_nox), 1)[0])
    print('m(Strode NOx, Strode O3) = %.4f ppbv ppbv-1' %
          np.polyfit((np.hstack(emfix_nox_neus) - np.hstack(std_nox_neus)), 
          (np.hstack(emfix_o3_neus) - np.hstack(std_o3_neus)), 1)[0])
    print('m(NOx, O3) = %.4f ppbv ppbv-1' %
          np.polyfit((delta_emiss_nox - delta_mr2_nox), 
          (delta_emiss_o3 - delta_mr2_o3), 1)[0])
    return std_inventory_daily, np.hstack(std_o3_neus)
# # # # # # # # # # # # #
def boxplot_cemsnox_castneto3_neus(neus_castnet, std_inventory_daily, 
    std_o3_neus):
    """plots boxplots of regionally-summed daily NOx emissions from CEMS in 
    the Northeastern United States for JJA 2000 - 2013 and boxplots of 
    regionally-averaged daily O3 from 1300-1400 hours (local time) for the 
    same period and region. The summertime mean emission inventory NO and O3
    from the Std. simulation in Strode et al. (2015) for 2000-2010 is plotted.
    
    Parameters
    ----------
    neus_castnet : list
        Site IDs for CASTNet sites in Northeastern U.S., [no. sites in region,]
    std_inventory_daily : numpy.ndarray   
        Daily Northeast-averaged NO emissions from Strode et al. [2015] 
        emissions inventory Std simulation for the measuring period 2000 - 
        2010, units of kg grid cell-1 s-1, [time,]
    std_o3_neus : numpy.ndarray
        Daily Northeast-averaged O3 from Strode et al. [2015] 
        emissions inventory Std simulation for the measuring period 2000 - 
        2010, units of ppbv, [time,]        
        
    Returns
    ----------
    None
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('/Users/ghkerr/phd/emissions/')
    import AQSCEMSobs      
    sampling_months = [6, 7, 8]
    years = [2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 
             2011, 2012, 2013]
    states_ab = ['CT', 'DC', 'DE', 'MA', 'MD', 'ME', 'NH', 'NJ', 'NY', 'PA', 
                 'RI', 'VA', 'VT', 'WV']
    afternoon_times = [13, 14]
    # extract CASTNet observations from 2000-2012
    castnet = find_conus_castnet(years)
    castnet['DATE_TIME'] = pd.to_datetime(castnet['DATE_TIME'])
    # load CEMS NOx emissions in NEUS
    nox_state, nox_lat, nox_lon = AQSCEMSobs.cems_specifystates_dailymean(
            '/Volumes/GAIGEKERR/emissions/CEMS/', states_ab, sampling_months)
    nox = nox_state['%d-0%d-01'%(years[0], sampling_months[0]):
                    '%d-0%d-31'%(years[-1], sampling_months[-1])].values
    # define date range
    date_idx = []
    for year in np.arange(years[0], years[-1]+1, 1):
        date_idx.append(pd.date_range('06-01-%d' %year, '08-31-%d' %year))
    # aggregate over years
    date_idx = np.hstack(date_idx)
    # lists to be filled with data from CASTNet sites
    lat_all, lon_all, o3_all, sites_all = [], [], [], []
    for siteid in np.unique(castnet['SITE_ID'].values):
        castnet_atsite = castnet.loc[castnet['SITE_ID'] == siteid]
        # only proceed if there are observations at site 
        if castnet_atsite['OZONE'].isnull().all() != True:
            lat_atsite = np.nanmean(castnet_atsite['Latitude'].values)
            lon_atsite = np.nanmean(castnet_atsite['Longitude'].values)
            # find CASTNet observations in afternoon 
            castnet_atsite = castnet_atsite.loc[
                    castnet_atsite['DATE_TIME'].dt.hour.isin(afternoon_times)]
            # sort by time
            castnet_atsite = castnet_atsite.sort_values(by = 'DATE_TIME')    
            # find daily average
            castnet_atsite = castnet_atsite.groupby(castnet_atsite['DATE_TIME'].dt.date).mean()
            # add missing values             
            if castnet_atsite.shape[0] != 276:
                castnet_atsite = castnet_atsite.reindex(date_idx, 
                    fill_value = np.nan)
            # add daily values at each sites to lists
            o3_all.append(castnet_atsite['OZONE'].values)
            sites_all.append(siteid)
            lat_all.append(lat_atsite)
            lon_all.append(lon_atsite)
    # reshape results from Strode et al. (2015) to a single value for each 
    # summer
    strode_inventory_sm = np.reshape(strode_inventory, (11, 92))
    strode_inventory_sm = np.mean(strode_inventory_sm, axis = 1)
    strode_o3_sm = np.reshape(strode_o3, (11, 92))
    strode_o3_sm = np.mean(strode_o3_sm, axis = 1)
    # find sites in Northeast
    where_neus = np.in1d(sites_all, neus_castnet)
    where_neus = np.where(where_neus == True)[0]
    # O3 at CASTNet sites in region
    o3_neus = np.array(o3_all)[where_neus]
    # find daily regional average
    o3_neus = np.nanmean(o3_neus, axis = 0)
    # reshape such that first dimension is the year and the second dimension
    # is the number of days in JJA (92)
    o3_neus = np.reshape(o3_neus, (len(years), 92))
    nox = np.reshape(nox, (len(years), 92))
    # initialize figure, axes
    fig = plt.figure(figsize = (11, 4))
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    for i in np.arange(nox.shape[0]):  
        nox_box = ax1.boxplot(nox[i], positions = [i], patch_artist = True, 
                              whis = 1.5, vert = True, sym = '', widths = 0.6)
        ozone_box = ax2.boxplot(o3_neus[i], positions = [i], patch_artist = True, 
                                whis = 1.5, vert = True, sym = '', widths = 0.6)
        # for NOx boxplots
        plt.setp(nox_box['boxes'], color = '#d95f02', facecolor = '#d95f02')
        plt.setp(nox_box['whiskers'], linestyle = '-', color = '#d95f02', 
                 linewidth = 2)
        plt.setp(nox_box['medians'], color = 'w')           
        plt.setp(nox_box['caps'], color = 'k', alpha = 0.0)  
        # for O3 boxplots
        plt.setp(ozone_box['boxes'], color = '#d95f02', facecolor = '#d95f02')
        plt.setp(ozone_box['whiskers'], linestyle = '-', color = '#d95f02',
                 linewidth = 2)                 
        plt.setp(ozone_box['medians'], color = 'w')           
        plt.setp(ozone_box['caps'], color = 'k', alpha = 0.0)          
    # titles
    ax1.set_title('(a)', fontsize = 16, x = 0.1, y = 1.03)
    ax2.set_title('(b)', fontsize = 16, x = 0.1, y = 1.03)
    # x-axis
    ax1.set_xticks(np.arange(-1, 14.5))
    ax2.set_xticks(np.arange(-1, 14.5))
    ax1.set_xticklabels(['', '2000', '', '2002', '', '2004', '', '2006', '', 
                         '2008', '', '2010', '', '2012', ''], fontsize = 12)
    ax2.set_xticklabels(['', '2000', '', '2002', '', '2004', '', '2006', '', 
                         '2008', '', '2010', '', '2012', ''], fontsize = 12)
    # y-axis    
    ax1.set_ylabel('CEMS '+
                   '[tons day$^{\mathregular{-1}}$]', fontsize=16, 
                   color='#d95f02')
    ax1.get_yaxis().set_label_coords(-0.15, 0.53)    
    ax2.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize=16, color='#d95f02')
    ax2.get_yaxis().set_label_coords(-0.11, 0.50)               
    ax1.tick_params(right = True, left = True, top = True, 
                    labelright = False, labelleft = True)
    for t in ax1.get_yticklabels():
        t.set_fontsize(12)    
        t.set_color('#d95f02')
    ax2.tick_params(right = True, left = True, top = True)
    for t in ax2.get_yticklabels():
        t.set_fontsize(12) 
        t.set_color('#d95f02')       
    # add results from Strode et al. (2015) for NO emission inventory
    ax1b = ax1.twinx()
    ax1b.plot(np.arange(0, 11, 1), strode_inventory_sm, 'o-', lw=2.,  
              markersize=6, color='#666666')
    for t in ax1b.get_yticklabels():
        t.set_fontsize(12) 
        t.set_color('#666666')          
    ax1b.set_ylim([2.5, 5.0])
    ax1b.set_ylabel('CTM '+
                    '[kg s$^{\mathregular{-1}}$ grid cell'+
                    ' s$^{\mathregular{-1}}$]', color = '#666666', fontsize = 16, 
                    rotation = 270)
    ax1b.get_yaxis().set_label_coords(1.21, 0.50)
    # add results from Strode et al. (2015) for O3
    ax2b = ax2.twinx()
    ax2b.plot(np.arange(0, 11, 1), strode_o3_sm, 'o-', lw=2., markersize=6, 
              color='#666666')
    ax2b.set_ylim([45, 65])
    ax2b.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize=16, color='#666666',
                    rotation=270)
    for t in ax2b.get_yticklabels():
        t.set_fontsize(12) 
        t.set_color('#666666')    
    ax2b.get_yaxis().set_label_coords(1.24, 0.53)
    plt.subplots_adjust(bottom=0.3, wspace=0.45)
    # create fake boxplot to represent this study's data
    fdat = np.random.normal(loc=0., scale=0.5, size=500)
    # add axis for legend and customize
    axleg = fig.add_axes([0.33, 0.12, 0.3, 0.1])
    # remove spines
    for side in ['top', 'bottom', 'right', 'left']:
        axleg.spines[side].set_visible(False)
    leg = axleg.boxplot([fdat], positions = [0], patch_artist=True, whis=1.5,
                        vert=False, sym='', widths=[0.45])
    axleg.plot([3], [0], 'o', lw=2., markersize=6, color='#666666')
    axleg.plot([2.5, 3.5], [0, 0], '-', lw=2., markersize=6, color='#666666')
    axleg.set_xticks([0, 3])
    axleg.set_xticklabels(['Observations', 'Strode et al. (2015)'], 
                           fontsize = 16)
    axleg.tick_params(axis='x', colors='w', labelcolor='k')
    axleg.set_yticks([])
    axleg.set_yticklabels([])
    plt.setp(leg['boxes'], color = '#d95f02', facecolor = '#d95f02')
    plt.setp(leg['whiskers'], linestyle = '-', color = '#d95f02',
             linewidth = 2)                 
    plt.setp(leg['medians'], color = 'w')           
    plt.setp(leg['caps'], color = 'k', alpha = 0.0)  
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'boxplot_cemsnox_castneto3_neus.eps', dpi = 300)
    return
# # # # # # # # # # # # #    
def map_allgmio3_p90p10(dat_o3, mr2_o3, emiss_o3, emiss_r, gmi_lat, gmi_lon,
    neus): 
    """plot maps of the  difference between the 10th and 90th percentiles from 
    the Transport and + Chemistry simulations and median O3 concentrations from
    the + Emissions simulation. A third map shows the ratio between these 
    event enhancements. Regions where the O3-T correlation (from + Emissions 
    simulation) falls below r = 0.5 are outlined. 

    Parameters
    ----------       
    dat_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, for model case 
        HindcastMR2-DiurnalAvgT, units of volume mixing ratio, [time, lat, lon] 
    mr2_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, for model case 
        HindcastMR2, units of volume mixing ratio, [time, lat, lon] 
    emiss_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, for model case 
        GHKerr-DailyEmiss, units of volume mixing ratio, [time, lat, lon]         
    emiss_r : numpy.ndarray
        The Pearson product-moment correlation coefficient calculated between 
        MERRA-2 2-meter temperatures and ozone at each GMI grid cell from the 
        + Emissions simulation, [lat, lon]
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
    neus : numpy.ndarray
        CTM grid where grid cells within NEUS have a value of 1 and grid 
        cells not in region have a value of NaN, [lat, lon]             
        
    Returns
    ----------      
    None                
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from matplotlib.colors import Normalize
    from matplotlib.colorbar import ColorbarBase    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # extreme O3 events in each simulation 
    dat_p90 = np.percentile(dat_o3, 90, axis = 0) * 1e9
    mr2_p90 = np.percentile(mr2_o3, 90, axis = 0) * 1e9
    dat_p10 = np.percentile(dat_o3, 10, axis = 0) * 1e9
    mr2_p10 = np.percentile(mr2_o3, 10, axis = 0) * 1e9    
    med_o3 = np.percentile(emiss_o3, 50, axis = 0)  * 1e9
    m = Basemap(projection = 'merc', llcrnrlon = -126., llcrnrlat = 24., 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'h', 
                area_thresh = 1000)
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)    
    # # # # for extreme high events
    vmin = 5.; vmax = 15.
    cmap = plt.get_cmap('PuBu', 10)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    # contributions from Transport simulation
    fig = plt.figure(figsize = (4, 7))    
    ax1 = plt.subplot2grid((3, 1), (0, 0))
    ax1.set_title('(a) Transport', ha = 'left', fontsize = 16, x = 0.03, 
                  y = 1.03)
    m.contourf(x, y, dat_p90 - med_o3, clevs, cmap = cmap, norm = norm, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax1, m)    
    # + Chemistry simulation
    ax2 = plt.subplot2grid((3, 1), (1, 0))
    ax2.set_title('(b) $\mathregular{+}$Chemistry', ha = 'left', fontsize = 16, x = 0.03, 
                  y = 1.03)
    m.contourf(x, y, mr2_p90 - med_o3, clevs, cmap = cmap, norm = norm, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax2, m)    
    plt.subplots_adjust(right = 0.7)
    cbar_ax = fig.add_axes([0.75, 0.43, 0.05, 0.4])
    cb = ColorbarBase(cbar_ax, cmap = cmap, norm = norm, extend = 'both',
                      orientation = 'vertical')
    cb.set_ticks(np.linspace(vmin, vmax, 6))
    cb.set_label(label = 'O$_{\mathregular{3, P}_{\mathregular{90}}}$ - '+
                 'O$_{\mathregular{3, P}_{\mathregular{50}}}$ [ppbv]', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)    
    # ratio from simulations
    ax3 = plt.subplot2grid((3, 1), (2, 0))
    ax3.set_title('(c) Transport/$\mathregular{+}$Chemistry', ha = 'left', fontsize = 16, 
                  x = 0.03, y = 1.03)
    vmin = 75.; vmax = 125.
    cmap = plt.get_cmap('RdBu_r', 10)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m.contourf(x, y, ((dat_p90 - med_o3)/(mr2_p90 - med_o3)) * 100., clevs, 
               cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax3, m)    
    m.contour(x, y, emiss_r, levels = [0.4], colors = ('orange',),
              linestyles = ('--',), linewidths = (1,), zorder = 1)     
    cbar_ax = fig.add_axes([0.75, 0.16, 0.05, 0.16])
    cb = ColorbarBase(cbar_ax, cmap = cmap, norm = norm, extend = 'both',
                      orientation = 'vertical')
    cb.set_ticks(np.linspace(vmin, vmax, 6))
    cb.set_label(label = '[%]', size = 16)
    cb.ax.tick_params(labelsize = 12)    
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_allgmio3_p90.eps', dpi = 300)
    plt.show()
    # # # # for extreme low events
    vmin = -15.; vmax = -5.
    cmap = plt.get_cmap('PuBu_r', 10)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    # contributions from Transport simulation
    fig = plt.figure(figsize = (4, 7))    
    ax1 = plt.subplot2grid((3, 1), (0, 0))
    ax1.set_title('(a) Transport', ha = 'left', fontsize = 16, x = 0.03, 
                  y = 1.03)
    m.contourf(x, y, (dat_p10 - med_o3), clevs, cmap = cmap, norm = norm, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax1, m)    
    # + Chemistry simulation
    ax2 = plt.subplot2grid((3, 1), (1, 0))
    ax2.set_title('(b) $\mathregular{+}$Chemistry', ha = 'left', fontsize = 16, x = 0.03, 
                  y = 1.03)
    m.contourf(x, y, (mr2_p10 - med_o3), clevs, cmap = cmap, norm = norm, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax2, m)    
    plt.subplots_adjust(right = 0.7)
    cbar_ax = fig.add_axes([0.75, 0.43, 0.05, 0.4])
    cb = ColorbarBase(cbar_ax, cmap = cmap, norm = norm, extend = 'both',
                      orientation = 'vertical')
    cb.set_ticks(np.linspace(vmin, vmax, 6))
    cb.set_label(label = 'O$_{\mathregular{3, P}_{\mathregular{10}}}$ - '+
                 'O$_{\mathregular{3, P}_{\mathregular{50}}}$ [ppbv]', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)    
    # ratio from simulations
    ax3 = plt.subplot2grid((3, 1), (2, 0))
    ax3.set_title('(c) Transport/$\mathregular{+}$Chemistry', ha = 'left', fontsize = 16, 
                  x = 0.03, y = 1.03)
    vmin = 75.; vmax = 125.
    cmap = plt.get_cmap('RdBu_r', 10)
    norm = Normalize(vmin = vmin, vmax = vmax)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    m.contourf(x, y, ((dat_p10 - med_o3)/(mr2_p10 - med_o3)) * 100., clevs, 
               cmap = cmap, extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.0)
    m.drawcoastlines(color = 'k', linewidth = 1.0)    
    fill_oceans(ax3, m)    
    m.contour(x, y, emiss_r, levels = [0.4], colors = ('orange',),
              linestyles = ('--',), linewidths = (1,), zorder = 1)     
    cbar_ax = fig.add_axes([0.75, 0.16, 0.05, 0.16])
    cb = ColorbarBase(cbar_ax, cmap = cmap, norm = norm, extend = 'both',
                      orientation = 'vertical')
    cb.set_ticks(np.linspace(vmin, vmax, 6))
    cb.set_label(label = '[%]', size = 16)
    cb.ax.tick_params(labelsize = 12)    
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_allgmio3_p10.eps', dpi = 300)
    return
# # # # # # # # # # # # #
def timeseries_map_hourlyvsoverpass_neus(o3, castnet_o3_neus, emiss_o3_neus,
    gmi_lat, gmi_lon, year, years):
    """for the desired summer, function plots a timeseries of regionally-
    averaged O3 in the Northeast from (1) CASTNet observations, (2) daily 
    average over 1300-1400 hours (local time) from hourly output at GMI sites,
    and (3) the Northeastern-averaged from output at overpass2 time from 
    gridded GMI. The bottom subplot shows the location of GMI sites with hourly 
    output and the mean difference between hourly output and co-located 
    overpass2 output. 

    Parameters
    ----------       
    o3 : numpy.ndarray
        GMI CTM daily O3 concentrations for model case GHKerr-DailyEmiss, units 
        of ppbv, [lat, lon, time,]
    castnet_o3_neus : numpy.ndarray
        Daily regionally-averaged O3 from CASTNet, units of ppbv, [time,]         
    emiss_o3_neus : numpy.ndarray
        Daily regionally-averaged O3 from GHKerr-DailyEmiss simulation, units 
        of ppbv, [time,]             
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]          
    year : int
        Year of interest
    years : list
        Years in measuring period
        
    Returns
    ----------      
    None             
    """
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.colorbar import ColorbarBase
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx 
    import commensurability
    years = [2008, 2009, 2010]            
    sampling_months = [6, 7, 8]
    sampling_hours = [15, 16, 17, 18, 19, 20]
    llcrnrlon = -83.
    llcrnrlat = 37.0
    urcrnrlon = -66.3
    urcrnrlat = 48.
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
        llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat, 
        resolution = 'h', area_thresh = 1000)
    castnet_sites_neus = ['ASH', 'WFM', 'WST', 'APT', 'SAR', 'CTH', 'WSP',
                        'ARE', 'BEL', 'PSU', 'LRL', 'PAR', 'PED', 'SHN', 
                        'CDR', 'VPI', 'MKG', 'KEF']
    timezone = 'GMT+4'
    # CASTNet sites, n.b. MR2 CASTNet data (with no regional average) needs 
    # to be loaded
    mr2_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
    commensurability.commensurate_castnet_gmi(castnet_sites_neus, 'HindcastMR2', 
        years, sampling_months, sampling_hours, timezone)
    # drop GMI hourly output that's missing
    mask = np.all(np.isnan(mr2_o3), axis=tuple((0,2)))
    mr2_o3 = mr2_o3[:,~mask]
    # # # # open daily mean CASTNet O3 and GMI trace gases
    meandaily = commensurability.commensurate_castnet_gmi_ra(castnet_sites_neus, 
        years, sampling_months, sampling_hours, timezone)
    castnet_lats, castnet_lons, gmi_lats, gmi_lons = \
        commensurability.commensurate_castnet_gmi_locations(castnet_sites_neus, 
        sampling_months, sampling_hours, timezone)
    whereyear = np.where(np.array(years)==year)[0][0]
    meano3_overpass,meano3_hourly = [],[]
    # Initialize figure, axes 
    fig = plt.figure()
    axt = plt.subplot2grid((2,1),(0,0))
    axt.set_title('(a)', fontsize = 16, x = -0.07, y = 1.03)
    axt.plot(castnet_o3_neus[whereyear*92:(whereyear+1)*92], '-k', lw=2, 
                             label='CASTNet')
    axt.plot(emiss_o3_neus[whereyear*92:(whereyear+1)*92], '-', color='#d95f02',
                         lw=2, label='Overpass')
    axt.plot(meandaily['EGU_T O3'][whereyear], '-', color='#1b9e77', lw=2, 
             label='Hourly')
    axt.set_xlim([0, 91])
    axt.set_ylim([30, 72])
    # axis labels
    axt.set_xticks([0, 14, 30, 44, 61, 75])
    axt.set_xticklabels(['1 June %d' %year, '', '1 July', '', '1 Aug', ''], 
                       ha = 'center', fontsize = 12)
    axt.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    plt.legend(loc=1, frameon = False, fontsize = 16, ncol=3,
               bbox_to_anchor = (1.15, 1.4))
    # Map subplot
    axb = plt.subplot2grid((2,1),(1,0))
    axb.set_ylabel('(b)', fontsize = 16, rotation=0)
    axb.yaxis.set_label_coords(-0.22,0.75)
    vmin = -5; vmax = 5
    cmap = plt.get_cmap('bwr', 12)
    clevs = np.linspace(vmin, vmax, 11, endpoint = True)
    x_gmi, y_gmi = m(gmi_lons, gmi_lats)
    # find overpass2 concentrations at each hourly profile site 
    i = 0
    for ilat, ilon in zip(gmi_lats, gmi_lons):
        gilat = geo_idx(ilat,gmi_lat)
        gilon = geo_idx(ilon,gmi_lon)
        io3 = o3[:,gilat,gilon]
        meano3_overpass.append(np.mean(io3))
        meano3_hourly.append(np.mean(mr2_o3[:,i]))
        i = i + 1
    diff = np.array(meano3_overpass) - np.array(meano3_hourly)
    sc = m.scatter(x_gmi, y_gmi, c = diff, s = 20, cmap = cmap, vmin = vmin, 
                   vmax = vmax, zorder = 30, linewidth = 1., edgecolor = 'k')    
    fill_oceans(axb, m)    
    # add countries, states, etc here so they sit atop the gridlines
    m.drawstates(color = 'k', linewidth = 0.5, zorder = 10)
    m.drawcountries(color = 'k', linewidth = 1.0, zorder = 10)
    m.drawcoastlines(color = 'k', linewidth = 1.0, zorder = 10)
    cax = fig.add_axes([0.30, 0.15, 0.43, 0.05])
    cb = fig.colorbar(sc, cax=cax, orientation='horizontal', extend = 'both')
    cb.set_ticks(np.linspace(vmin, vmax, 5, endpoint = True))
    cb.set_label(label = '$\mathregular{O}_\mathregular{3,\:Overpass}$ '+
        '$-$ ${\mathregular{O}_\mathregular{3,\:Hourly}}$ [ppbv]', size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(hspace=0.33, bottom=0.22)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/'+
                'timeseries_map_hourlyvsoverpass_neus.eps', dpi=300)
    return    
# # # # # # # # # # # # #  
def scatter_noxo3(mr2_o3, emiss_o3, mr2_no, emiss_no, mr2_no2, emiss_no2, 
    neus_states, gmi_lat, gmi_lon):
    """function finds the difference between modeled NOx and O3 in the 
    Northeastern U.S. for (1) Std and EmFix simulations from Strode et al. 
    (2015) and (2) the + Chemistry and + Emissions simulations in current work. 
    These daily differences are plotted as scatterpoints and the lines of best 
    fit are calculated. 
    
    Parameters
    ----------  
    mr2_o3 : numpy.ndarray    
        GMI CTM surface-level ozone at overpass time from the + Chemistry 
        simulation, units of volume mixing ratio, [time, lat, lon]
    emiss_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from the + Emissions 
        simulation, units of volume mixing ratio, [time, lat, lon]
    mr2_no : numpy.ndarray    
        GMI CTM surface-level NO at overpass time from the + Chemistry 
        simulation, units of volume mixing ratio, [time, lat, lon]
    emiss_no : numpy.ndarray
        GMI CTM surface-level NO at overpass time from the + Emissions 
        simulation, units of volume mixing ratio, [time, lat, lon]        
    mr2_no2 : numpy.ndarray    
        GMI CTM surface-level NO2 at overpass time from the + Chemistry 
        simulation, units of volume mixing ratio, [time, lat, lon]
    emiss_no2 : numpy.ndarray
        GMI CTM surface-level NO2 at overpass time from the + Emissions 
        simulation, units of volume mixing ratio, [time, lat, lon]        
    neus_states : list
        State names (string format) in region      
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,] 

    Returns
    ----------          
    None
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/GMI/')
    import commensurability 
    # load Std
    (gmi_lat_c, gmi_lon_c, eta_c, times_c, std_co, std_no, std_no2, std_o3, 
     std_cloudfraction, std_gridboxheight) = \
    commensurability.open_overpass2('HindcastFFIgac2', [2000, 2001, 2002, 
        2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010])
    # load EmFix
    (gmi_lat_c, gmi_lon_c, eta_c, times_c, emfix_co, emfix_no, emfix_no2, emfix_o3, 
     std_cloudfraction, std_gridboxheight) = \
    commensurability.open_overpass2('Hindcast3Igac2', [2000, 2001, 2002, 
        2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010])
    gmi_lon_c = np.mod(gmi_lon_c - 180.0, 360.0) - 180.0    
    # emissions inventory from + Emissions simulation
    emiss_no_inventory, lon_inventory, lat_inventory = \
     commensurability.open_perturbed_emissions()
    # emissions inventory from + Chemistry simulation 
    mr2_no_inventory, lon_inventory, lat_inventory = \
     commensurability.open_unperturbed_emissions('1x1.25_IAVanthGFED4')
    lon_inventory = np.mod(lon_inventory[0] - 180.0, 360.0) - 180.0
    # emissions inventory from Std/EmFix simulations
    std_no_inventory, lon_inventory_c, lat_inventory_c = \
     commensurability.open_unperturbed_emissions('2x2.5_IAVanthGFED3gcEF')
    lon_inventory_c = np.mod(lon_inventory_c[0] - 180.0, 360.0) - 180.0
    # find grid cells in Northeast for different simulations/inventories
    m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 24.0, 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'c', 
                area_thresh = 1000)
    neus_states = pollutants_constants.NORTHEAST_STATES
    neus = find_grid_in_region(m, neus_states, gmi_lat, gmi_lon)
    neus_c = find_grid_in_region(m, neus_states, gmi_lat_c, gmi_lon_c) 
    neus_inventory = find_grid_in_region(m, neus_states, lat_inventory[0], 
                                         lon_inventory)  
    neus_inventory_c = find_grid_in_region(m, neus_states, lat_inventory_c[0], 
                                           lon_inventory_c)       
    # # # # find Northeast-averaged inventories, NOx, O3 
    # from Strode et al. [2015] simulations
    std_no_inventory_c_neus = np.array(std_no_inventory) * neus_inventory_c
    std_no_inventory_c_neus = np.nanmean(std_no_inventory_c_neus, axis = tuple((2, 3)))
    emiss_no_inventory_neus = np.vstack(emiss_no_inventory) * neus_inventory
    emiss_no_inventory_neus = np.nanmean(emiss_no_inventory_neus, axis = tuple((1, 2)))             
    mr2_no_inventory_neus = np.vstack(mr2_no_inventory) * neus_inventory
    mr2_no_inventory_neus = np.nanmean(mr2_no_inventory_neus, axis = tuple((1, 2)))    
    # since the inventory from MR2 only varies from month-to-month, repeat 
    # monthly values for every day in month
    mr2_inventory_daily = []
    for year in (mr2_no_inventory_neus[0:3], 
                 mr2_no_inventory_neus[3:6],
                 mr2_no_inventory_neus[6:]): 
        mr2_inventory_daily.append(np.repeat(year[0], 30))
        mr2_inventory_daily.append(np.repeat(year[1], 31))
        mr2_inventory_daily.append(np.repeat(year[2], 31))
    mr2_no_inventory_neus = np.hstack(mr2_inventory_daily)
    emiss_o3_neus = np.nanmean((emiss_o3 * neus), axis = tuple((1, 2))) * 1e9
    mr2_o3_neus = np.nanmean((mr2_o3 * neus), axis = tuple((1, 2))) * 1e9
    # for emissions inventories, n.b. repeat monthly mean values for every day and 
    # save off 2000 values for EmFix
    std_inventory_daily, emfix_inventory_daily = [], []  
    emfix_inventory_mm = []  
    for year in std_no_inventory_c_neus: 
        std_inventory_daily.append(np.repeat(year[0], 30))
        std_inventory_daily.append(np.repeat(year[1], 31))
        std_inventory_daily.append(np.repeat(year[2], 31))
        emfix_inventory_daily.append(np.repeat(std_no_inventory_c_neus[0][0], 30))
        emfix_inventory_daily.append(np.repeat(std_no_inventory_c_neus[0][1], 31))
        emfix_inventory_daily.append(np.repeat(std_no_inventory_c_neus[0][2], 31))
        emfix_inventory_mm.append(std_no_inventory_c_neus[0][0])
        emfix_inventory_mm.append(std_no_inventory_c_neus[0][1])
        emfix_inventory_mm.append(std_no_inventory_c_neus[0][2])    
    # stack arrays
    mr2_inventory_daily = np.hstack(mr2_inventory_daily)
    std_inventory_daily = np.hstack(std_inventory_daily)
    emfix_inventory_daily = np.hstack(emfix_inventory_daily)
    emiss_no_inventory_neus = np.hstack(emiss_no_inventory_neus)    
    # find grid cells in Northeast for different simulations/inventories
    m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 24.0, 
                urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'c', 
                area_thresh = 1000)
    neus_states = pollutants_constants.NORTHEAST_STATES
    neus = find_grid_in_region(m, neus_states, gmi_lat, gmi_lon)
    neus_c = find_grid_in_region(m, neus_states, gmi_lat_c, gmi_lon_c) 
    # # # # find Northeast-averaged NOx, O3 
    # from Strode et al. [2015] simulations
    std_nox_c_neus = np.array(std_no + std_no2) * neus_c
    std_nox_c_neus = np.nanmean(std_nox_c_neus, axis = tuple((1, 2))) * 1e9
    emfix_nox_c_neus = np.array(emfix_no + emfix_no2) * neus_c
    emfix_nox_c_neus = np.nanmean(emfix_nox_c_neus, axis = tuple((1, 2))) * 1e9
    std_o3_neus = np.nanmean((std_o3 * neus_c), axis = tuple((1, 2))) * 1e9
    emfix_o3_neus = np.nanmean((emfix_o3 * neus_c), axis = tuple((1, 2))) * 1e9
    # from + Emissions/+ Chemistry simulations
    emiss_nox_neus = (emiss_no + emiss_no2) * neus
    emiss_nox_neus = np.nanmean(emiss_nox_neus, axis = tuple((1, 2)))* 1e9           
    mr2_nox_neus = (mr2_no + mr2_no2) * neus
    mr2_nox_neus = np.nanmean(mr2_nox_neus, axis = tuple((1, 2))) * 1e9
    emiss_o3_neus = (emiss_o3) * neus
    emiss_o3_neus = np.nanmean(emiss_o3_neus, axis = tuple((1, 2)))* 1e9           
    mr2_o3_neus = (mr2_o3) * neus
    mr2_o3_neus = np.nanmean(mr2_o3_neus, axis = tuple((1, 2))) * 1e9
    kerr_delta_no = (emiss_nox_neus - mr2_nox_neus)
    kerr_delta_o3 = (emiss_o3_neus-mr2_o3_neus)
    strode_delta_no = (emfix_nox_c_neus - std_nox_c_neus)
    strode_delta_o3 = (emfix_o3_neus - std_o3_neus)
    # Plotting
    fig = plt.figure()
    ax1 = plt.subplot2grid((1,1),(0,0))
    ax1.plot(kerr_delta_no, kerr_delta_o3, 'o', 
             color = '#d95f02', markersize = 3, alpha = 1.)
    ax1.plot(strode_delta_no, strode_delta_o3,  'o', 
             color = '#666666', markersize = 3, alpha = 1.)
    # Plot lines of best fit 
    kerrx = np.arange(kerr_delta_no.min()-0.2, kerr_delta_no.max()+0.2, 0.01)
    strodex = np.arange(strode_delta_no[92:].min()-0.25, 
                        strode_delta_no[92:].max(), 0.01) 
    kerrm, kerrb = np.polyfit(kerr_delta_no, kerr_delta_o3, 1)
    strodem, strodeb = np.polyfit(strode_delta_no[92:], 
                                  strode_delta_o3[92:], 1)         
    ax1.plot(kerrx, kerrb + kerrm * kerrx, '--', lw=2., color='#d95f02', 
             zorder =10)
    ax1.plot(strodex, strodeb + strodem * strodex, '--', lw=2., 
             color='#666666')    
    # Indicate slope of best fit lines (i.e., the sensitivity of O3 to changes
    # in NOx emissions) 
    ax1.text(0.05, 0.73, 'This study:\n' + r'$\frac{\mathregular{\Delta\:O}_'
             r'{\mathregular{3}}}{\mathregular{\Delta\:NO}_{x}}$ = %.2f '%kerrm+
             r'$\frac{\mathregular{ppbv}}{\mathregular{ppbv}}$', 
             transform=ax1.transAxes, color = '#d95f02', fontsize=16)
    ax1.text(0.54, 0.06, 'Strode et al. (2015):\n' + r'$\frac{\mathregular'+
             r'{\Delta\:O}_{\mathregular{3}}}{\mathregular{\Delta\:NO}'+
             r'_{x}}$ = %.2f $\frac{\mathregular{ppbv}}'%strodem +
             r'{\mathregular{ppbv}}$',
             transform=ax1.transAxes, color = '#666666', fontsize=16)
    ax1.set_xlabel('$\mathregular{\Delta}$ NO$_x$ '+
                   '[ppbv]', 
                   fontsize = 16)  
    ax1.set_ylabel('$\mathregular{\Delta}$ O$_\mathregular{3}$ [ppbv]',
                   fontsize = 16)         
    for t in ax1.get_xticklabels():
        t.set_fontsize(12)
    for t in ax1.get_yticklabels():
        t.set_fontsize(12)    
    plt.subplots_adjust(bottom=0.15)      
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_noxo3.eps', dpi = 300)        
    return std_inventory_daily, std_o3_neus
# # # # # # # # # # # # #    
def timeseries_map_hourlyvsoverpass_neus_nomap(castnet_o3_neus, emiss_o3_neus,
    year, years):
    """same as function 'timeseries_map_hourlyvsoverpass_neus' but with no 
    map showing the bias at individual GMI sites with hourly output. 
    
    castnet_o3_neus : numpy.ndarray
        Daily regionally-averaged O3 from CASTNet, units of ppbv, [time,]         
    emiss_o3_neus : numpy.ndarray
        Daily regionally-averaged O3 from GHKerr-DailyEmiss simulation, units 
        of ppbv, [time,]             
    year : int
        Year of interest
    years : list
        Years in measuring period
        
    Returns
    ----------      
    None             
    """
    import numpy as np
    from matplotlib.lines import Line2D
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    import commensurability
    years = [2008, 2009, 2010]            
    sampling_months = [6, 7, 8]
    sampling_hours = [15, 16, 17, 18, 19, 20]
    castnet_sites_neus = ['ASH', 'WFM', 'WST', 'APT', 'SAR', 'CTH', 'WSP',
                        'ARE', 'BEL', 'PSU', 'LRL', 'PAR', 'PED', 'SHN', 
                        'CDR', 'VPI', 'MKG', 'KEF']
    timezone = 'GMT+4'
    # CASTNet sites, n.b. MR2 CASTNet data (with no regional average) needs 
    # to be loaded
    mr2_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
    commensurability.commensurate_castnet_gmi(castnet_sites_neus, 'HindcastMR2', 
        years, sampling_months, sampling_hours, timezone)
    # drop GMI hourly output that's missing
    mask = np.all(np.isnan(mr2_o3), axis=tuple((0,2)))
    mr2_o3 = mr2_o3[:,~mask]
    # # # # open daily mean CASTNet O3 and GMI trace gases
    meandaily = commensurability.commensurate_castnet_gmi_ra(castnet_sites_neus, 
        years, sampling_months, sampling_hours, timezone)
    castnet_lats, castnet_lons, gmi_lats, gmi_lons = \
        commensurability.commensurate_castnet_gmi_locations(castnet_sites_neus, 
        sampling_months, sampling_hours, timezone)
    whereyear = np.where(np.array(years)==year)[0][0]
    # Initialize figure, axes 
    fig = plt.figure(figsize=(6,3))
    axt = plt.subplot2grid((1,1),(0,0))
    #axt.set_title('(a)', fontsize = 16, x = -0.07, y = 1.03)
    # Overpass approach
    axt.plot(castnet_o3_neus[whereyear*92:(whereyear+1)*92], '-', 
             color='#d95f02', lw=2)
    axt.plot(emiss_o3_neus[whereyear*92:(whereyear+1)*92], '--', 
             color='#d95f02', lw=0.75)
    # Hourly approach
    axt.plot(np.nanmean(mr2_castnet[whereyear], axis=0), '-', color='#1b9e77',
                         lw=2, zorder = 10)
    axt.plot(meandaily['EGU_T O3'][whereyear], '--', color='#1b9e77', lw=0.75, 
             zorder = 10)
    axt.set_xlim([0, 91])
    axt.set_ylim([30, 70])
    axt.set_yticks([30, 40, 50, 60, 70])
    axt.set_yticklabels(['30', '40', '50', '60', '70'], fontsize=12)
    # axis labels
    axt.set_xticks([0, 14, 30, 44, 61, 75])
    axt.set_xticklabels(['1 June %d' %year, '', '1 July', '', '1 Aug', ''], 
                       ha = 'center', fontsize = 12)
    axt.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    axt.tick_params(right = True, left = True, top = True, bottom = True,
                    labelbottom = True, labeltop = False)
    # create legend
    colors = ['black', 'black', '#1b9e77', '#d95f02']
    sizes = [2, 0.75, 2, 2]
    linestyles = ['-', '--', '-', '-']
    labels = ['CASTNet', 'CTM', 'Hourly', 'Overpass']
    lines = []
    for c, s, l in zip(colors, sizes, linestyles):
        lines.append(Line2D([0], [0], color=c, linewidth=s, linestyle=l))
    plt.legend(lines, labels, loc = 9, frameon = False, fontsize = 16, 
               ncol = 2, bbox_to_anchor = (0.5, 1.48))
    plt.subplots_adjust(top=0.7)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/'+
                'timeseries_map_hourlyvsoverpass_neus_nomap.eps', dpi=300)
    return 
# # # # # # # # # # # # #    
def timeseries_transportchemistryo3_atpoints(t2m_overpass, o3, dat_o3, gmi_lat, 
    gmi_lon):
    """same as function 'timeseries_mr2o3dato3t2m_atpoint' but shows the 
    O3 from the Transport and + Chemistry simulations and co-located 
    temperatures at four grid cells (Montana, Ohio, California, Mississippi)
    
    Parameters
    ----------  
    t2m_atoverpass : numpy.ndarry
        The MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]    
    o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from MR2/+ Chemistry 
        simulation, units of volume mixing ratio, [time, lat, lon]        
    dat_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time from Diurnal-AvgT/
        Transport simulation, units of volume mixing ratio, [time, lat, lon]      
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
        
    Returns
    ----------     
    None         
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # indices of grid cell in South
    slon_idx_south = 33
    slat_idx_south = 8 + 2
    print('(%.3f, %.3f) is chosen coordinate point in South...' 
          %(gmi_lat[slat_idx_south], gmi_lon[slon_idx_south]))   
    # select values at point
    o3_atpoint = o3[184:, slat_idx_south, slon_idx_south] * 1e9
    dato3_atpoint = dat_o3[184:, slat_idx_south, slon_idx_south] * 1e9
    t2m_atpoint = t2m_overpass[184:, slat_idx_south, slon_idx_south]
    # calculate correlation coefficients/sensitivity
    print('r(Transport O3, + Chemistry O3) = %.3f' 
          %(np.corrcoef(o3_atpoint, dato3_atpoint)[0, 1]))
    print('r(+ Chemistry O3, T2m) = %.3f' 
          %(np.corrcoef(o3_atpoint, t2m_atpoint)[0, 1]))
    print('d(+ Chemistry O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_atpoint, o3_atpoint, 1)[0]))
    print('d(Transport O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_atpoint, dato3_atpoint, 1)[0]))
    print('Mean/Std Transport O3 = %.3f, %.3f' 
          %(np.mean(dato3_atpoint), np.std(dato3_atpoint)))
    print('Mean + Chemistry O3 = %.3f, %.3f' 
          %(np.mean(o3_atpoint), np.std(o3_atpoint)))
    slon_idx_north = 39
    slat_idx_north = 18
    print('(%.3f, %.3f) is chosen coordinate point in North...' 
          %(gmi_lat[slat_idx_north], gmi_lon[slon_idx_north]))   
    # select values at point
    o3_atpoint_north = o3[184:, slat_idx_north, slon_idx_north] * 1e9
    dato3_atpoint_north = dat_o3[184:, slat_idx_north, slon_idx_north] * 1e9
    t2m_atpoint_north = t2m_overpass[184:, slat_idx_north, slon_idx_north]
    # calculate correlation coefficients/sensitivity
    print('r(Transport O3, + Chemistry O3) = %.3f' 
          %(np.corrcoef(o3_atpoint_north, dato3_atpoint_north)[0, 1]))
    print('r(+ Chemistry O3, T2m) = %.3f' 
          %(np.corrcoef(o3_atpoint_north, t2m_atpoint_north)[0, 1]))
    print('d(+ Chemistry O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_atpoint_north, o3_atpoint_north, 1)[0]))
    print('d(Transport O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_atpoint_north, dato3_atpoint_north, 1)[0]))
    print('Mean/Std Transport O3 = %.3f, %.3f' 
          %(np.mean(dato3_atpoint_north), np.std(dato3_atpoint_north)))
    print('Mean + Chemistry O3 = %.3f, %.3f' 
          %(np.mean(o3_atpoint_north), np.std(o3_atpoint_north)))
    # indices of grid cell in Montana
    slon_idx_mt = 20
    slat_idx_mt = 22
    print('(%.3f, %.3f) is chosen coordinate point in Montana...' 
          %(gmi_lat[slat_idx_mt], gmi_lon[slon_idx_mt]))   
    # select values at point
    o3_mt = o3[184:, slat_idx_mt, slon_idx_mt] * 1e9
    dato3_mt = dat_o3[184:, slat_idx_mt, slon_idx_mt] * 1e9
    t2m_mt = t2m_overpass[184:, slat_idx_mt, slon_idx_mt]
    # calculate correlation coefficients/sensitivity
    print('r(Transport O3, + Chemistry O3) = %.3f' 
          %(np.corrcoef(o3_mt, dato3_mt)[0, 1]))
    print('r(+ Chemistry O3, T2m) = %.3f' 
          %(np.corrcoef(o3_mt, t2m_mt)[0, 1]))
    print('d(+ Chemistry O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_mt, o3_mt, 1)[0]))
    print('d(Transport O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_mt, dato3_mt, 1)[0]))
    print('Mean/Std Transport O3 = %.3f, %.3f' 
          %(np.mean(dato3_mt), np.std(dato3_mt)))
    print('Mean + Chemistry O3 = %.3f, %.3f' 
          %(np.mean(o3_mt), np.std(o3_mt)))
    # indices of grid cell in California
    slon_idx_ca = 7
    slat_idx_ca = 15
    print('(%.3f, %.3f) is chosen coordinate point in California...' 
          %(gmi_lat[slat_idx_ca], gmi_lon[slon_idx_ca]))   
    # select values at point
    o3_ca = o3[184:, slat_idx_ca, slon_idx_ca] * 1e9
    dato3_ca = dat_o3[184:, slat_idx_ca, slon_idx_ca] * 1e9
    t2m_ca = t2m_overpass[184:, slat_idx_ca, slon_idx_ca]
    # calculate correlation coefficients/sensitivity
    print('r(Transport O3, + Chemistry O3) = %.3f' 
          %(np.corrcoef(o3_ca, dato3_ca)[0, 1]))
    print('r(+ Chemistry O3, T2m) = %.3f' 
          %(np.corrcoef(o3_ca, t2m_ca)[0, 1]))
    print('d(+ Chemistry O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_ca, o3_ca, 1)[0]))
    print('d(Transport O3)/d(T2m) = %.3f' 
          %(np.polyfit(t2m_ca, dato3_ca, 1)[0]))
    print('Mean/Std Transport O3 = %.3f, %.3f' 
          %(np.mean(dato3_ca), np.std(dato3_ca)))
    print('Mean + Chemistry O3 = %.3f, %.3f' 
          %(np.mean(o3_ca), np.std(o3_ca)))
    # Plotting
    fig = plt.figure(figsize=(10,4))
    axlt = plt.subplot2grid((2, 2), (0, 0))
    axlb = plt.subplot2grid((2, 2), (1, 0))
    axrt = plt.subplot2grid((2, 2), (0, 1))
    axrb = plt.subplot2grid((2, 2), (1, 1))
    # Montana
    axlt.set_title('(a) Montana', fontsize = 16, x = 0.2, y = 1.03)
    lns1 = axlt.plot(o3_mt, lw = 2., 
                   c = pollutants_constants.COLOR_CHEMISTRY, 
                   label = '$\mathregular{+}$Chemistry', zorder = 10)
    lns2 = axlt.plot(dato3_mt, lw = 2., 
                   c = pollutants_constants.COLOR_TRANSPORT,
                   label = 'Transport', zorder = 5)
    axlt.text(0.96, 0.05, '%.2f, %.2f, %.2f'
              %(np.corrcoef(o3_mt, t2m_mt)[0, 1],
                np.polyfit(t2m_mt, o3_mt, 1)[0],
                np.polyfit(t2m_mt, dato3_mt, 1)[0]), 
              transform=axlt.transAxes,fontsize=14, ha='right')        
    for t in axlt.get_yticklabels():
        t.set_fontsize(12)   
    axlt.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    axltb = axlt.twinx()
    lns3 = axltb.plot(t2m_mt, lw = 2., color = '#ff7f00',
             label = 'T', zorder = 1)
    axltb.set_xticks([0, 14, 29, 45, 61, 76])
    axltb.set_xticklabels([])
    axltb.set_yticklabels([])
    # California
    axlb.set_title('(c) California', fontsize = 16, x = 0.2, y = 1.03)
    lns1 = axlb.plot(o3_ca, lw = 2., c = pollutants_constants.COLOR_CHEMISTRY, 
                     zorder = 10)
    lns2 = axlb.plot(dato3_ca, lw = 2., 
                     c = pollutants_constants.COLOR_TRANSPORT,
                     zorder = 5)
    axlb.text(0.96, 0.05, '%.2f, %.2f, %.2f'
              %(np.corrcoef(o3_ca, t2m_ca)[0, 1],
                np.polyfit(t2m_ca, o3_ca, 1)[0],
                np.polyfit(t2m_ca, dato3_ca, 1)[0]), 
              transform=axlb.transAxes,fontsize=14, ha='right')    
    for t in axlb.get_yticklabels():
        t.set_fontsize(12)   
    for t in axlb.get_xticklabels():
        t.set_fontsize(12)       
    axlb.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    axlbb = axlb.twinx()
    lns3 = axlbb.plot(t2m_ca, lw = 2., color = '#ff7f00',
             label = 'T', zorder = 1)
    axlbb.set_xticks([0, 14, 29, 45, 61, 76])
    axlbb.set_xticklabels(['1 June 2010', '', '1 July', '', '1 August', ''])
    axlbb.set_yticklabels([])
    # Ohio
    axrt.set_title('(b) Ohio', fontsize = 16, x = 0.98, ha = 'right', y = 1.03)
    lns1 = axrt.plot(o3_atpoint_north, lw = 2., 
                   c = pollutants_constants.COLOR_CHEMISTRY, 
                   label = 'Transport', zorder = 10)
    lns2 = axrt.plot(dato3_atpoint_north, lw = 2., 
                   c = pollutants_constants.COLOR_TRANSPORT,
                   label = '+ Chemistry', zorder = 5)
    axrt.text(0.96, 0.05, '%.2f, %.2f, %.2f'
              %(np.corrcoef(o3_atpoint_north, t2m_atpoint_north)[0, 1],
                np.polyfit(t2m_atpoint_north, o3_atpoint_north, 1)[0],
                np.polyfit(t2m_atpoint_north, dato3_atpoint_north, 1)[0]), 
              transform=axrt.transAxes,fontsize=14, ha='right')        
    axrt.set_yticklabels([])
    axrtb = axrt.twinx()
    lns3 = axrtb.plot(t2m_atpoint_north, lw = 2., color = '#ff7f00',
             label = 'T', zorder = 1)
    axrtb.set_xticks([0, 14, 29, 45, 61, 76])
    axrtb.set_xticklabels([])
    axrt.tick_params(right = False, left = True, top = True, bottom = False,
                    labelbottom = False, labeltop = True)
    axrtb.set_ylabel('T [K]', rotation = 270, fontsize = 16, color = '#ff7f00')
    axrtb.get_yaxis().set_label_coords(1.2, 0.53)
    for t in axrtb.get_yticklabels():
        t.set_fontsize(12)    
        t.set_color('#ff7f00')
    # Mississippi
    axrb.set_title('(d) Mississippi', fontsize = 16, x = 0.98, ha = 'right', 
                   y = 1.03)
    lns1 = axrb.plot(o3_atpoint, lw = 2., 
                   c = pollutants_constants.COLOR_CHEMISTRY, zorder = 10)
    lns2 = axrb.plot(dato3_atpoint, lw = 2., 
                   c = pollutants_constants.COLOR_TRANSPORT, zorder = 5)
    axrb.set_yticklabels([])
    axrbb = axrb.twinx()
    lns3 = axrbb.plot(t2m_atpoint, lw = 2., color = '#ff7f00', zorder = 1)
    axrbb.set_xticks([0, 14, 29, 45, 61, 76])
    axrbb.set_xticklabels(['1 June 2010', '', '1 July', '', '1 August', ''], 
                         fontsize = 16)
    axrbb.set_ylabel('T [K]', rotation = 270, fontsize = 16, color = '#ff7f00')
    axrbb.get_yaxis().set_label_coords(1.2, 0.53)
    axrb.text(0.96, 0.05, '%.2f, %.2f, %.2f'
              %(np.corrcoef(o3_atpoint, t2m_atpoint)[0, 1],
                np.polyfit(t2m_atpoint, o3_atpoint, 1)[0],
                np.polyfit(t2m_atpoint, dato3_atpoint, 1)[0]), 
              transform=axrb.transAxes,fontsize=14, ha='right')            
    for t in axrb.get_xticklabels():
        t.set_fontsize(12)
    for t in axrbb.get_yticklabels():
        t.set_fontsize(12)    
        t.set_color('#ff7f00')
    # ensure that all axes have the same aspect/scale
    for ax in [axlt, axlb, axrt, axrb]:
        ax.set_ylim([20, 80])    
    for ax in [axltb, axlbb, axrtb, axrbb]:
        ax.set_xlim([0, 91]) 
        ax.set_ylim([285, 315])               
    # Adjust subplots so that a map can be fit in the middle
    plt.subplots_adjust(hspace = 0.3, wspace = 0.5)
    # add inset map showing location of cells
    left, bottom, width, height = [0.445, 0.36, 0.135, 0.26]
    axi = fig.add_axes([left, bottom, width, height])        
    m_small = Basemap(projection = 'merc', llcrnrlon = -126., 
                      llcrnrlat = 24., urcrnrlon = -66.3, 
                      urcrnrlat = 50., resolution = 'h', area_thresh = 1000)
    m_small.drawstates(color = 'k', linewidth = 0.5)
    m_small.drawcountries(color = 'k', linewidth = 1.0)
    m_small.drawcoastlines(color = 'k', linewidth = 1.0)   
    fill_oceans(axi, m_small)   
    x, y = m_small(gmi_lon[slon_idx_south], gmi_lat[slat_idx_south])
    m_small.scatter(x, y, 40, color = 'r', marker = '.', edgecolor = 'r', 
              zorder = 20)
    x, y = m_small(gmi_lon[slon_idx_north], gmi_lat[slat_idx_north])
    m_small.scatter(x, y, 40, color = 'r', marker = '.', edgecolor = 'r', 
              zorder = 20)
    x, y = m_small(gmi_lon[slon_idx_mt], gmi_lat[slat_idx_mt])
    m_small.scatter(x, y, 40, color = 'r', marker = '.', edgecolor = 'r', 
              zorder = 20)
    x, y = m_small(gmi_lon[slon_idx_ca], gmi_lat[slat_idx_ca])
    m_small.scatter(x, y, 40, color = 'r', marker = '.', edgecolor = 'r', 
              zorder = 20)
    leg = axlt.legend(loc = 9, bbox_to_anchor = (1.25, 1.43), fontsize = 16,
                    ncol = 2)
    leg.get_frame().set_linewidth(0.0)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/'+
                'timeseries_transportchemistryo3_atpoints.eps', dpi = 300)
# # # # # # # # # # # # #
def map_allgmio3_do3dt_byprecip(merra_lat, merra_lon, gmi_lat, gmi_lon, 
    t2m_overpass, dat_o3, mr2_o3, emiss_sens, emiss_r):
    """find precipitation in NEUS from MERRA and then calculate dO3/dT 
    on days with high (> 50%-ile) and low (< 50%-ile) regionally-averaged 
    precipitation.

    Parameters
    ----------  
    merra_lat : numpy.ndarray
        MERRA-2 latitude coordinates, units of degrees north, [lat,]
    merra_lon : numpy.ndarray
        MERRA-2 longitude coordinates, units of degrees east, [lon,]       
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
    t2m_overpass : numpy.ndarray
        MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]
    dat_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, for model case 
        HindcastMR2-DiurnalAvgT, units of volume mixing ratio, [time, lat, lon] 
    mr2_o3 : numpy.ndarray
        GMI CTM surface-level ozone at overpass time, for model case 
        HindcastMR2, units of volume mixing ratio, [time, lat, lon] 
    emiss_sens : numpy.ndarray    
        The O3-climate penalty at each GMI grid cell from the + Emissions 
        simulation, units of ppbv K^-1, [lat, lon]         
    emiss_r : numpy.ndarray
        The Pearson product-moment correlation coefficient calculated between 
        MERRA-2 2-meter temperatures and ozone at each GMI grid cell from the 
        + Emissions simulation, [lat, lon]

    Returns
    ----------         
    None    
    """
    import datetime
    import numpy as np
    import netCDF4 as nc
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/globalo3/')
    import globalo3_open
    ncfile = nc.Dataset(pollutants_constants.PATH_CLIMATOLOGIES + 
                        'climatology_2mmy.nc')
    pr_all = ncfile.variables['pr'][:]
    lat = ncfile.variables['lat'][:]
    lon = ncfile.variables['lon'][:]
    # create list of days ('%YYYY-%mm-%dd') in measuring period 
    daysmy = []
    # function uses a generator function to abstract the iteration over the 
    # specified range of dates
    def daterange(start_date, end_date):
        for n in range(int ((end_date - start_date).days)):
            yield start_date + datetime.timedelta(n)
    for year in np.arange(pollutants_constants.START_YEAR, 
                          pollutants_constants.END_YEAR + 1, 1):
        start_date = datetime.date(year, 1, 1)
        end_date = datetime.date(year, 12, 31)
        days_ty = []
        for single_date in daterange(start_date, end_date + 
                                     datetime.timedelta(days = 1)):
            days_ty.append(single_date.strftime("%Y-%m-%d"))
        daysmy.append(days_ty)
    # convert dates from string format to datetime.date type 
    dtimes = []
    dtimes = [datetime.datetime.strptime(date, '%Y-%m-%d').date() for 
              date in np.hstack(daysmy)]
    # months of dtime datetime.date objects
    month_idxs = [day.month for day in dtimes]
    month_idxs = np.array(month_idxs)
    # find indices of summertime months (i.e. June, July, August) and summer 
    # days, years
    summonth_idx = np.where((month_idxs == 6) | (month_idxs == 7) | 
            (month_idxs == 8))
    sumdays = np.array(dtimes)[summonth_idx]    
    sumdays_years = [day.year for day in np.array(dtimes)[summonth_idx]]
    sumdays_years = np.array(sumdays_years)    
    # Find days in measuring period 2008-2010
    wheremp = np.where(np.in1d(sumdays_years, np.array([2008, 2009, 2010])) ==
                       True)[0]
    prmp = pr_all[wheremp]
    # Reduce precipitation to same spatial domain as GMI results 
    prmp_gmi = prmp[:, np.where(lat==merra_lat[0])[0][0]:
        np.where(lat==merra_lat[-1])[0][0]+1, np.where(lon==merra_lon[0])[0][0]:
        np.where(lon==merra_lon[-1])[0][0]+1]
    prmp_lat = lat[np.where(lat==merra_lat[0])[0][0]:
        np.where(lat==merra_lat[-1])[0][0]+1]
    prmp_lon = lon[np.where(lon==merra_lon[0])[0][0]:
        np.where(lon==merra_lon[-1])[0][0]+1]
    # Interpolate MERRA-1 precipitation to GMI resolution 
    prmp_gmi = globalo3_open.interpolate_merra_to_ctmresolution(gmi_lat, 
        gmi_lon, prmp_lat, prmp_lon, prmp_gmi, checkplot='yes') 
    # Recalculate dO3/dT for high precipitation days
    dat_do3dt2m_high = np.empty(shape = t2m_overpass.shape[1:])
    dat_do3dt2m_high[:] = np.nan
    for i, ilat in enumerate(gmi_lat):
        for j, jlon in enumerate(gmi_lon):
            prmp_high = np.where(prmp_gmi[:,i,j] > np.percentile(
                    prmp_gmi[:,i,j], 50))[0]
            if prmp_high.shape[0] == 0:
                dat_do3dt2m_high[i, j] = np.polyfit(t2m_overpass[
                    :, i, j], (dat_o3*1e9)[:, i, j], 1)[0]
            else:                 
                dat_do3dt2m_high[i, j] = np.polyfit(t2m_overpass[
                    prmp_high, i, j], (dat_o3*1e9)[prmp_high, i, j], 1)[0]
    mr2_do3dt2m_high = np.empty(shape = t2m_overpass.shape[1:])
    mr2_do3dt2m_high[:] = np.nan
    for i, ilat in enumerate(gmi_lat):
        for j, jlon in enumerate(gmi_lon):
            prmp_high = np.where(prmp_gmi[:,i,j] > np.percentile(
                    prmp_gmi[:,i,j], 50))[0]   
            if prmp_high.shape[0] == 0:
                mr2_do3dt2m_high[i, j] = np.polyfit(t2m_overpass[
                        :, i, j], (mr2_o3*1e9)[:, i, j], 1)[0]                
            else: 
                mr2_do3dt2m_high[i, j] = np.polyfit(t2m_overpass[
                        prmp_high, i, j], (mr2_o3*1e9)[prmp_high, i, j], 1)[0]
    # For low precipitation days
    dat_do3dt2m_low = np.empty(shape = t2m_overpass.shape[1:])
    dat_do3dt2m_low[:] = np.nan
    for i, ilat in enumerate(gmi_lat):
        for j, jlon in enumerate(gmi_lon):
            prmp_low = np.where(prmp_gmi[:,i,j] < np.percentile(
                    prmp_gmi[:,i,j], 50))[0]
            if prmp_low.shape[0] == 0:
                dat_do3dt2m_low[i, j] = np.polyfit(t2m_overpass[
                        :, i, j], (dat_o3*1e9)[:, i, j], 1)[0]                
            else:
                dat_do3dt2m_low[i, j] = np.polyfit(t2m_overpass[
                        prmp_low, i, j], (dat_o3*1e9)[prmp_low, i, j], 1)[0]
    mr2_do3dt2m_low = np.empty(shape = t2m_overpass.shape[1:])
    mr2_do3dt2m_low[:] = np.nan
    for i, ilat in enumerate(gmi_lat):
        for j, jlon in enumerate(gmi_lon):
            prmp_low = np.where(prmp_gmi[:,i,j] < np.percentile(
                    prmp_gmi[:,i,j], 50))[0]            
            if prmp_low.shape[0] == 0:
                mr2_do3dt2m_low[i, j] = np.polyfit(t2m_overpass[
                        :, i, j], (mr2_o3*1e9)[:, i, j], 1)[0]                
            else: 
                mr2_do3dt2m_low[i, j] = np.polyfit(t2m_overpass[
                        prmp_low, i, j], (mr2_o3*1e9)[prmp_low, i, j], 1)[0]
    # Plotting 
    map_allgmio3_do3dt(dat_do3dt2m_high, mr2_do3dt2m_high, 
        emiss_sens, emiss_r, gmi_lat, gmi_lon, '_highprecip')
    map_allgmio3_do3dt(dat_do3dt2m_low, mr2_do3dt2m_low, 
        emiss_sens, emiss_r, gmi_lat, gmi_lon, '_lowprecip')
    return 
# # # # # # # # # # # # #    
#import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#import sys
#sys.path.append('/Users/ghkerr/phd/')
#import pollutants_constants
#sys.path.append('/Users/ghkerr/phd/GMI/')
#import commensurability 
#years = [2008, 2009, 2010]
## # # # load CASTNet O3 
#castnet = find_conus_castnet(years)
## # # # load MERRA-2 meteorology
#t2m, t10m, u2m, u10m, v2m, v10m, ps, merra_lat, merra_lon, times_all = \
#commensurability.load_MERRA2(years)
## # # # load CTM simulations 
## from Transport simulation
#(gmi_lat, gmi_lon, eta, times, dat_co, dat_no, dat_no2, dat_o3, 
# dat_cloudfraction, dat_gridboxheight) = \
#commensurability.open_overpass2('HindcastMR2-DiurnalAvgT', years)
## from + Chemistry simulation 
#(gmi_lat, gmi_lon, eta, times, mr2_co, mr2_no, mr2_no2, mr2_o3, 
# mr2_cloudfraction, mr2_gridboxheight) = \
# commensurability.open_overpass2('HindcastMR2', years)
## from + Emissions simulation
#(gmi_lat, gmi_lon, eta, times, emiss_co, emiss_no, emiss_no2, emiss_o3, 
# emiss_cloudfraction, emiss_gridboxheight) = \
#commensurability.open_overpass2('GHKerr-DailyEmiss', years)
#gmi_lon = np.mod(gmi_lon - 180.0, 360.0) - 180.0
## # # # determine ozone-temperature sensitivity/correlations
## from CASTNet sites
#(r_castnet, do3dt2m_castnet, lat_castnet, lon_castnet, t_castnet, o3_castnet, 
# sites_castnet) = castnet_r_do3d2t(castnet, t2m, merra_lat, merra_lon, 
# times_all)
## from Transport simulation
#dat_sens, dat_tls, dat_r, dat_t2m_overpass, dat_ps_overpass, dat_p = \
#calculate_gmi_r_do3dt2m(merra_lat, merra_lon, gmi_lat, gmi_lon, t2m, dat_o3, 
#    ps)
## from + Chemistry simulation 
#mr2_sens, mr2_tls, mr2_r, mr2_t2m_overpass, mr2_ps_overpass, mr2_p = \
#calculate_gmi_r_do3dt2m(merra_lat, merra_lon, gmi_lat, gmi_lon, t2m, mr2_o3, 
#    ps)
## from + Emissions simulation
#emiss_sens, emiss_tls, emiss_r, emiss_t2m_overpass, emiss_ps_overpass, emiss_p = \
#calculate_gmi_r_do3dt2m(merra_lat, merra_lon, gmi_lat, gmi_lon, t2m, emiss_o3, 
#    ps)
## # # # calculate regionally-averaged fields
#m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 24.0, 
#            urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'c', 
#            area_thresh = 1000)
#neus_states = pollutants_constants.NORTHEAST_STATES
#neus = find_grid_in_region(m, neus_states, gmi_lat, gmi_lon)
## for CASTNet
#neus_castnet = ['ASH', 'HOW', 'ACA', 'WST', 'HWF', 'ABT', 'WSP', 'CTH', 
#                'MKG', 'KEF', 'PSU', 'ARE', 'LRL', 'CDR', 'PAR', 'VPI', 
#                'PED', 'SHN', 'BWR', 'BEL']
#castnet_sens_neus, castnet_r_neus, castnet_t2m_neus, castnet_o3_neus = \
#calculate_castnet_r_do3dt2m_regionmean(t_castnet, o3_castnet, 
#    sites_castnet, neus_castnet)
## for Transport simulation
#dat_sens_neus, dat_r_neus, dat_t2m_neus, dat_o3_neus = \
#calculate_gmi_r_do3dt2m_regionmean(dat_t2m_overpass, dat_o3, neus, 'Transport')
## for + Chemistry simulation
#mr2_sens_neus, mr2_r_neus, mr2_t2m_neus, mr2_o3_neus = \
#calculate_gmi_r_do3dt2m_regionmean(mr2_t2m_overpass, mr2_o3, neus, '+ Chemistry')
## for + Emissions simulation
#emiss_sens_neus, emiss_r_neus, emiss_t2m_neus, emiss_o3_neus = \
#calculate_gmi_r_do3dt2m_regionmean(emiss_t2m_overpass, emiss_o3, neus, '+ Emissions')
## # # # load AQS MDA8 O3
#sc = list(pollutants_constants.EPA_DICT.values()    
#ozone_mean_mda8, ozone_nomean_mda8, ozone_mda8 = find_conus_aqsmda8(sc)
## correlation coefficients, O3-T2m sensitivity at AQS sites
#r, do3dt2m = aqs_mda8_r_do3dt2m(ozone_nomean_mda8, merra_lat, merra_lon)
## # # # # # # # # # # # #
## visualizations
# focus region, CTM resolution, and mean CASTNet/CTM + Chemistry O3
#map_meanmr2o3meancastneto3_conus(emiss_o3, o3_castnet, gmi_lat, gmi_lon,
#    lat_castnet, lon_castnet, neus_castnet, sites_castnet)
# modeled dO3-dT2m and correlation coefficients
#map_ro3t2m_do3dt2m_conus_gmi(gmi_lat, gmi_lon, dat_sens, dat_r, lat_castnet, 
#    lon_castnet, r_castnet, do3dt2m_castnet, 'Diurnal-AvgT')
#map_ro3t2m_do3dt2m_conus_gmi(gmi_lat, gmi_lon, mr2_sens, mr2_r, lat_castnet, 
#    lon_castnet, r_castnet, do3dt2m_castnet, 'MR2')
#map_ro3t2m_do3dt2m_conus_gmi(gmi_lat, gmi_lon, emiss_sens, emiss_r, 
#    lat_castnet, lon_castnet, r_castnet, do3dt2m_castnet, 'GHKerr-DailyEmiss')
## ratio of O3-T2m sensitivities from Transport and + Chemistry simulations
#map_sensitivityratio_conus(gmi_lat, gmi_lon, dat_r, dat_sens, mr2_sens)
## standard deviations of various fields
#map_std_90ptile(mr2_t2m_overpass, o3, dat_o3, gmi_lat, gmi_lon)
## timeseries and scatterplots of O3, T2m at sites with maximum, median, 
## and minimum sensitivities
#timeseriesscatter_castnetdo3dt2m_cases(do3dt2m_castnet, t_castnet, o3_castnet)
## correlation between O3 from Transport and + Chemistry simulations
#map_rdato3mr2o3_conus(o3, dat_o3, gmi_lat, gmi_lon)
## correlation between O3 from + Chemistry and CTM grid box height
#map_rmr2o3gridboxheight_conus(o3, gridboxheight, gmi_lat, gmi_lon) 
## scatterplots of regionally-averged O3 versus 2-meter temperature and slopes 
## from different regressions
#scatter_castnett2mo3_gmit2mo3_slopes(castnet_o3_neus, castnet_t2m_neus, 
#    emiss_o3_neus, emiss_t2m_neus, 'neus')
## timeseries of regionally-averaged O3 from observations and simulations
#timeseries_castneto3allgmio3(dat_o3_neus, mr2_o3_neus, emiss_o3_neus, 
#    castnet_o3_neus, 'neus', 2010, [2008, 2009, 2010])
# scatterplots and KDES of 2-meter temperature and O3 from different 
# simulations
#scatterhist_castneto3allgmio3(dat_o3_neus, mr2_o3_neus, emiss_o3_neus, 
#    castnet_o3_neus, mr2_t2m_neus, castnet_t2m_neus, 'neus')
## scatterplot of the O3-climate penalty from the Transport vs. + Chemistry 
## simulations with colors corresponding to the correlation coefficient between
## 2-meter temperatures and O3 from the + Emissions simulation
#scatter_dmr2o3dt2m_ddato3dt2m_rt2memisso3_conus(mr2_sens, dat_sens, emiss_r, 
#    gmi_lat, gmi_lon)
## plot percentage contribution from simulations on hot and cold days
#map_allgmio3_hotcold(dat_o3, mr2_o3, emiss_o3, mr2_t2m_overpass, 
#    gmi_lat, gmi_lon, neus)
## timeseries of 2-meter temperatures, CASTNet O3, CEMS NOx
#timeseries_t2m_castneto3_cemsnox(castnet_o3_neus, castnet_t2m_neus, 
#    dat_o3_neus, mr2_o3_neus, emiss_o3_neus, 2010, years, 'neus')
## plot emissions inventory from + Chemistry and + Emissions simulations
## at single grid cell
#NO_inventory_atpoint(mr2_t2m_overpass, 39.2904, 283.387799, 2010, gmi_lat, 
#    gmi_lon)
## plot coefficient of determination between O3 and T2m
#map_r2o3t2m_conus(lat_castnet, lon_castnet, r_castnet, gmi_lat, gmi_lon, 
#    emiss_r, 'GHKerr-DailyEmiss')
# maps of the 90th percentile of O3
#map_allgmio3_p90p10(dat_o3, mr2_o3, emiss_o3, emiss_r, gmi_lat, gmi_lon,
#    neus)
## maps of the O3-climate penalty 
#map_allgmio3_do3dt(dat_sens, mr2_sens, emiss_sens, emiss_r, gmi_lat, 
#    gmi_lon, '')
# plot timeseries of O3, T at grid cells with high and low O3-T correlations
#timeseries_mr2o3dato3t2m_atpoint(mr2_t2m_overpass, dat_o3, mr2_o3, gmi_lat, 
#    gmi_lon)
#timeseries_transportchemistryo3_atpoints(emiss_t2m_overpass, mr2_o3, dat_o3, 
#    gmi_lat, gmi_lon)    
## plot Strode et al. (2015) simulations to compare with + Emissions/
## + Chemistry simulations
#strode_inventory, strode_o3 = scatter_noxo3(mr2_o3, emiss_o3, mr2_no, 
#    emiss_no, mr2_no2, emiss_no2, neus_states, gmi_lat, gmi_lon)
## plot boxplots of CEMS NOx and CASTNet O3
#boxplot_cemsnox_castneto3_neus(neus_castnet, strode_inventory, 
#    strode_o3)
## compare different methods of regional averaging
#compare_regional_averaging(neus, t_castnet, o3_castnet, castnet_t2m_neus, 
#    castnet_o3_neus, sites_castnet, neus_castnet, emiss_t2m_overpass, 
#    emiss_o3, emiss_t2m_neus, emiss_o3_neus, emiss_sens, do3dt2m_castnet)    
## plot comparison of hourly GMI output with overpass2 output
#timeseries_map_hourlyvsoverpass_neus(emiss_o3*1e9, castnet_o3_neus, 
#    emiss_o3_neus, gmi_lat, gmi_lon, 2010, years)   
#timeseries_map_hourlyvsoverpass_neus_nomap(castnet_o3_neus, emiss_o3_neus,
#    2010, years)
#scatter_inventorynoo3(mr2_o3, emiss_o3, neus_states, gmi_lat, gmi_lon)
# # # # Reviewer comments   
# Plot maps of dO3/dT on days with above versus below average precipitation 
#map_allgmio3_do3dt_byprecip(merra_lat, merra_lon, gmi_lat, gmi_lon, 
#    emiss_t2m_overpass, dat_o3, mr2_o3, emiss_sens, emiss_r)