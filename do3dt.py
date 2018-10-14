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
                 zorder = 11, linewidth = 2.))      
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
    None         
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
    where_region[:] = 0.
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
def find_conus_castnet(): 
    """function opens CASTNet observations for the measuring period 2008-2010
    and selects summer (JJA) observations 
    
    Parameters
    ----------    
    None           

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
    for year in [2008, 2009, 2010]:
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
    O3-T2m sensitivity at each CASTNet site. 

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
    afternoon_times = [13]
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
    r : numpy.ndarray
        The Pearson product-moment correlation coefficient calculated between 
        MERRA-2 2-meter temperatures and ozone at each GMI grid cell [lat, lon]
    t2m_atoverpass_interp : numpy.ndarray
        MERRA-2 2-meter temperatures at overpass2 time interpolated to the 
        resolution of the CTM, units of K, [time, lat, lon]
    ps_atoverpass_interp : numpy.ndarray
        MERRA-2 sea level pressure at overpass2 time interpolated to the 
        resolution of the CTM, units of Pa, [time, lat, lon]
    """
    import numpy as np
    import scipy.interpolate
    import datetime
    import pytz
    import timezonefinder
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
                                                      delta_degree = 10)
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
                (XI, YI))
        ps_atoverpass_interp[z] = scipy.interpolate.griddata(
                (X.flatten(), Y.flatten()), ps_atoverpass[z].flatten(), 
                (XI, YI))        
    # calculate the O3-temperature sensitivity at each CTM grid cell 
    print('calculating dO3-dT2m...')
    do3dt2m = np.empty(shape = t2m_atoverpass_interp.shape[1:])
    do3dt2m[:] = np.nan
    for i, ilat in enumerate(gmi_lat):
        for j, jlon in enumerate(gmi_lon):
            do3dt2m[i, j] = np.polyfit(t2m_atoverpass_interp[:, i, j], 
                                       o3[:, i, j], 1)[0]
    # calculate the O3-temperature sensitivity at each CTM grid cell 
    print('calculating correlation coefficients...')
    r = np.empty(shape = t2m_atoverpass_interp.shape[1:])
    r[:] = np.nan
    for i, ilat in enumerate(gmi_lat):
        for j, jlon in enumerate(gmi_lon):
            r[i, j] = np.corrcoef(t2m_atoverpass_interp[:, i, j], 
                                       o3[:, i, j], 1)[0, 1]
    return do3dt2m, r, t2m_atoverpass_interp, ps_atoverpass_interp
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
    llcrnrlat = 25.0
    urcrnrlon = -66.3
    urcrnrlat = 50.
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
        llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat, 
        resolution = 'h', area_thresh = 1000)
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)    
    fig = plt.figure(figsize = (11, 5))
    # O3-T2m correlation map
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
                llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, 
                urcrnrlat = urcrnrlat, resolution = 'h', area_thresh = 1000)
    ax = plt.subplot2grid((1, 2), (0, 0))
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
    cax = fig.add_axes([0.125, 0.25, 0.352, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = r'$\mathregular{\rho}$(O$_{\mathregular{3}}$' + 
                 ', T$_{\mathregular{2 m}}$)', size = 16)    
    cb.ax.tick_params(labelsize = 12)
    cb.set_ticks(np.linspace(vmin, vmax, 6))    
    # O3-T2m sensitivity map         
    ax2 = plt.subplot2grid((1, 2), (0, 1))
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
    cax = fig.add_axes([0.5477, 0.25, 0.352, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_ticks(np.linspace(vmin, vmax, 6))
    cb.set_label(label = '$\mathregular{\partial}$O$_{\mathregular{3}}$ ' + 
                 '$\mathregular{\partial}$T$_{\mathregular{2 m}}^' +
                 '{\mathregular{-1}}$ [ppbv K$^{\mathregular{-1}}$]', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)         
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'map_ro3dt2m_do3dt2m_conus_gmi_%s.eps' %case, dpi = 350)    
    return
# # # # # # # # # # # # #
def map_sensitivityratio_conus(gmi_lat, gmi_lon, dat_sens, mr2_sens):
    """function plot a map of the ratio of the Transport ozone-temperature 
    sensitivity to the + Chemistry ozone-temperature sensitivity (i.e. dO3/
    dT2m). A value of 1 would mean that the sensitivity didn't change between 
    simulations whereas a value of 0.5 would mean that the sensitivity is 
    halved between simulations.
    
    Parameters
    ----------  
    gmi_lat : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]      
    gmi_lon : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lon,]  
    dat_sens : numpy.ndarray    
        The O3-T2m sensitivity at each GMI grid cell in the Transport 
        simulation [lat, lon]   
    dat_sens : numpy.ndarray    
        The O3-T2m sensitivity at each GMI grid cell in the + Chemistry 
        simulation [lat, lon]           

    Returns
    ----------     
    None    
    """
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.colorbar import ColorbarBase
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    # continental U.S. focus region map 
    llcrnrlon = -130.
    llcrnrlat = 24.8
    urcrnrlon = -66.3
    urcrnrlat = 50.
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
                llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, 
                urcrnrlat = urcrnrlat, resolution = 'h', area_thresh = 1000)
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    x, y = np.meshgrid(gmi_lon, gmi_lat)
    x, y = m(x, y)
    # set up colorbar
    cmap = plt.get_cmap('RdBu', 8)
    clevs = np.linspace(0, 2, 9, endpoint = True)
    norm = Normalize(vmin = 0, vmax = 2)
    m.contourf(x, y, (dat_sens/mr2_sens), clevs, cmap = cmap, 
               extend = 'both')
    m.drawstates(color = 'k', linewidth = 0.5)
    m.drawcountries(color = 'k', linewidth = 1.5)
    m.drawcoastlines(color = 'k', linewidth = 1.5)
    # add colorbar
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = '$\partial$O$_{3}$ ' + 
                 '$\partial$T$_{\mathregular{2 m}}^{\mathregular{-1}}$ ratio', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/map_sensitivityratio_conus.eps', 
                dpi = 300)    
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
                 '{\mathregular{O_{\mathregular{3, +\:Chemistry}}}}$', 
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
def timeseries_mr2o3dato3t2m_atpoint(t2m_overpass, o3, dat_o3, gmi_lat, 
    gmi_lon):
    """in the southeastern U.S. the correlations between O3 and temperature 
    and the O3-temperature sensitivity are low; however, the difference between
    the Transport and + Chemistry simulations still can give information about
    the role of transport and deposition on O3 variability. This functions
    plots a time series of O3 from the Transport and + Chemistry simulations
    and T2m for the summer of 2010 at a grid cell in the southeastern U.S. 
    and calculates relevant correlation coefficients and sensitivities. 
    
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
    # indices of 
    slon_idx = 33
    slat_idx = 8
    print('(%.1f, %.1f) is chosen coordinate point...' 
          %(gmi_lat[slat_idx], gmi_lon[slon_idx]))   
    # select values at point
    o3_atpoint = o3[184:, slat_idx, slon_idx] * 1e9
    dato3_atpoint = dat_o3[184:, slat_idx, slon_idx] * 1e9
    t2m_atpoint = t2m_overpass[184:, slat_idx, slon_idx]
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
    # initialize figure, axis 
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    lns1 = ax.plot(o3_atpoint, lw = 2., 
                   c = pollutants_constants.COLOR_CHEMISTRY, 
                   label = 'Transport', zorder = 10)
    lns2 = ax.plot(dato3_atpoint, lw = 2., 
                   c = pollutants_constants.COLOR_TRANSPORT,
                   label = '+ Chemistry', zorder = 5)
    ax.set_xlim([0, len(o3_atpoint) - 1])
    for t in ax.get_yticklabels():
        t.set_fontsize(12)   
    ax.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    # twin axis
    ax2 = ax.twinx()
    lns3 = ax2.plot(t2m_atpoint, lw = 2., color = '#377eb8',
             label = 'T$_{\mathregular{2\:m}}$', zorder = 1)
    ax2.set_xlim([0, len(t2m_atpoint) - 1])
    ax2.set_xticks([0, 14, 29, 45, 61, 76])
    ax2.set_xticklabels(['1 June', '', '1 July', '', '1 August 2010', ''])
    ax2.set_ylabel('T$_{\mathregular{2\:m}}$ [K]', rotation = 270, fontsize = 16)   
    ax2.get_yaxis().set_label_coords(1.15, 0.53)
    for t in ax.get_xticklabels():
        t.set_fontsize(12)
    for t in ax2.get_yticklabels():
        t.set_fontsize(12)    
    # add legend
    lns = lns1 + lns2 + lns3
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc = 9, ncol = 3, frameon = False, 
              fontsize = 16, bbox_to_anchor = (0.52, 1.17))
    plt.subplots_adjust(right = 0.85)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/'
                'timeseries_mr2o3dato3t2m_%dN_%dW_2010.eps' 
                %(gmi_lat[slat_idx], gmi_lon[slon_idx]), dpi = 300)
    return
# # # # # # # # # # # # #    
def map_meanmr2o3meancastneto3_conus(o3, o3_castnet, gmi_lat, gmi_lon,
    lat_castnet, lon_castnet):
    """function plots mean O3 fields from + Chemistry/MR2 simulation and 
    mean CASTNet observations for the measuring period 2008-2010. To illustrate
    the model resolution, parallels and meridians at the CTM's resolution 
    are also plotted. 
    
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
    # show mean JJA O3 as simulated by + Chemistry/HindcastMR2 simulation
    o3_mean = np.mean(o3, axis = 0)
    # from CASTNet
    o3_castnet_mean = np.nanmean(np.array(o3_castnet), axis = 1)
    m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 25., 
            urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'h', 
            area_thresh = 1000)
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
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
    cax = fig.add_axes([0.15, 0.16, 0.73, 0.05])
    cb = ColorbarBase(cax, cmap = cmap, norm = norm, 
                      orientation = 'horizontal', extend = 'both')
    cb.set_label(label = '<O$_{\mathregular{3, +\:Chemistry\:|\:CASTNet}}$>', 
                 size = 16)
    cb.ax.tick_params(labelsize = 12)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/'
                'map_meanmr2o3meancastneto3_conus.eps', dpi = 300)
# # # # # # # # # # # # #    
import numpy as np
import pandas as pd
import sys
sys.path.append('/Users/ghkerr/phd/')
import pollutants_constants
sys.path.append('/Users/ghkerr/phd/GMI/')
import commensurability 
years = [2008, 2009, 2010]
## load AQS MDA8 O3
#sc = list(pollutants_constants.EPA_DICT.values()    
#ozone_mean_mda8, ozone_nomean_mda8, ozone_mda8 = find_conus_aqsmda8(sc)
## load CASTNet O3 
#castnet = find_conus_castnet()
## load MERRA-2 2-meter temperatures
#t2m, t10m, u2m, u10m, v2m, v10m, ps, merra_lat, merra_lon, times_all = \
#commensurability.load_MERRA2(years)
## load CTM simulations from Transport and + Chemistry simulations
#(gmi_lat, gmi_lon, eta, times, co, no, no2, o3, cloudfraction, gridboxheight) = \
#commensurability.open_overpass2('HindcastMR2', years)
#(gmi_lat, gmi_lon, eta, times, dat_co, dat_no, dat_no2, dat_o3, dat_cloudfraction, dat_gridboxheight) = \
#commensurability.open_overpass2('HindcastMR2-DiurnalAvgT', years)
#gmi_lon = np.mod(gmi_lon - 180.0, 360.0) - 180.0
## determine ozone-temperature sensitivity from Transport and + Chemistry 
## simulations
#mr2_sens, mr2_r, mr2_t2m_overpass, mr2_ps_overpass = \
#calculate_gmi_r_do3dt2m(merra_lat, merra_lon, gmi_lat, gmi_lon, t2m, o3, ps)
#dat_sens, dat_r, dat_t2m_overpass, dat_ps_overpass = \
#calculate_gmi_r_do3dt2m(merra_lat, merra_lon, gmi_lat, gmi_lon, t2m, dat_o3, ps)
## correlation coefficients, O3-T2m sensitivity at AQS sites
#r, do3dt2m = aqs_mda8_r_do3dt2m(ozone_nomean_mda8, merra_lat, merra_lon)
## correlation coefficients, O3-T2m sensitivity at CASTNet sites
#r_castnet, do3dt2m_castnet, lat_castnet, lon_castnet, t_castnet, o3_castnet, sites_castnet = \
#castnet_r_do3d2t(castnet, t2m, merra_lat, merra_lon, times_all)
# # # # # # # # # # # # #
# visualizations
## modeled dO3-dT2m and correlation coefficients
#map_ro3t2m_do3dt2m_conus_gmi(gmi_lat, gmi_lon, mr2_sens, mr2_r, lat_castnet, 
#    lon_castnet, r_castnet, do3dt2m_castnet, 'MR2')
#map_ro3t2m_do3dt2m_conus_gmi(gmi_lat, gmi_lon, dat_sens, dat_r, lat_castnet, 
#    lon_castnet, r_castnet, do3dt2m_castnet, 'Diurnal-AvgT')
## ratio of O3-T2m sensitivities from Transport and + Chemistry simulations
#map_sensitivityratio_conus(gmi_lat, gmi_lon, dat_sens, mr2_sens)
## standard deviations of various fields
#map_std_90ptile(mr2_t2m_overpass, o3, dat_o3, gmi_lat, gmi_lon)
## timeseries and scatterplots of O3, T2m at sites with maximum, median, 
## and minimum sensitivities
#timeseriesscatter_castnetdo3dt2m_cases(do3dt2m_castnet, t_castnet, o3_castnet)
## correlation between O3 from Transport and + Chemistry simulations
#map_rdato3mr2o3_conus(o3, dat_o3, gmi_lat, gmi_lon)
## correlation between O3 from + Chemistry and CTM grid box height
#map_rmr2o3gridboxheight_conus(o3, gridboxheight, gmi_lat, gmi_lon) 
## time series of O3 and T2m at a grid point with negative O3-T2m correlation 
## and sensitivity
#timeseries_mr2o3dato3t2m_atpoint(mr2_t2m_overpass, o3, dat_o3, gmi_lat, 
#    gmi_lon)

#import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap    
#
#m = Basemap(projection = 'merc', llcrnrlon = -130., llcrnrlat = 24.8, 
#            urcrnrlon = -66.3, urcrnrlat = 50., resolution = 'c', 
#            area_thresh = 1000)
#neus_states = pollutants_constants.NORTHEAST_STATES
#
#
#neus = find_grid_in_region(m, neus_states, gmi_lat, gmi_lon)






