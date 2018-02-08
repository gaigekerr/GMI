#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAME
    gmi_validation.py
PURPOSE
    Script opens simulations and observations of atmospheric trace gases 
    for a specified region (here specified by GMI sites) and extracts 
    daytime (1100-1600 hour, local time) observations and produces daytime mean
    regionally-averaged concentrations for the measuring period 2005-2010
PROGRAMMER
    Gaige Hunter Kerr
REVISION HISTORY
    07022018 -- initial version created
    08022018 -- functions 'find_commensurate_t2m' and 'load_t2m' added
    """
def open_castnet_gmi(case, year, sampling_months, sampling_hours):
    """function opens hourly CASTNET ozone observations for the specified 
    year and finds commensurate GMI column ozone concentrations
    
    Parameters
    ----------    
    case : str
        Hindcast family (i.e. MR2, MR2-CCMI, FFIgac2, etc.)
    year : int
        Year of interest
    sampling_months : list 
        Months of interest
    sampling_hours : list 
        Hours during which ozone concentrations are fetched

    Returns
    ----------      
    castnet : pandas.core.frame.DataFrame
        Hourly ozone concentrations at CASTNET stations for the specified 
        years, months, and hours, [no. obs, 5]
    o3 : numpy.ndarray
        Hourly column ozone concentrations at GMI stations, units of volume
        mixing ratio, [no. GMI stations, sampling hours * days in period, 
        pressure]    
    co : numpy.ndarray
        Same as o3 but for carbon monoxide (CO)
    no : numpy.ndarray
        Same as o3 but for nitric oxide (NO)
    no2 : numpy.ndarray         
        Same as o3 but for nitrogen dioxide (NO2)
    gmi_sites : numpy.ndarray
        Three letter codes of GMI station names, [no. GMI stations,]
    gmi_lat : numpy.ndarray
        Latitude of GMI stations, degrees north, [no. GMI stations,]
    gmi_lon : numpy.ndarray
        Longitude of GMI stations, degrees east, [no. GMI stations,]
    times_ty : numpy.ndarray        
        datetime.datetime objects for hours in measuring period, [sampling 
        hours * days in period]
    """
    import numpy as np
    import pandas as pd
    from datetime import datetime, timedelta
    from netCDF4 import Dataset    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # # # # open CASTNET O3 observations for year of interest
    castnet = pd.read_csv((pollutants_constants.PATH_CASTNET + 
                            '%s/ozone_%s.csv' %(year, year)), header = 0,
                            low_memory = False)
    # columns in CASTNET dataframe
    castnet_cols = ['SITE_ID',
                    'DATE_TIME',
                    'OZONE',
                    'QA_CODE',
                    'UPDATE_DATE']
    castnet = pd.DataFrame(castnet, columns = castnet_cols)
    # convert 'DATE_TIME' column from pandas.core.series.Series to a DateTime 
    # object 
    castnet['DATE_TIME'] = pd.to_datetime(castnet['DATE_TIME'])
    # consider only summertime ozone observations (JJA)
    castnet = castnet.loc[castnet['DATE_TIME'].dt.month.isin(sampling_months)]
    # sample CASTNET observations only during sampling_hours (i.e. in military
    # time format); in this case it's during the daytime (11-16 LST), when 
    # the boundary layer is well-mixed
    castnet = castnet.loc[castnet['DATE_TIME'].dt.hour.isin([x - 4 
                          for x in sampling_hours])]
    # remove the last three letters of SITE_ID column for direct comparison 
    # with Cooper/GMI station names 
    castnet['SITE_ID'] = castnet['SITE_ID'].astype(str).str[:-3].astype(np.str)
    castnet['DATE_TIME'] = pd.DatetimeIndex(castnet['DATE_TIME'])
    # # # # open GMI hourly column O3 concentrations for year of interest 
    infile = Dataset(pollutants_constants.PATH_GMI + '%s/%d/' %(case, year)
         + '%s_%d.nc' %(case, year), 'r')
    o3 = infile.variables['o3'][:]
    co = infile.variables['co'][:]
    no = infile.variables['no'][:]
    no2 = infile.variables['no2'][:]
    gmi_sites = infile.variables['sites'][:]
    gmi_lat = infile.variables['lat'][:]
    gmi_lon = infile.variables['lon'][:]
    # # # # find indicies of relevant dates/hours
    def perdelta(start, end, delta):
        curr = start
        while curr < end:
            yield curr
            curr += delta
    times_ty = []
    # only summer model output was downloaded for Replay simulations, so 
    # adjust times_ty to reflect just JJA 
    if case == 'M2G_c90' or case == 'EGU_T':
        for result in perdelta(datetime(year, 6, 1, 0), 
                               datetime(year, 9, 1, 0), 
                               timedelta(hours = 1)):
            times_ty.append(result)     
    # loop through days in measuring period (year of interest) with an hourly
    # timestep
    else: 
        for result in perdelta(datetime(year, 1, 1, 0), 
                               datetime(year + 1, 1, 1, 0), 
                               timedelta(hours = 1)):
            times_ty.append(result)                 
    # find positions which overlap with the sampling months and sampling hours
    # of interest for commensurate comparison with CASTNET observations
    month_idxs = [day.month for day in times_ty]
    month_idxs = np.array(month_idxs)
    hour_idxs = [day.hour for day in times_ty]
    hour_idxs = np.array(hour_idxs)
    poi1 = np.in1d(month_idxs, sampling_months).reshape(month_idxs.shape)
    poi2 = np.in1d(hour_idxs, sampling_hours).reshape(hour_idxs.shape)
    # poi is an array of indicies where hours are within sampling_hours 
    # and months are within sampling_months specified
    poi = np.where((poi1 == True) & (poi2 == True))[0]
    times_ty = np.array(times_ty)[poi]
    # index O3 concentrations for relevant dates/times
    if case == 'M2G_c90' or case == 'EGU_T':
        o3 = o3[:, poi]
        co = co[:, poi]
        no = no[:, poi]
        no2 = no2[:, poi]
    else:
        o3 = o3[:, poi, :] 
        co = co[:, poi, :] 
        no = no[:, poi, :] 
        no2 = no2[:, poi, :] 
    return castnet, o3, co, no, no2, gmi_sites, gmi_lat, gmi_lon, times_ty
# # # # # # # # # # # # #
def gmi_castnet_dailyavg_singleyear(gmi_siteid_fr, case, year, regavg):
    """function opens hourly output (CO, NO, NO2, and O3) from a specified 
    GMI CTM hindcast simulation during a specified summer and finds co-located 
    O3 concentrations from CASTNet (if more than one CASTNet station is 
    associated with a particular GMI site, these observations are averaged). 
    From the hourly output, ozone is fetched for the afternoon (11am - 4pm) 
    when the boundry layer is well-mixed, and these hours are averaged to 
    produce a daily average. Also returned are daily, regionally-averaged 
    timeseries of CO, NO, and NO2 at all GMI sites in focus region (not just 
    at co-located CASTNet sites). 

    Parameters
    ----------    
    gmi_sites_fr : list
        GMI site names in region
    case : str
        Hindcast family (i.e. MR2, MR2-CCMI, FFIgac2, etc.)
    year : int
        Year of interest        
    regavg : str
        If 'yes', then a regional average of trace gases will be calculated 
        over all GMI sites in focus region 

    Returns
    ----------       
    castnet_o3 : list
        All available complete JJA ozone measurements from CASTNet sites in 
        'castnet_sites_fr' during the specified summer averaged over the 6 hour
        daily sampling window, units of ppbv, [stations with complete records, 
        sampling hours * days in summer]       
    gmi_o3 : list
        JJA GMI O3 concentrations at NDACC sites which are co-loated with 
        CASTNET sites with complete records averaged over the 6 hour daily 
        sampling window, units of ppbv, [stations with complete records, 
        sampling hours * days in summer]
    gmi_co : list
        Same as gmi_o3 but for GMI CO, units of ppbv
    gmi_no : list
        Same as gmi_o3 but for GMI NO, units of ppbv
    gmi_no2 : list
        Same as gmi_o3 but for GMI NO2, units of ppbv
    gmi_lat_at : list
        Latitude of GMI sites with corresponding CASTNet observations with 
        complete observations for summer of interest, 
        [stations with complete records,]
    gmi_lon_at : list 
        Longitude of GMI sites with corresponding CASTNet observations with 
        complete observations for summer of interest, 
        [stations with complete records,]
    castnet_lat_at : list
        Latitude of CASTNet sites with complete observations for the specified
        summer, [stations with complete records,]
    castnet_lon_at : list 
        Longitude of CASTNet sites with complete observations for the specified
        summer, [stations with complete records,]  
    """
    import numpy as np
    import pandas as pd
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # # # # GMI/CASTNET site match-up information created with Strode et al. 
    # [2015] supplementary information (Table S1) 
    siteinfo = np.genfromtxt(pollutants_constants.PATH_GMISTATIONINFO + 
                             'gmi_castnet_matchup.csv', dtype = 'str', 
                             skip_header = 1, delimiter = ',')
    # columns for site information are 
    # [:, 0] - GMI sites (from GMI file naming convention)
    # [:, 1] - Corresponding CASTNet site names
    # [:, 2] - Timezone (string, e.g. EST, MST) 
    # [:, 3] - Sitenum 
    # [:, 4] - GMIaltbox (where to sample idaily, etc GMI grid w.r.t. pressure)
    castnet_sites = list(siteinfo[:, 1])
    gmi_sites = list(siteinfo[:, 0])     
    # # # # CASTNet site info
    dirpath = pollutants_constants.PATH_CASTNET
    castnet_siteinfo = pd.read_csv((dirpath + 'site.csv'), header = 0)
    cols = ['SITE_ID',
            'SITE_NUM',
            'SITE_NAME',
            'ACTIVE',
            'INACTIVE',
            'AGENCY',
            'STATE',
            'COUNTY',
            'TIME_ZONE',
            'LATITUDE',
            'LONGITUDE',
            'ELEVATION',
            'MAPID',
            'LAND_USE',
            'TERRAIN',
            'MLM',
            'NADP_ID',
            'NADP_DISTANCE',
            'UPDATE_DATE']
    castnet_siteinfo = pd.DataFrame(castnet_siteinfo, columns = cols)
    # remove trailing numbers from SITE_ID
    castnet_siteinfo['SITE_ID'] = castnet_siteinfo['SITE_ID'].map(lambda x: str(x)[:3])
    # list creation for concentrations, site locations
    gmi_co, gmi_no, gmi_no2, gmi_o3 = [], [], [], []    
    castnet_o3 = []
    gmi_lat_at, gmi_lon_at = [], []
    castnet_lat_at, castnet_lon_at = [], []
    # open files
    sampling_months = [6, 7, 8]
    sampling_hours = [15, 16, 17, 18, 19, 20]
    castnet_ty, o3_ty, co_ty, no_ty, no2_ty, sites, gmi_lat, gmi_lon, \
    times_ty = open_castnet_gmi(case, year, sampling_months, sampling_hours)
    # loop through GMI sites in focus region 
    for site in gmi_siteid_fr: 
        # month index for monthly mean GMI lookup 
        month_idxs = [day.month for day in times_ty]
        month_idxs = np.array(month_idxs)
        month_where = np.where((month_idxs == 6) |
                               (month_idxs == 7) |
                               (month_idxs == 8))[0]   
        # find index/indices of GMI site in 'gmi_castnet_matchup.csv'
        gmi_site_pos = np.where(np.array(gmi_sites) == site)[0]
        # find the GMI altitude grid box for GMI site; if GMI site corresponds to 
        # more than one CASTNet site, then find average altitude grid box (rounded
        # to nearest integer)
        gmi_site_alt = siteinfo[gmi_site_pos, 4]   
        gmi_site_alt = gmi_site_alt.tolist()
        gmi_site_alt = [int(x) for x in gmi_site_alt]
        gmi_site_alt = int(np.mean(gmi_site_alt))
        # extract CO, NO, NO2, and O3 (n.b. for CO, NO, and NO2 these time 
        # series are saved off for all GMI in focus region since we're not 
        # concerned about GMI being co-located with GMI sites since most sites 
        # don't measure these trace gases) and convert to ppbv
        # for CO
        co_ty_as = co_ty[gmi_site_pos[0], :, gmi_site_alt]
        co_ty_as = co_ty_as[month_where]
        # average every 6 hours (i.e. 11 - 16 local time) for a daily average    
        co_ty_as = np.mean(np.array(co_ty_as).reshape(-1, 6), axis = 1)
        # append to yearly files
        gmi_co.append(co_ty_as * 1e9)    
        # for NO
        no_ty_as = no_ty[gmi_site_pos[0], :, gmi_site_alt]
        no_ty_as = no_ty_as[month_where]
        no_ty_as = np.mean(np.array(no_ty_as).reshape(-1, 6), axis = 1)    
        gmi_no.append(no_ty_as * 1e9)  
        # for NO2
        no2_ty_as = no2_ty[gmi_site_pos[0], :, gmi_site_alt]
        no2_ty_as = no2_ty_as[month_where]
        no2_ty_as = np.mean(np.array(no2_ty_as).reshape(-1, 6), axis = 1)    
        gmi_no2.append(no2_ty_as * 1e9)  
        # for O3
        o3_ty_as = o3_ty[gmi_site_pos[0], :, gmi_site_alt]
        o3_ty_as = o3_ty_as[month_where]          
        # find CASTNet sites ID(s) co-located with GMI site
        castnet_SITE_ID = np.array(castnet_sites)[gmi_site_pos]
        # find CASTNet observations corresponding to CASTNet site ID(s)
        castnet_ty_atgmi = castnet_ty.loc[castnet_ty['SITE_ID'].isin(castnet_SITE_ID)]
        # extract CASTNet lat/lon
        if castnet_ty_atgmi.shape[0] > 551: 
            castnet_siteinfo_atgmi = castnet_siteinfo.loc[
                    castnet_siteinfo['SITE_ID'].isin(castnet_SITE_ID)][['LATITUDE', 'LONGITUDE']]
            castnet_siteinfo_atgmi = castnet_siteinfo_atgmi.drop_duplicates()
            castnet_lat_at.append(castnet_siteinfo_atgmi['LATITUDE'].values)
            castnet_lon_at.append(castnet_siteinfo_atgmi['LONGITUDE'].values)   
        # sort and find daily average, if multiple CASTNet sites correspond
        # to a single GMI site
        castnet_ty_atgmi = castnet_ty_atgmi.sort_values(by = 'DATE_TIME')
        castnet_ty_atgmi = castnet_ty_atgmi.groupby(['DATE_TIME']).mean()
        # only save off CASTNet and GMI O3 for a particular site if complete 
        # observations exist; extract lat/lon of GMI sites considered for year
        # of interest
        if o3_ty_as.shape[0] == castnet_ty_atgmi.shape[0]:
            gmi_lat_at.append(gmi_lat[gmi_site_pos])
            gmi_lon_at.append(gmi_lon[gmi_site_pos])
            castnet_ty_atgmi = np.mean(np.array(castnet_ty_atgmi['OZONE'].values).reshape(-1, 6), axis = 1)
            o3_ty_as = np.mean(np.array(o3_ty_as).reshape(-1, 6), axis = 1) 
            castnet_o3.append(castnet_ty_atgmi)
            gmi_o3.append(o3_ty_as * 1e9)
#    # # # # OPTIONAL, plot location of GMI, CASTNet sites
#    import matplotlib.pyplot as plt
#    from mpl_toolkits.basemap import Basemap
#    fig = plt.figure()
#    ax = plt.subplot2grid((1, 1), (0, 0), rowspan = 1, colspan = 1)
#    ax.set_title(year)
#    m = Basemap(projection = 'cass', llcrnrlon = -85., llcrnrlat = 36, urcrnrlon = -63.,
#                urcrnrlat = 50., resolution = 'i', lat_0 = 20, lon_0 = -88, 
#                area_thresh = 10000)      
#    x_castnet, y_castnet = m(np.hstack(castnet_lon_at), np.hstack(castnet_lat_at))
#    x_gmi, y_gmi = m(np.hstack(gmi_lon_at), np.hstack(gmi_lat_at))
#    mgmi = m.plot(x_gmi, y_gmi, 'ro', markersize = 12, label = 'GMI CTM')  
#    mcastnet = m.plot(x_castnet, y_castnet, 'kx', markersize = 12, label = 'CASTNet') 
#    plt.legend()
#    m.drawstates()
#    m.drawcountries()
#    m.drawcoastlines()
#    # # # # 
    if regavg == 'yes': 
        castnet_o3 = np.nanmean(np.array(castnet_o3), axis = 0)
        gmi_o3 = np.mean(np.array(gmi_o3), axis = 0)
        gmi_co = np.mean(np.array(gmi_co), axis = 0)
        gmi_no = np.mean(np.array(gmi_no), axis = 0)
        gmi_no2 = np.mean(np.array(gmi_no2), axis = 0)
    return castnet_o3, gmi_o3, gmi_co, gmi_no, gmi_no2, gmi_lat_at, gmi_lon_at, castnet_lat_at, castnet_lon_at
# # # # # # # # # # # # #
def gmi_castnet_dailyavg_multiyear(case, gmi_sites_fr, regavg):
    """uses gmi_castnet_dailyavg_singleyear to open observed and modeled 
    trace gases for a specified range of years.
    
    Parameters
    ----------    
    case : str
        Hindcast family (i.e. MR2, MR2-CCMI, FFIgac2, etc.)
    gmi_sites_fr : list
        GMI site names in region
    regavg : str
        If 'yes', then a regional average of trace gases will be calculated 
        over all GMI sites in focus region 

    Returns
    ----------       
    castnet_o3_all : list
        CASTNet stations with complete JJA O3 measurements in the specified 
        focus region for the specified years, units of ppbv, if regavg is 'yes'
        [number of years, days in summer] and if 'no' [number of years, 
        number of GMI sites, days in summer]
    gmi_o3_all : list 
        Same as castnet_o3_all but for GMI O3 at GMI sites with complete 
        CASTNet observations, units of ppbv
    gmi_co_all : list 
        Same as gmi_o3_all but for GMI CO at all GMI sites in focus region, 
        units of ppbv    
    gmi_no_all : list 
        Same as gmi_o3_all but for GMI NO at all GMI sites in focus region, 
        units of ppbv    
    gmi_no2_all : list 
        Same as gmi_o3_all but for GMI NO2 at all GMI sites in focus region, 
        units of ppbv      
    gmi_lat_at_all : list
        Latitude of GMI sites with corresponding CASTNet observations with 
        complete observations
    gmi_lon_at_all : list
        Longitude of GMI sites with corresponding CASTNet observations with 
        complete observations
    castnet_lat_at_all : list
        Latitude of CASTNet sites with complete observations
    castnet_lon_at_all : list  
        Longitude of CASTNet sites with complete observations    
    """                              
    # creation of multi-year lists 
    castnet_o3_all = []
    gmi_o3_all = []
    gmi_co_all = []
    gmi_no_all = []
    gmi_no2_all = []
    gmi_lat_at_all = []
    gmi_lon_at_all = []
    castnet_lat_at_all = []
    castnet_lon_at_all = []
    # loop through years in measuring period
    for year in [2005, 2006, 2007, 2008, 2009, 2010]: 
        print('extracting trace gases for %s!' %year)
        # open observed/modeled trace gases for year
        castnet_o3, gmi_o3, gmi_co, gmi_no, gmi_no2, gmi_lat_at, \
        gmi_lon_at, castnet_lat_at, castnet_lon_at = gmi_castnet_dailyavg_singleyear(gmi_siteid_fr, case, year, 'yes')
        # append to multi-year list
        castnet_o3_all.append(castnet_o3)
        gmi_o3_all.append(gmi_o3)
        gmi_co_all.append(gmi_co)
        gmi_no_all.append(gmi_no)
        gmi_no2_all.append(gmi_no2)
        gmi_lat_at_all.append(gmi_lat_at)
        gmi_lon_at_all.append(gmi_lon_at)
        castnet_lat_at_all.append(castnet_lat_at)
        castnet_lon_at_all.append(castnet_lon_at)        
    return castnet_o3_all, gmi_o3_all, gmi_co_all, gmi_no_all, gmi_no2_all, \
    gmi_lat_at_all, gmi_lon_at_all, castnet_lat_at_all, castnet_lon_at_all
# # # # # # # # # # # # #
def find_summer_days():
    """function finds all summer (1 June - 31 August) days in the measuring 
    period 2000 - 2014.

    Parameters
    ----------  
    None
    
    Returns
    ----------    
    daysmy : pandas.tseries.index.DatetimeIndex
        Summer (1 June - 31 August) days in measuring period
    """
    import datetime
    import pandas as pd
    import numpy as np
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants     
    # create list of days in measuring period 
    daysmy = []
    # function uses a generator function to abstract the iteration over the specified range of dates
    def daterange(start_date, end_date):
        for n in range(int ((end_date - start_date).days)):
            yield start_date + datetime.timedelta(n)
    for year in np.arange(pollutants_constants.START_YEAR, 
                          pollutants_constants.END_YEAR + 1, 1):
        start_date = datetime.date(year, 1, 1)
        end_date = datetime.date(year, 12, 31)
        days_ty = []
        for single_date in daterange(start_date, end_date + datetime.timedelta(days = 1)):
            days_ty.append(single_date)
        daysmy.append(days_ty)
    daysmy = pd.to_datetime(np.hstack(daysmy))
    daysmy = daysmy[np.where((daysmy.month == 6) | (daysmy.month == 7) | (daysmy.month == 8))]
    return daysmy    
# # # # # # # # # # # # #   
def load_t2m(): 
    """function opens 2m (surface-level) MERRA-1 temperature climatologies 
    (daily mean temperature)
    
    Parameters
    ----------    
    None

    Returns
    ----------    
    t2m_all : numpy.ndarray
        MERRA-1 daily mean 2m temperature over the North American domain, 
        units of K, [time, lat, lon]
    lat : numpy.ndarray
        MERRA-1 latitude spine, units of degrees north, [lat,]    
    lon : numpy.ndarray
        MERRA-1 longitude spine, units of degrees east, [lon,]
    """
    # open 2 m and surface-level climatologies
    import numpy as np
    import netCDF4 as nc
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants    
    ncfile = nc.Dataset(pollutants_constants.PATH_CLIMATOLOGIES + 
                        'climatology_2mmy.nc')
    t2m_all = ncfile.variables['t2m'][:]
    lat = ncfile.variables['lat'][:]
    lon = ncfile.variables['lon'][:]
    # reduce domain to be consisent with stagnation 
    llat_reduced_idx, ulat_reduced_idx = np.abs(lat - 20.).argmin(), np.abs(lat - 55.).argmin()
    rlon_reduced_idx, llon_reduced_idx = np.abs(lon + 130.).argmin(), np.abs(lon + 60.).argmin()
    lat = lat[llat_reduced_idx:ulat_reduced_idx]
    lon = lon[rlon_reduced_idx: llon_reduced_idx]
    t2m_all = t2m_all[:, llat_reduced_idx:ulat_reduced_idx, 
                     rlon_reduced_idx: llon_reduced_idx]
    ncfile.close()
    # find summer days in measuring period which spans T2m record in 
    # climatology file
    daysmy = find_summer_days()
    # n.b. We have GMI CTM output from 2005 - 2010, select only these years
    in_gmi_period = np.where((daysmy.year == 2005) |
                             (daysmy.year == 2006) |
                             (daysmy.year == 2007) |
                             (daysmy.year == 2008) |
                             (daysmy.year == 2009) |
                             (daysmy.year == 2010))[0]
    t2m_all = t2m_all[in_gmi_period]
    return t2m_all, lat, lon
# # # # # # # # # # # # #
def find_commensurate_t2m(castnet_lon_at_all, castnet_lat_at_all, regavg): 
    """function loads MERRA 2 meter temperatures for measuring period 2005-
    2010 and finds MERRA grid cells which are co-located with CASTNet sites
    in focus region. Two meter temperature at these grid cells is then found
    and, if desired, averaged to produce a regional average temperature.
    
    Parameters
    ----------    
    castnet_lon_at_all : list  
        Longitudes of CASTNet sites with complete observations for measuring 
        period, [years in measuring period] 
    castnet_lat_at_all : list
        Latitudes of CASTNet sites with complete observations for measuring 
        period, [years in measuring period
    regavg : str
        If 'yes', then a regional average of trace gases will be calculated 
        over all GMI sites in focus region 

    Returns
    ----------
    t2m_fr : numpy.ndarray
        MERRA 2 meter temperatures nearest to every CASTNet site in focus 
        region or regionally-averaged
    """
    import numpy as np
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    # load MERRA 2 meter temperature
    t2m_all, lat, lon = load_t2m()
    # stack all years of CASTNet lat/lon together 
    castnet_lon_at_all = np.hstack(castnet_lon_at_all)
    castnet_lon_at_all = np.hstack(castnet_lon_at_all)
    castnet_lat_at_all = np.hstack(castnet_lat_at_all)
    castnet_lat_at_all = np.hstack(castnet_lat_at_all)    
    lat_idxS, lon_idxS = [], []
    # find MERRA find grid cells closet to CASTNet sites
    for ilat, ilon in zip(list(castnet_lat_at_all), list(castnet_lon_at_all)): 
            lat_idx = geo_idx(ilat, lat)
            lon_idx = geo_idx(ilon, lon)
            lat_idxS.append(lat_idx)
            lon_idxS.append(lon_idx)
    # merge lists of lat/lon indices, zip returns species_l, a list of tuples
    species_l = zip(lat_idxS, lon_idxS)
    # get unique values in species_l by converting the list to a set
    species_l = list(set(species_l))
    species_l = np.array(species_l)
    # select MERRA data closest to CASTNet sites
    t2m_fr = t2m_all[:, species_l[:, 0], species_l[:, 1]]
    lat_fr = lat[species_l[:, 0]]
    lon_fr = lon[species_l[:, 1]]
#    # # # # OPTIONAL, plot location of CASTNet sites, MERRA grid cells
#    import matplotlib.pyplot as plt
#    from mpl_toolkits.basemap import Basemap
#    fig = plt.figure()
#    ax = plt.subplot2grid((1, 1), (0, 0), rowspan = 1, colspan = 1)
#    m = Basemap(projection = 'cass', llcrnrlon = -85., llcrnrlat = 36, urcrnrlon = -63.,
#                urcrnrlat = 50., resolution = 'i', lat_0 = 20, lon_0 = -88, 
#                area_thresh = 10000)      
#    x_castnet, y_castnet = m(np.hstack(castnet_lon_at_all), 
#                             np.hstack(castnet_lat_at_all))
#    x_merra, y_merra = m(lon_fr, lat_fr)
#    mmmera = m.plot(x_merra, y_merra, 'rs', markersize = 12, label = 'MERRA')  
#    mcastnet = m.plot(x_castnet, y_castnet, 'kx', markersize = 12, label = 'CASTNet') 
#    plt.legend()
#    m.drawstates()
#    m.drawcountries()
#    m.drawcoastlines()
#    # # # # 
    # calculate regional average
    if regavg == 'yes': 
        t2m_fr = np.mean(t2m_fr, axis = 1)
    return t2m_fr
# # # # # # # # # # # # #