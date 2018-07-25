#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAME
    commensurability.py
PURPOSE
    Script opens hourly CTM output, trace gas observations, and meteorology
    and finds commensurate concentrations
PROGRAMMER
    Gaige Hunter Kerr
REVISION HISTORY
    05042018 -- initial version created
    06042018 -- functions 'open_castnetsiting' and 'commensurate_castnet_gmi' 
    08042018 -- functions 'find_summer_days', 'load_t2m', and 'commensurate_t2m'
                added
    09042018 -- functions 'commensurate_aqsno2co', 'get_merged_csv', and 
                'load_aqsno2co' added
    13042018 -- functions 'commensurate_castnet_gmi_diurnal' and 
                'commensurate_aqsno2co_diurnal' added to extract diurnal cycle 
                from observed and modeled trace gas concentrations 
    14042018 -- function 'commensurate_castnet_gmi_locations' added and 
                'commensurate_aqsno2co_diurnal' edited such that the locations 
                (latitude/longitude coordinates) of AQS stations which measure
                NO2 and CO are extracted and returned
    16042018 -- function 'commensurate_castnet_gmi_diurnal_ra' added; 
                this function uses 'commensurate_castnet_gmi_locations'
                to open observed/modelled CASTNet/GMI concentrations from all 
                GMI CTM model configurations. Also function 
                'commensurate_castnet_gmi_ra' added
    22042018 -- functions 'open_perturbed_emissions' and 
                'commensurate_emiss_perturbed' added; these functions open 
                perturbed emissions inventories for GMI CTM and find emissions 
                at grid cells co-located with CASTNet sites
    02052018 -- issue discovered with CASTNet data as some stations/years 
                report multiple measurements per hourly observation. This is 
                addressed in function 'commensurate_castnet_gmi' by averaging 
                these duplicate measurements
    02052018 -- bug discovered in function 'open_castnet_singyear', when
                function read in yearly CASTNet files it failed to account for
                the column 'OZONE_F.' This should not have modified any 
                results, but screwed up pandas' reading of the dataframe's 
                final two columns (i.e. 'QA_OZONE' and 'UPDATE_DATE')
    09052018 -- function 'commensurate_geos_gmi' added to extract timeseries 
                of trace gases from Katie Travis' GEOS-Chem simulations at 
                CASTNet sites
    10052018 -- bug discovered in 'commensurate_aqsno2co;' before it was not
                sampling hourly vales at hours of interest, rather it was 
                taking the mean over all hours in diurnal cycle. 
    16052018 -- functions 'commensurate_geos_gmi' and 
                'commensurate_geos_gmi_diurnal' added
    18052018 -- function 'load_aqsno2co' modified so that (1) it only reads in 
                hourly observations for years of interest to save computing 
                time and (2) also read in O3. Name of function changed to 
                'load_aqshourly'
    01062018 -- function 'commensurate_aqstracegas_diurnal' modified so that 
                only unique latitude/longitude coordinates pairs for AQS NO2, 
                CO, and O3 sites are returned. Previously all unique latitude
                and longitudes (not pairs) had been returned
    04062018 -- previously function 'open_castnet_singyear' had converted 
                local time to UTC/GMT by subtracting off the UTC offset from 
                the local times; however this caused an error for hours between
                O UTC and the UTC offset (i.e. if UTC offset was -4, then 
                O-4 local time would become -4 to 0 UTC). This was changed 
                so that times were converted using pytz
    24062018 -- function 'commensurate_aqstracegas_siting' added
    05072018 -- functions 'open_gridded_idailyCTM' and 
                'commensurate_aqstracegas_gridded' added to examine gridded
                12Z HindcastMR2 output
    06072018 -- function 'commensurate_aqstracegas_gridded' changed to take 
                siting environment (i.e. urban, rural, suburban) into account
    10072018 -- added function 'selectenv_gridded_idailyCTM'
    20072018 -- function 'open_inst6_3d_ana_Np' added to open MERRA-2 6-hourly
                meteorology
    25072018 -- functions added to open GMI profile information
"""
# # # # # # # # # # # # # 
def open_gmi_singyear(case, year, sampling_months, sampling_hours):
    """function opens hourly GMI trace gas concentrations (O3, CO, NO, and 
    NO2) for the specified model case for a single year. Return concentrations
    represent hourly model output at specified months and hours. Note: 
    hours contained in variable 'sampling_hours' represent hours in coordinated
    universal time. 
    
    Parameters
    ----------    
    case : str
        Hindcast family (i.e. MR2, MR2-CCMI, FFIgac2, etc.)
    year : int
        Year of interest
    sampling_months : list 
        Months of interest
    sampling_hours : list 
        Hours during which trace gas concentrations are returned

    Returns
    ----------      
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
    from datetime import datetime, timedelta
    from netCDF4 import Dataset    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # open GMI hourly column O3 concentrations for year of interest 
    infile = Dataset(pollutants_constants.PATH_GMI + '%s/' %(case)
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
    # loop through days in year of interest with an hourly
    # timestep; for Replay simulations (i.e. case 'M2G_c90') and perturbed 
    # emissions simulations (i.e. case 'EGU_T'), only summer (JJA) model 
    # output was downloaded, created time dimension must be commensurate with 
    # case  
    if case == 'M2G_c90' or case == 'EGU_T':
        for result in perdelta(datetime(year, 6, 1, 0), 
                               datetime(year, 9, 1, 0), 
                               timedelta(hours = 1)):
            times_ty.append(result)     
    else: 
        for result in perdelta(datetime(year, 1, 1, 0), 
                               datetime(year + 1, 1, 1, 0), 
                               timedelta(hours = 1)):
            times_ty.append(result)          
    # find indices/positions of model's time dimension in months/hours 
    # of interest; poi is an array of indicies where hours are within 
    # variable 'sampling_hours' and months within variable 'sampling_months'
    month_idxs = [day.month for day in times_ty]
    month_idxs = np.array(month_idxs)
    hour_idxs = [day.hour for day in times_ty]
    hour_idxs = np.array(hour_idxs)
    poi1 = np.in1d(month_idxs, sampling_months).reshape(month_idxs.shape)
    poi2 = np.in1d(hour_idxs, sampling_hours).reshape(hour_idxs.shape)
    poi = np.where((poi1 == True) & (poi2 == True))[0]
    times_ty = np.array(times_ty)[poi]
    # index trace gas concentrations within 'sampling_hours' and 
    # 'sampling_months'
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
    print('CTM data for %s loaded!' %year)        
    return o3, co, no, no2, gmi_sites, gmi_lat, gmi_lon, times_ty
# # # # # # # # # # # # #
def open_castnet_singyear(year, sampling_months, sampling_hours):
    """function opens yearly CASTNET observations and finds observations taken
    at the specified hours during specified months. 
    
    Parameters
    ----------    
    year : int
        Year of interest
    sampling_months : list 
        Months of interest
    sampling_hours : list 
        Hours during which trace gas concentrations are returned

    Returns
    ----------        
    castnet : pandas.core.frame.DataFrame
        Hourly ozone concentrations at CASTNET stations for the specified 
        year, months, and hours, [no. obs, 5]    
    """
    import numpy as np
    import pytz
    import pandas as pd
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
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
    # syntax: timezone after .tz_localize is timezone that the observations
    # are in, timezone after .tz_convert is the timezone to convert to; 
    # options for timezones can be found by pytz.all_timezones. Code is set 
    # up to handle any UTC offset (i.e. 'GMT+4', 'GMT+5', etc.)
    timezone = 'GMT+4'
    timezone = pytz.timezone('Etc/%s' %timezone)
    castnet['DATE_TIME'] = castnet['DATE_TIME'].dt.values.tz_localize(timezone).tz_convert(pytz.utc)   
    # strip off trailing +00:00 from index
    # sample CASTNET observations only during hours in variable 
    # 'sampling_hours'; n.b. values in 'sampling_hours' are in UTC
    castnet = castnet.loc[castnet['DATE_TIME'].dt.hour.isin(sampling_hours)]
    # consider only ozone observations taken during months contained in 
    # 'sampling_months'
    castnet = castnet.loc[castnet['DATE_TIME'].dt.month.isin(sampling_months)]
    # remove the last three letters of SITE_ID column for direct comparison 
    # with Cooper/GMI station names 
    castnet['SITE_ID'] = castnet['SITE_ID'].astype(str).str[:-3].astype(np.str)
    castnet['DATE_TIME'] = castnet['DATE_TIME'].astype(str).str[:-6]    
    print('CASTNet data for %s loaded!' %year)
    return castnet
# # # # # # # # # # # # #
def open_castnetsiting():
    """function opens CASTNet site information file

    Parameters
    ----------    
    None
    
    Returns
    ----------        
    csites : pandas.core.frame.DataFrame
        CASTNet site information, [147, 19]
    """
    import pandas as pd
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants    
    # load CASTNet siting info
    dirpath = pollutants_constants.PATH_CASTNET
    csites = pd.read_csv((dirpath + 'site.csv'), header = 0)
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
    csites = pd.DataFrame(csites, columns = cols)
    # remove trailing numbers from SITE_ID
    csites['SITE_ID'] = csites['SITE_ID'].map(lambda x: str(x)[:3])      
    return csites    
# # # # # # # # # # # # # 
def commensurate_castnet_gmi(castnet_sites_fr, case, years, 
                             sampling_months, sampling_hours):
    """function opens CASTNet and GMI O3 concentrations during the years, hours
    and months of interest and finds colocated (or nearly colocated) 
    concentrations using the locations of CASTNet sites contained in 
    'castnet_sites_fr'. n.b. function produces a daily average of hourly 
    values assuming that len(sampling_hours) = 6. Change code if this isn't the
    the case. 

    Parameters
    ----------    
    castnet_sites_fr : list
        CASTNET site names in focus region
    case : str
        Hindcast family (i.e. MR2, MR2-CCMI, FFIgac2, etc.)
    years : list
        Years of interest
    sampling_months : list 
        Months of interest
    sampling_hours : list 
        Hours during which trace gas concentrations are returned

    Returns
    ----------      
    comm_castnet : numpy.ndarray
        Existing CASTNet O3 observations at each station contained in variable 
        'castnet_sites_fr' for years in measuring period, units of ppbv, 
        [years in measuring period, stations in 'castnet_sites_fr', 
        days in months in 'sampling_months']
    comm_gmi_o3 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per billion volume, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months']
    comm_gmi_no : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per billion volume, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months']
    comm_gmi_no2 : numpy.ndarray
        GMI CTM NO2 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per billion volume, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months']
    comm_gmi_co : numpy.ndarray     
        GMI CTM CO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per million volume, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months']    
    gmi_sites_fr : list 
        GMI sites colocated (or nearly colocated) with CASTNet sites in 
        'castnet_sites_fr'
    """
    import numpy as np
    import pandas as pd
    from calendar import monthrange
    import warnings    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # GMI/CASTNET site information Table S1 in Strode et al. [2015]
    gmi_castnet_matchup = np.genfromtxt(pollutants_constants.PATH_GMISTATIONINFO + 
        'gmi_castnet_matchup.csv', dtype = 'str', delimiter = ',', 
        skip_header = 1)
    # columns for 'gmi_castnet_matchup' are 
    # [:, 0] - GMI station name (for GMI file naming convention)
    # [:, 1] - CASTNET station name (used in Cooper et al. [2012]?)
    # [:, 2] - TIMEZONE (string, e.g. EST, MST) 
    # [:, 3] - SITENUM 
    # [:, 4] - GMIaltbox (where to sample idaily, etc GMI grid wrt pressure)
    gmi_sites = list(gmi_castnet_matchup[:, 0])     
    castnet_sites = list(gmi_castnet_matchup[:, 1])
    # find number of hours in variable 'sampling_months', n.b. this assumes that
    # the number of hours in 'sampling_months' is constant from year to year,
    # so it wouldn't work for a leap year 
    sos = []
    sos += [(monthrange(years[0], x)[1]) for x in sampling_months]    
    sos = np.sum(sos)
    # create empty arrays to be filled with commensurate CASTNet and GMI O3 
    # concentrations for multi-year period
    comm_castnet = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_castnet[:] = np.nan
    comm_gmi_o3 = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_gmi_o3[:] = np.nan
    comm_gmi_no = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_gmi_no[:] = np.nan
    comm_gmi_no2 = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_gmi_no2[:] = np.nan
    comm_gmi_co = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_gmi_co[:] = np.nan    
    # loop through years in measuring period
    for year, counter1 in zip(years, np.arange(0, len(years), 1)): 
        # load GMI CTM output
        o3, co, no, no2, gmi_sites_i, gmi_lat, gmi_lon, times_ty = \
        open_gmi_singyear(case, year, sampling_months, sampling_hours)
        # times in year of interest with hourly timestep, retrieve only 
        # hours in variable 'sampling_hours' in ET
        dates = pd.date_range('%s-01-%s' %(sampling_months[0], year), 
            '%s-01-%s' %(sampling_months[-1] + 1, year), freq = '1H')[:-1]                
        dates = dates[np.where(dates.hour.isin([x for x in sampling_hours]))[0]]
        # ensure input data matches reference table
        if set(gmi_sites) == set(gmi_sites_i):
            del gmi_sites_i
        else: 
            'Input file sites do not match Strode list!'
        # load CASTNet output
        castnet = open_castnet_singyear(year, sampling_months, sampling_hours)
        # find GMI sites commensurate with CASTNet sites
        gmi_sites_fr = []
        for castnet_site in castnet_sites_fr:
            idx = np.where(np.array(castnet_sites) == castnet_site)
            gmi_sites_fr.append(np.array(gmi_sites)[idx][0]) 
        # loop through GMI and CASTNet sites and a counter in focus region 
        for gmi_site, castnet_site, counter2 in zip(gmi_sites_fr, castnet_sites_fr, 
                                                   np.arange(0, len(castnet_sites_fr), 1)):
            sitepos = np.where(np.array(gmi_sites) == gmi_site)[0]
            # find the average altitude (GMI grid box)
            gmi_sitealt = gmi_castnet_matchup[sitepos, 4]   
            gmi_sitealt = gmi_sitealt.tolist()
            gmi_sitealt = [int(x) for x in gmi_sitealt]
            gmi_sitealt = int(np.mean(gmi_sitealt))   
            # find GMI O3 and other trace gases
            gmi_o3 = o3[sitepos[0], :, gmi_sitealt] 
            gmi_no = no[sitepos[0], :, gmi_sitealt] 
            gmi_no2 = no2[sitepos[0], :, gmi_sitealt] 
            gmi_co = co[sitepos[0], :, gmi_sitealt] 
            # find CASTNet observations
            castnet_atsite = castnet.loc[castnet['SITE_ID'].isin([castnet_site])]
            castnet_atsite = castnet_atsite.sort_values(by = 'DATE_TIME')
            # if CASTNet site has observations for the month(s) and year of 
            # interest, fetch these observations
            if castnet_atsite.shape[0] > 0: 
                # n.b. 2 May 2018: discovered that some CASTNet sites have 
                # more than one entry for each hour in the summer time, if this 
                # is the case, average over hourly values for each duplicate
                # hourly observation (i.e. CASTNet site 'ASH' has 1104 
                # observations for JJA 2000 and it should have 552 
                # observations)
                if castnet_atsite.shape[0] == 1104: 
                    castnet_atsite = castnet_atsite.groupby(['SITE_ID','DATE_TIME']).mean()
                # add missing values                 
                else: 
                    castnet_atsite.index = pd.DatetimeIndex(castnet_atsite['DATE_TIME'])
                    castnet_atsite = castnet_atsite.reindex(dates, 
                        fill_value = np.nan)
                # average every 6 hours (i.e. 11 - 16 local time) for a 
                # daily average and supress warnings
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message="Mean of empty slice")
                    comm_castnet[counter1, counter2] = np.nanmean(np.array(
                        castnet_atsite['OZONE'].values).reshape(-1, 6), axis = 1)
                    comm_gmi_o3[counter1, counter2] = np.nanmean(np.array(
                        gmi_o3 * 1e9).reshape(-1, 6), axis = 1)      
                    comm_gmi_no[counter1, counter2] = np.nanmean(np.array(
                        gmi_no * 1e9).reshape(-1, 6), axis = 1)   
                    comm_gmi_no2[counter1, counter2] = np.nanmean(np.array(
                        gmi_no2 * 1e9).reshape(-1, 6), axis = 1)   
                    comm_gmi_co[counter1, counter2] = np.nanmean(np.array(
                        gmi_co * 1e6).reshape(-1, 6), axis = 1)                   
    return comm_castnet, comm_gmi_o3, comm_gmi_no, comm_gmi_no2, comm_gmi_co, gmi_sites_fr
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
    return t2m_all, lat, lon
# # # # # # # # # # # # #   
def commensurate_t2m(comm_castnet, castnet_sites_fr, years, sampling_months):
    """function loads MERRA 2 meter temperatures and identifes MERRA grid cells 
    co-located with CASTNet sites which have data for the years of interest
    in focus region.
    
    Parameters
    ----------    
    comm_castnet : numpy.ndarray
        Existing CASTNet O3 observations at each station contained in variable 
        'castnet_sites_fr' for years in measuring period, units of ppbv, 
        [years in measuring period, stations in 'castnet_sites_fr', 
        days in months in 'sampling_months']
    castnet_sites_fr : list
        CASTNET site names in focus region
    years : list
        Years of interest
    sampling_months : list 
        Months of interest

    Returns
    ----------
    comm_t2m : numpy.ndarray
        MERRA 2 meter temperatures co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of K, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']
    merra_lats_fr : list  
        Longitudes of MERRA grid cells which are co-located (or nearly co-
        located with CASTNet observations during years of interest
    merra_lons_fr : list
        Latitudes of MERRA grid cells which are co-located (or nearly co-
        located with CASTNet observations during years of interest        
    """
    import numpy as np
    from calendar import monthrange
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    # load MERRA 2 meter temperature
    t2m_all, lat, lon = load_t2m()
    # find summer days in measuring period which spans T2m record in 
    # climatology file
    daysmy = find_summer_days()
    # find lat/lon of CASTNet sites in focus region
    csites = open_castnetsiting()
    # create empty arrays to be filled with temperatures commensurate 
    # to CASTNet observations    
    sos = []
    sos += [(monthrange(years[0], x)[1]) for x in sampling_months]    
    sos = np.sum(sos)
    comm_t2m = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_t2m[:] = np.nan
    # create lists to be filled with MERRA grid cell locations
    merra_lats_fr, merra_lons_fr = [], []        
    # loop through years in measuring period 
    for year, counter1 in zip(years, np.arange(0, len(years), 1)): 
        # in a given year, loop through CASTNet sites in focus region and 
        # determine if data exists
        for site, castnet_data, counter2 in zip(castnet_sites_fr, 
                                                comm_castnet[counter1], 
                                                np.arange(0, len(castnet_sites_fr), 1)):
            # if data exists, then find MERRA T2m commensurate with CASTNet site
            if np.where(np.isnan(castnet_data) == True)[0].shape[0] != sos:
                # find longitude/latitude of CASTNET station 
                csites_atsite = csites.loc[csites['SITE_ID'].isin([site])]
                lon_atsite = csites_atsite['LONGITUDE'].values[0]
                lat_atsite = csites_atsite['LATITUDE'].values[0]
                # find closest MERRA grid cell 
                lat_idx = geo_idx(lat_atsite, lat)
                lon_idx = geo_idx(lon_atsite, lon)
                t2m_atsite = t2m_all[np.where(daysmy.year.isin([year]) == True), 
                                     lat_idx, lon_idx]
                comm_t2m[counter1, counter2] = t2m_atsite
                merra_lats_fr.append(lat[lat_idx]) 
                merra_lons_fr.append(lon[lon_idx])
        print('T2m data for %s loaded!' %year)
    return comm_t2m, merra_lats_fr, merra_lons_fr
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
def load_aqshourly(years): 
    """function loads hourly observations of AQS NO2, CO, and O3
    
    Parameters
    ---------- 
    years : list
        Years of interest 
        
    Returns
    ----------    
    COdf : pandas.core.frame.DataFrame
        Raw AQS hourly summary CO observations for U.S. for year(s) of 
        interest, (no. obs, 24)
    NO2df : pandas.core.frame.DataFrame
        Raw AQS hourly summary NO2 observations for U.S. for year(s) of 
        interest, (no. obs, 24)
    O3df : pandas.core.frame.DataFrame
        Raw AQS hourly summary O3 observations for U.S. for year(s) of 
        interest, (no. obs, 24)
    """
    import numpy as np
    import sys, os
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    #cols = ["State Code", "County Code", "Site Num", "Parameter Code", 
    #        "POC", "Latitude", "Longitude", "Datum", "Parameter Name",
    #        "Date Local", "Time Local", "Date GMT", "Time GMT", 
    #        "Sample Measurement", "Units of Measure", "MDL", "Uncertainty",
    #        "Qualifier", "Method Type", "Method Code", "Method Name",
    #        "State Name", "County Name", "Date of Last Change"]
    cols = ["Latitude", "Longitude", "Date GMT", "Time GMT", 
            "Sample Measurement"]
    # specificy data types for columns to avoid DtypeWarning being raised
    #dtype = {'State Code' : np.str, 'County Code' : np.str, 'Site Num' : np.str, 
    #         'Parameter Code' : np.str, 'POC' : np.int64, 'Latitude' : np.float64, 
    #         'Longitude' : np.float64, 'Datum' : np.str, 'Parameter Name' : np.str,
    #         'Date Local' : np.str, 'Time Local' : np.str, 'Date GMT' : np.str, 
    #         'Time GMT' : np.str, 'Sample Measurement' : np.float64,
    #         'Units of Measure' : np.str,  'MDL' : np.float64, 
    #         'Uncertainty' : np.str, 'Qualifier' : np.str, 
    #         'Method Type' : np.str, 'Method Code' : np.str, 
    #         'Method Name' : np.str,'State Name' : np.str, 
    #         "County Name" : np.str, 'Date of Last Change' : np.str}  
    dtype = {'Latitude' : np.float64, 'Longitude' : np.float64, 
             'Date GMT' : np.str, 'Time GMT' : np.str, 
             'Sample Measurement' : np.float64}      
    # fetch file names for measuring period
    COfiles, NO2files, O3files = [], [], []
    for year in years:
        COfiles.append(os.path.join(pollutants_constants.PATH_AQS_CO, 
                                    'hourly_42101_%s_reduced.csv' %year))
        NO2files.append(os.path.join(pollutants_constants.PATH_AQS_NO2, 
                                     'hourly_42602_%s_reduced.csv' %year))
        O3files.append(os.path.join(pollutants_constants.PATH_AQS_O3, 
                                     'hourly_44201_%s_reduced.csv' %year))
    # read multiple CSV files (yearly) into Pandas dataframe for hourly CO, 
    # NO2, and O3        
    COdf = get_merged_csv(COfiles, dtype = dtype, index_col = None, 
                          usecols = cols)
    NO2df = get_merged_csv(NO2files, dtype = dtype, index_col = None, 
                           usecols = cols)  
    O3df = get_merged_csv(O3files, dtype = dtype, index_col = None, 
                          usecols = cols) 
    return COdf, NO2df, O3df 
# # # # # # # # # # # # # 
def commensurate_aqstracegas(castnet_sites_fr, years, 
                             sampling_months, sampling_hours): 
    """function loads hourly AQS CO, NO2, and O3 observations and identifies 
    observations near CASTNet sites which have data for the years of interest
    in the focus region. In this case the variable 'searchrad' constitutes 
    a bounding box in which to search for AQS sites near CASTNet sites. 
    Returned values represent daily average values averaged over the hours 
    in variable 'sampling_hours'
    
    Parameters
    ----------    
    castnet_sites_fr : list
        CASTNET site names in focus region
    years : list
        Years of interest
    sampling_months : list 
        Months of interest
    sampling_hours : list
        Hours of interest
        
    Returns
    ----------
    comm_co : numpy.ndarray
        AQS CO observations within a bounding box defined by variable 
        'searchrad' of the location of CASTNet station averaged over the hours
        in variable 'sampling_hours'. For all AQS CO observations within 
        bounding box, measurements are averaged, units of parts per million, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months']
    comm_no2 : numpy.ndarray
        AQS NO2 observations within a bounding box defined by variable 
        'searchrad' of the location of CASTNet station averaged over the hours
        in variable 'sampling_hours'. For all AQS NO2 observations within 
        bounding box, measurements are averaged, units of parts per billion, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months']
    comm_o3 : numpy.ndarray
        AQS O3 observations within a bounding box defined by variable 
        'searchrad' of the location of CASTNet station averaged over the hours
        in variable 'sampling_hours'. For all AQS O3 observations within 
        bounding box, measurements are averaged, units of parts per billion, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months']        
    aqs_co_coords : list
        Coordinates of unique AQS stations measuring CO within bounding boxes 
        defined by CASTNet stations
    aqs_no2_coords : list
        Coordinates of unique AQS stations measuring NO2 within bounding boxes 
        defined by CASTNet stations
    aqs_o3_coords : list
        Coordinates of unique AQS stations measuring O3 within bounding boxes 
        defined by CASTNet stations     
    """
    import numpy as np
    import pandas as pd
    from calendar import monthrange
    # open AQS NO2/CO observations 
    CO, NO2, O3 = load_aqshourly(years)
    # find lat/lon of CASTNet sites in focus region
    csites = open_castnetsiting()
    # total number of days in each year
    sos = []
    sos += [(monthrange(years[0], x)[1]) for x in sampling_months]    
    sos = np.sum(sos)
    # create empty arrays to be filled with CO, NO2 observations commensurate 
    # to CASTNet observations    
    comm_co = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_co[:] = np.nan
    comm_no2 = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_no2[:] = np.nan
    comm_o3 = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_o3[:] = np.nan    
    # open CASTNet observations for particular model configuration (in this 
    # case 'HindcastMR2') to see which stations have observations for the 
    # measuring period
    comm_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
    commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2', years, 
                             sampling_months, sampling_hours)    
    # convert time
    CO['Date GMT'] = pd.to_datetime(CO['Date GMT'])
    NO2['Date GMT'] = pd.to_datetime(NO2['Date GMT'])
    O3['Date GMT'] = pd.to_datetime(O3['Date GMT'])       
    # define search radius, all AQS stations within this distance of CASTNet
    # site will be fetched
    searchrad = 1.0
    # list creation, will be filled with lat/lon coordinates of AQS stations
    aqs_co_coords = []
    aqs_no2_coords = []
    aqs_o3_coords = []    
    # loop through years in measuring period 
    for year, counter1 in zip(years, np.arange(0, len(years), 1)): 
        # find AQS NO2, CO, and O3 observations for year and in variables 
        # 'sampling_months' and 'sampling_hours'
        NO2_year = NO2.loc[NO2['Date GMT'].dt.year.isin([year])]
        NO2_year = NO2_year.loc[NO2_year['Date GMT'].dt.month.isin(sampling_months)]  
        NO2_year = NO2_year.loc[NO2_year['Time GMT'].isin(
                ['%d:00' %(x) for x in sampling_hours])] 
        CO_year = CO.loc[CO['Date GMT'].dt.year.isin([year])]
        CO_year = CO_year.loc[CO_year['Date GMT'].dt.month.isin(sampling_months)]  
        CO_year = CO_year.loc[CO_year['Time GMT'].isin(
                ['%d:00' %(x) for x in sampling_hours])] 
        O3_year = O3.loc[O3['Date GMT'].dt.year.isin([year])]
        O3_year = O3_year.loc[O3_year['Date GMT'].dt.month.isin(sampling_months)]  
        O3_year = O3_year.loc[O3_year['Time GMT'].isin(
                ['%d:00' %(x) for x in sampling_hours])]             
        # in a given year, loop through CASTNet sites in focus region and
        # determine if data exists  
        for site, castnet_data, counter2 in zip(castnet_sites_fr, 
                                                comm_castnet[counter1], 
                                                np.arange(0, len(castnet_sites_fr), 1)):
            # if data exists, the find AQS with the distance from the CASTNet site
            # defined by the search radius
            if np.where(np.isnan(castnet_data) == True)[0].shape[0] != sos:
                # find longitude/latitude of CASTNET station 
                csites_atsite = csites.loc[csites['SITE_ID'].isin([site])]
                lon_atsite = csites_atsite['LONGITUDE'].values[0]
                lat_atsite = csites_atsite['LATITUDE'].values[0]
                # define bounding box            
                top = lat_atsite + searchrad
                left = lon_atsite - searchrad
                bottom = lat_atsite - searchrad
                right = lon_atsite + searchrad
                # NO2/CO observations in bounding box
                CO_atsite = CO_year[(CO_year['Latitude'] > bottom) & 
                                    (CO_year['Latitude'] <= top) &
                                    (CO_year['Longitude'] > left) &
                                    (CO_year['Longitude'] <= right)]   
                NO2_atsite = NO2_year[(NO2_year['Latitude'] > bottom) & 
                                      (NO2_year['Latitude'] <= top) &
                                      (NO2_year['Longitude'] > left) &
                                      (NO2_year['Longitude'] <= right)]
                O3_atsite = O3_year[(O3_year['Latitude'] > bottom) & 
                                    (O3_year['Latitude'] <= top) &
                                    (O3_year['Longitude'] > left) &
                                    (O3_year['Longitude'] <= right)]                
                # save off unique latitudes and longitudes 
                if CO_atsite.empty != True:
                    aqs_co_coords.append(CO_atsite[['Latitude', 
                        'Longitude']].drop_duplicates().values)
                if NO2_atsite.empty != True:
                    aqs_no2_coords.append(NO2_atsite[['Latitude', 
                        'Longitude']].drop_duplicates().values)
                if O3_atsite.empty != True:
                    aqs_o3_coords.append(O3_atsite[['Latitude', 
                        'Longitude']].drop_duplicates().values)                    
                # find hours in sampling_hours 
                CO_atsite = CO_atsite.loc[CO_atsite['Time GMT'].isin(
                        ['%s:00' %x for x in sampling_hours])] 
                NO2_atsite = NO2_atsite.loc[NO2_atsite['Time GMT'].isin(
                        ['%s:00' %x for x in sampling_hours])]                     
                O3_atsite = O3_atsite.loc[O3_atsite['Time GMT'].isin(
                        ['%s:00' %x for x in sampling_hours])]                                         
                # group by date and average  
                CO_atsite = CO_atsite.groupby(['Date GMT']).mean()
                CO_atsite = CO_atsite['Sample Measurement']
                NO2_atsite = NO2_atsite.groupby(['Date GMT']).mean()
                NO2_atsite = NO2_atsite['Sample Measurement']
                O3_atsite = O3_atsite.groupby(['Date GMT']).mean()
                O3_atsite = O3_atsite['Sample Measurement']                
                # dates in year/months of interest
                if sampling_months == [6, 7, 8]:
                    idx = pd.date_range('%s-01-%s' %(sampling_months[0], year),                                         
                                        '%s-31-%s' %(sampling_months[-1], year))
                if sampling_months == [8, 9]:
                    idx = pd.date_range('%s-01-%s' %(sampling_months[0], year), 
                                        '%s-30-%s' %(sampling_months[-1], year))                    
                # add missing dates to pandas dataframe from 
                # https://stackoverflow.com/questions/19324453/add-missing-dates-to-pandas-dataframe
                CO_atsite.index = pd.DatetimeIndex(CO_atsite.index)
                CO_atsite = CO_atsite.reindex(idx, fill_value = np.nan).values
                NO2_atsite.index = pd.DatetimeIndex(NO2_atsite.index)
                NO2_atsite = NO2_atsite.reindex(idx, fill_value = np.nan).values   
                O3_atsite.index = pd.DatetimeIndex(O3_atsite.index)
                O3_atsite = O3_atsite.reindex(idx, fill_value = np.nan).values                   
                # add to array
                comm_co[counter1, counter2] = CO_atsite
                comm_no2[counter1, counter2] = NO2_atsite
                comm_o3[counter1, counter2] = O3_atsite
        print('AQS CO/NO2 data for %s loaded!' %year)   
    return (comm_co, comm_no2, comm_o3, aqs_co_coords, aqs_no2_coords, 
            aqs_o3_coords)
# # # # # # # # # # # # #
def commensurate_aqstracegas_gridded(df, gmi, times, lat, lon, environment):
    """function uses CTM grid and, for each day with 12Z CTM output, finds 
    all available AQS observations in each grid cell if all siting environments
    are needed. If only rural, suburban, and/or urban observations are needed, 
    then function retrieves just observations in these environments. If no AQS 
    stations are located within the bounds of a CTM grid cell, a nan value is 
    returned. If > 1 AQS station exists in a grid cell, the daily values at 
    these stations are averaged. 
    
    Parameters
    ----------        
    df : pandas.core.frame.DataFrame
        Raw AQS hourly summary trace gas observations for U.S. for years of 
        interest, [no. obs, 5]
    gmi : numpy.ndarray
        Gridded CTM output, units of volume mixing ratio, [time, pressure, lat, 
        lon]    
    times : numpy.ndarray
        datetime.datetime objects corresponding to 12Z daily CTM output, 
        [time,] 
    lat : numpy.ndarray
        Latitude coordinates of focus region, units of degrees north, [lat,]
    lon : numpy.ndarray
        Longitude coordinates of focus region, units of degrees east, [lon,]
    environment : list
        Environment (i.e., rural, urban, and/or suburban) in which AQS 
        observations are retrieved; n.b. if all observations are wanted, 
        pass [] (faster) or ['URBAN AND CENTER CITY', 'RURAL', 'SUBURBAN']
        (slower). Otherwise, pass a list of individual environments (i.e., 
        ['URBAN'])
        
    Returns
    ----------
    station_coordinates : list
        List of arrays containing the unique locations (latitude/longitude 
        coordinates) of AQS stations in CTM grid cells
    colocated_tg : numpy.ndarray
        Daily 12Z AQS trace gas measurements for trace gas of interest (i.e., 
        NO2, CO, or O3) within CTM grid cells and environment of interest. If 
        no AQS stations are located within the bounds of a CTM grid cell, a nan 
        value is returned. If > 1 AQS station exists in a grid cell, the daily 
        values at these stations are averaged, [time, lat, lon]
    """
    import numpy as np
    import pandas as pd
    # read AQS site file containing information about the siting of AQS
    # sites
    datapath = '/Users/ghkerr/phd/aqs_station_siting/data/'
    aqs_sites = pd.read_csv((datapath + 'aqs_sites.csv'), header = 0)
    site_cols = ['State Code', 'County Code',	 'Site Number', 'Latitude',
                 'Longitude', 'Datum', 'Elevation', 'Land Use',	
                 'Location Setting', 'Site Established Date', 
                 'Site Closed Date', 'Met Site State Code', 
                 'Met Site County Code', 'Met Site Site Number', 
                 'Met Site Type', 'Met Site Distance', 'Met Site Direction', 
                 'GMT Offset', 'Owning Agency', 'Local Site Name', 'Address', 
                 'Zip Code', 'State Name', 'County Name', 'City Name', 
                 'CBSA Name', 'Tribe Name', 'Extraction Date']
    aqs_sites = pd.DataFrame(aqs_sites, columns = site_cols)     
    # list to be filled with latitude and longitude coordinates of AQS 
    # stations
    station_coordinates = []
    new = []
    # select only 12Z observations
    df = df.loc[df['Time GMT'].isin(['12:00'])]
    # convert input parameter 'times' into 'YYYY-MM-DD' format 
    times_ymd = [pd.to_datetime(d).strftime('%Y-%m-%d') for d in times]
    # select 12Z observations in JJA
    df = df.loc[df['Date GMT'].isin(times_ymd)]
    # create grid with same resolution as CTM to be filled with AQS trace 
    # gas measurements, grid is 3D (not 4D)
    colocated_tg = np.empty((gmi.shape[0], gmi.shape[2], gmi.shape[3]))
    colocated_tg[:] = np.nan
    # CTM resolution 
    latres = np.diff(lat).mean()
    lonres = np.diff(lon).mean()
    # loop through GMI latitude and longitude coordinatesion
    for i, ilat in enumerate(lat):
        for j, jlon in enumerate(lon):
            # convert longitude from (0-360) to (-180 to 180)
            jlon = np.mod(jlon - 180.0, 360.0) - 180.0
            # find AQS trace gas measurements in grid cell 
            df_ingrid = df.loc[(df['Latitude'] > ilat - latres/2.) & 
                               (df['Latitude'] <= ilat + latres/2.) &
                               (df['Longitude'] > jlon - lonres/2.) &
                               (df['Longitude'] <= jlon + lonres/2.)]
            if df_ingrid.shape[0] > 0: 
                # if AQS observations from all stations (rural, urban, and sub-
                # urban) are needed, following lines are executed
                if environment == []:
                    # save off unique latitudes and longitudes
                    station_coordinates.append(df_ingrid[
                            ['Latitude', 'Longitude']].drop_duplicates().values)                    
                    # in the case that >1 AQS site exists in grid cell, average 
                    # over all sites in cell 
                    df_ingrid = df_ingrid.groupby(['Date GMT']).mean()
                    # add missing values                 
                    df_ingrid = df_ingrid.reindex(times_ymd, fill_value = np.nan)
                    # add daily grid cell averages to grid
                    colocated_tg[:, i, j] = df_ingrid['Sample Measurement'].values
                # if AQS observations from a particular environment (i.e., 
                # rural, urban, or suburban) are needed
                else: 
                    # list of latitude, longitudes in environment of interest
                    lats_in_env, lons_in_env = [], []
                    for (slat, slon) in df_ingrid[['Latitude', 
                        'Longitude']].drop_duplicates().values:
                        # find index in AQS site record corresponding to 
                        # station's observations
                        siting_lat = np.where(np.around(
                                aqs_sites['Latitude'].values, 4) == np.around(
                                        slat, 4))
                        siting_lon = np.where(np.around(
                                aqs_sites['Longitude'].values, 4) == np.around(
                                        slon, 4))
                        siting = np.intersect1d(siting_lat, siting_lon)
                        siting = aqs_sites.iloc[siting[0]]['Location Setting']                    
                        # only save off observations from environment of 
                        # interest
                        if siting in environment: 
                            lats_in_env.append(slat)
                            lons_in_env.append(slon)
                    # find observations taken at latitudes/longitudes in 
                    # environment
                    df_ingrid = df_ingrid.loc[df_ingrid['Latitude'].isin(
                            lats_in_env) & df_ingrid['Longitude'].isin(
                            lons_in_env)]
                    # same as above 
                    station_coordinates.append(df_ingrid[
                            ['Latitude', 'Longitude']].drop_duplicates().values)                        
                    df_ingrid = df_ingrid.groupby(['Date GMT']).mean()
                    df_ingrid = df_ingrid.reindex(times_ymd, fill_value = np.nan)
                    colocated_tg[:, i, j] = df_ingrid['Sample Measurement'].values
    return station_coordinates, colocated_tg
# # # # # # # # # # # # #    
def commensurate_aqstracegas_diurnal(comm_castnet, castnet_sites_fr, years, 
                                     sampling_months):
    """function loads hourly AQS CO and NO2 observations and identifies 
    observations near CASTNet sites which have data for the years of interest
    in the focus region. In this case the variable 'searchrad' constitutes 
    a bounding box in which to search for AQS sites near CASTNet sites. 
    All hourly observations are retrieved and output as to better understand
    the diurnal cycle of NO2 and CO. 
    
    Parameters
    ----------    
    comm_castnet : numpy.ndarray
        Existing CASTNet O3 observations at each station contained in variable 
        'castnet_sites_fr' for years in measuring period, units of ppbv, 
        [years in measuring period, stations in 'castnet_sites_fr', 
        days in months in 'sampling_months']
    castnet_sites_fr : list
        CASTNET site names in focus region
    years : list
        Years of interest
    sampling_months : list 
        Months of interest
        
    Returns
    ----------
    comm_co : numpy.ndarray
        Hourly AQS CO observations within a bounding box defined by variable 
        'searchrad' of the location of CASTNet station averaged over the hours
        in variable 'sampling_hours'. For all AQS CO observations within 
        bounding box, measurements are averaged, units of parts per million, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months', 24]
    comm_no2 : numpy.ndarray
        Hourly AQS NO2 observations within a bounding box defined by variable 
        'searchrad' of the location of CASTNet station averaged over the hours
        in variable 'sampling_hours'. For all AQS NO2 observations within 
        bounding box, measurements are averaged, units of parts per billion, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months', 24]   
    comm_o3 : numpy.ndarray
        Hourly AQS NO2 observations within a bounding box defined by variable 
        'searchrad' of the location of CASTNet station averaged over the hours
        in variable 'sampling_hours'. For all AQS O3 observations within 
        bounding box, measurements are averaged, units of parts per billion, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months', 24]         
    aqs_co_coords : list
        Coordinates of unique AQS stations measuring CO within bounding boxes 
        defined by CASTNet stations
    aqs_no2_coords : list
        Coordinates of unique AQS stations measuring NO2 within bounding boxes 
        defined by CASTNet stations
    aqs_o3_coords : list
        Coordinates of unique AQS stations measuring O3 within bounding boxes 
        defined by CASTNet stations      
    """
    import numpy as np
    import pandas as pd
    from calendar import monthrange
    from datetime import datetime    
    # open AQS NO2/CO observations 
    CO, NO2, O3 = load_aqshourly(years)
    # find lat/lon of CASTNet sites in focus region
    csites = open_castnetsiting()
    # total number of days in each year
    sos = []
    sos += [(monthrange(years[0], x)[1]) for x in sampling_months]    
    sos = np.sum(sos)
    # create empty arrays to be filled with CO, NO2 observations commensurate 
    # to CASTNet observations; final dimension (24) is for each hour in day
    comm_co = np.empty([len(years), len(castnet_sites_fr), sos, 24], 
                             dtype = float)
    comm_co[:] = np.nan
    comm_no2 = np.empty([len(years), len(castnet_sites_fr), sos, 24], 
                             dtype = float)
    comm_no2[:] = np.nan
    comm_o3 = np.empty([len(years), len(castnet_sites_fr), sos, 24], 
                             dtype = float)
    comm_o3[:] = np.nan    
    # convert time
    CO['Date GMT'] = pd.to_datetime(CO['Date GMT'])
    NO2['Date GMT'] = pd.to_datetime(NO2['Date GMT'])
    O3['Date GMT'] = pd.to_datetime(O3['Date GMT'])
    # define search radius, all AQS stations within this distance of CASTNet
    # site will be fetched
    searchrad = 1.0
    # list creation; will be filled with unique longitudes and latitudes of 
    # AQS stations nearby CASTNet sites        
    aqs_co_coords = []
    aqs_no2_coords = []
    aqs_o3_coords = []
    # loop through years in measuring period 
    for year, counter1 in zip(years, np.arange(0, len(years), 1)): 
        # find AQS NO2/CO observations for year and in variables 
        # 'sampling_months' and do not restrict hours 
        NO2_year = NO2.loc[NO2['Date GMT'].dt.year.isin([year])]
        NO2_year = NO2_year.loc[NO2_year['Date GMT'].dt.month.isin(sampling_months)]  
        CO_year = CO.loc[CO['Date GMT'].dt.year.isin([year])]
        CO_year = CO_year.loc[CO_year['Date GMT'].dt.month.isin(sampling_months)] 
        O3_year = O3.loc[O3['Date GMT'].dt.year.isin([year])]
        O3_year = O3_year.loc[O3_year['Date GMT'].dt.month.isin(sampling_months)]         
        # dates in year/months of interest
        if sampling_months == [6, 7, 8]:
            idx = pd.date_range('%s-01-%s' %(sampling_months[0], year),                                         
                                '%s-31-%s' %(sampling_months[-1], year))
        if sampling_months == [8, 9]:
            idx = pd.date_range('%s-01-%s' %(sampling_months[0], year), 
                                '%s-30-%s' %(sampling_months[-1], year))   
        # in a given year, loop through CASTNet sites in focus region and
        # determine if data exists  
        for site, castnet_data, counter2 in zip(castnet_sites_fr, 
                                                comm_castnet[counter1], 
                                                np.arange(0, len(castnet_sites_fr), 1)):
            # if data exists, the find AQS with the distance from the CASTNet site
            # defined by the search radius
            if np.where(np.isnan(castnet_data) == True)[0].shape[0] != sos:
                # find longitude/latitude of CASTNET station 
                csites_atsite = csites.loc[csites['SITE_ID'].isin([site])]
                lon_atsite = csites_atsite['LONGITUDE'].values[0]
                lat_atsite = csites_atsite['LATITUDE'].values[0]
                # define bounding box            
                top = lat_atsite + searchrad
                left = lon_atsite - searchrad
                bottom = lat_atsite - searchrad
                right = lon_atsite + searchrad
                # NO2/CO observations in bounding box
                CO_atsite = CO_year[(CO_year['Latitude'] > bottom) & 
                                    (CO_year['Latitude'] <= top) &
                                    (CO_year['Longitude'] > left) &
                                    (CO_year['Longitude'] <= right)]   
                NO2_atsite = NO2_year[(NO2_year['Latitude'] > bottom) & 
                                      (NO2_year['Latitude'] <= top) &
                                      (NO2_year['Longitude'] > left) &
                                      (NO2_year['Longitude'] <= right)] 
                O3_atsite = O3_year[(O3_year['Latitude'] > bottom) & 
                                    (O3_year['Latitude'] <= top) &
                                    (O3_year['Longitude'] > left) &
                                    (O3_year['Longitude'] <= right)]                 
                # find latitudes of AQS stations within bounding box and 
                # thereafter find unique combinations 
                if CO_atsite.empty != True:
                    aqs_co_coords.append(CO_atsite[['Latitude', 
                        'Longitude']].drop_duplicates().values)
                if NO2_atsite.empty != True:
                    aqs_no2_coords.append(NO2_atsite[['Latitude', 
                        'Longitude']].drop_duplicates().values)
                if O3_atsite.empty != True:
                    aqs_o3_coords.append(O3_atsite[['Latitude', 
                        'Longitude']].drop_duplicates().values)   
                # loop through days in variable 'sampling_months' in year of 
                # interest
                for day, counter3 in zip(idx, np.arange(0, len(idx), 1)): 
                    all_times = pd.date_range(day, day + 1, freq = '1H')[:-1]
                    day = datetime.strftime(day, '%Y-%m-%d')
                    # find CO/NO2 on day
                    CO_atsite_onday = CO_atsite.loc[CO_atsite['Date GMT'].isin([day])]
                    NO2_atsite_onday = NO2_atsite.loc[NO2_atsite['Date GMT'].isin([day])]
                    O3_atsite_onday = O3_atsite.loc[O3_atsite['Date GMT'].isin([day])]                    
                    # average over multiple locations (if any)
                    CO_atsite_onday = CO_atsite_onday.groupby(['Time GMT']).mean()
                    CO_atsite_onday = CO_atsite_onday['Sample Measurement']
                    NO2_atsite_onday = NO2_atsite_onday.groupby(['Time GMT']).mean()
                    NO2_atsite_onday = NO2_atsite_onday['Sample Measurement']
                    O3_atsite_onday = O3_atsite_onday.groupby(['Time GMT']).mean()
                    O3_atsite_onday = O3_atsite_onday['Sample Measurement']                    
                    # fill in missing hours with NaNs
                    CO_atsite_onday.index = pd.DatetimeIndex(
                            day + ' ' + CO_atsite_onday.index)
                    CO_atsite_onday = CO_atsite_onday.reindex(all_times, 
                        fill_value = np.nan).values
                    NO2_atsite_onday.index = pd.DatetimeIndex(
                            day + ' ' + NO2_atsite_onday.index)
                    NO2_atsite_onday = NO2_atsite_onday.reindex(all_times, 
                        fill_value = np.nan).values   
                    O3_atsite_onday.index = pd.DatetimeIndex(
                            day + ' ' + O3_atsite_onday.index)
                    O3_atsite_onday = O3_atsite_onday.reindex(all_times, 
                        fill_value = np.nan).values                                                                   
                    # add to array
                    comm_co[counter1, counter2, counter3] = CO_atsite_onday
                    comm_no2[counter1, counter2, counter3] = NO2_atsite_onday
                    comm_o3[counter1, counter2, counter3] = O3_atsite_onday
    return (comm_co, comm_no2, comm_o3, aqs_co_coords, aqs_no2_coords, 
            aqs_o3_coords)
# # # # # # # # # # # # 
def commensurate_castnet_gmi_diurnal(castnet_sites_fr, case, years, 
                                     sampling_months):
    """function opens CASTNet and GMI trace gas concentrations during the years
    and months of interest and finds colocated (or nearly colocated) 
    concentrations using the locations of CASTNet sites contained in 
    'castnet_sites_fr'. The diurnal cycles (i.e. 0 - 23 UTC) are retrieved for 
    both GMI and the CASTNet observations.

    Parameters
    ----------    
    castnet_sites_fr : list
        CASTNET site names in focus region
    case : str
        Hindcast family (i.e. MR2, MR2-CCMI, FFIgac2, etc.)
    years : list
        Years of interest
    sampling_months : list 
        Months of interest

    Returns
    ----------      
    comm_castnet : numpy.ndarray
        Diurnal measurements of existing CASTNet O3 observations at each 
        station contained in variable 'castnet_sites_fr' for years in measuring 
        period, units of ppbv, [years in measuring period, stations in 
        'castnet_sites_fr', days in months in 'sampling_months', 24]    
    comm_gmi_o3 : numpy.ndarray
        Diurnal GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per billion volume,
        [years in measuring period, stations in 'castnet_sites_fr', 
        days in months in 'sampling_months', 24]    
    comm_gmi_no : numpy.ndarray
        Diurnal GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per billion volume,
        [years in measuring period, stations in 'castnet_sites_fr', 
        days in months in 'sampling_months', 24]    
    comm_gmi_no2 : numpy.ndarray
        Diurnal GMI CTM NO2 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per billion volume,
        [years in measuring period, stations in 'castnet_sites_fr', 
        days in months in 'sampling_months', 24]
    comm_gmi_co : numpy.ndarray
        Diurnal GMI CTM CO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per million volume,
        [years in measuring period, stations in 'castnet_sites_fr', 
        days in months in 'sampling_months', 24]    
    """
    import pandas as pd
    import numpy as np
    from calendar import monthrange
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # GMI/CASTNET site information Table S1 in Strode et al. [2015]
    gmi_castnet_matchup = np.genfromtxt(pollutants_constants.PATH_GMISTATIONINFO + 
        'gmi_castnet_matchup.csv', dtype = 'str', delimiter = ',', 
        skip_header = 1)
    # columns for 'gmi_castnet_matchup' are 
    # [:, 0] - GMI station name (for GMI file naming convention)
    # [:, 1] - CASTNET station name (used in Cooper et al. [2012]?)
    # [:, 2] - TIMEZONE (string, e.g. EST, MST) 
    # [:, 3] - SITENUM 
    # [:, 4] - GMIaltbox (where to sample idaily, etc GMI grid wrt pressure)
    gmi_sites = list(gmi_castnet_matchup[:, 0])     
    castnet_sites = list(gmi_castnet_matchup[:, 1])
    # find number of hours in variable 'sampling_months', n.b. this assumes that
    # the number of hours in 'sampling_months' is constant from year to year,
    # so it wouldn't work for a leap year 
    sos = []
    sos += [(monthrange(years[0], x)[1]) for x in sampling_months]    
    sos = np.sum(sos)
    # for retrieval of diurnal cycle
    sampling_hours = [x for x in np.arange(0, 24, 1)]
    # create empty arrays to be filled with commensurate CASTNet and GMI O3, NO2, 
    # and CO concentrations for multi-year period
    comm_castnet = np.empty([len(years), len(castnet_sites_fr), sos, 24], 
                             dtype = float)
    comm_castnet[:] = np.nan
    comm_gmi_o3 = np.empty([len(years), len(castnet_sites_fr), sos, 24], 
                             dtype = float)
    comm_gmi_o3[:] = np.nan
    comm_gmi_no = np.empty([len(years), len(castnet_sites_fr), sos, 24], 
                             dtype = float)
    comm_gmi_no[:] = np.nan
    comm_gmi_no2 = np.empty([len(years), len(castnet_sites_fr), sos, 24], 
                             dtype = float)
    comm_gmi_no2[:] = np.nan
    comm_gmi_co = np.empty([len(years), len(castnet_sites_fr), sos, 24], 
                             dtype = float)
    comm_gmi_co[:] = np.nan
    # loop through years in measuring period
    for year, counter1 in zip(years, np.arange(0, len(years), 1)): 
        # load GMI CTM output and retrieve all hours to understand diurnal cycle
        o3, co, no, no2, gmi_sites_i, gmi_lat, gmi_lon, times_ty = \
        open_gmi_singyear(case, year, sampling_months, sampling_hours)
        # times in year of interest with hourly timestep 
        dates = pd.date_range('%s-01-%s' %(sampling_months[0], year), 
            '%s-01-%s' %(sampling_months[-1] + 1, year), freq = '1H')[:-1]
        # ensure input data matches reference table
        if set(gmi_sites) == set(gmi_sites_i):
            del gmi_sites_i
        else: 
            'Input file sites do not match Strode list!'
        # load CASTNet output
        castnet = open_castnet_singyear(year, sampling_months, sampling_hours)
        # find GMI sites commensurate with CASTNet sites
        gmi_sites_fr = []
        for castnet_site in castnet_sites_fr:
            idx = np.where(np.array(castnet_sites) == castnet_site)
            gmi_sites_fr.append(np.array(gmi_sites)[idx][0]) 
        # loop through GMI and CASTNet sites and a counter in focus region 
        for gmi_site, castnet_site, counter2 in zip(gmi_sites_fr, castnet_sites_fr, 
                                                   np.arange(0, len(castnet_sites_fr), 1)):
                sitepos = np.where(np.array(gmi_sites) == gmi_site)[0]
                # find the average altitude (GMI grid box)
                gmi_sitealt = gmi_castnet_matchup[sitepos, 4]   
                gmi_sitealt = gmi_sitealt.tolist()
                gmi_sitealt = [int(x) for x in gmi_sitealt]
                gmi_sitealt = int(np.mean(gmi_sitealt))   
                # find GMI O3 
                gmi_o3 = o3[sitepos[0], :, gmi_sitealt] 
                gmi_no = no[sitepos[0], :, gmi_sitealt] 
                gmi_no2 = no2[sitepos[0], :, gmi_sitealt] 
                gmi_co = co[sitepos[0], :, gmi_sitealt] 
                # find CASTNet observations
                castnet_atsite = castnet.loc[castnet['SITE_ID'].isin([castnet_site])]
                castnet_atsite = castnet_atsite.sort_values(by = 'DATE_TIME')
                # if complete (i.e. observations for all every day contained in the months
                # in variable 'sampling_months' for all hours in 'sampling_hours'), then 
                # append
                if castnet_atsite.shape[0] > 0:
                    # fill missing CASTNet data with NaNs                    
                    if castnet_atsite.shape[0] > 2208: 
                        castnet_atsite = castnet_atsite.groupby(['DATE_TIME']).mean()
                        castnet_atsite = castnet_atsite.reindex(dates, fill_value = np.nan)                        
                    else: 
                        castnet_atsite.index = pd.DatetimeIndex(castnet_atsite['DATE_TIME'])
                        castnet_atsite = castnet_atsite.reindex(dates, 
                            fill_value = np.nan)
                    # add station's/model's diurnal cycle of trace gases to array
                    comm_castnet[counter1, counter2] = castnet_atsite['OZONE'].values.reshape(sos, 
                            len(sampling_hours)) 
                    comm_gmi_o3[counter1, counter2] = gmi_o3.reshape(sos, 
                            len(sampling_hours)) * 1e9
                    comm_gmi_no[counter1, counter2] = gmi_no.reshape(sos, 
                            len(sampling_hours)) * 1e9
                    comm_gmi_no2[counter1, counter2] = gmi_no2.reshape(sos, 
                            len(sampling_hours)) * 1e9
                    comm_gmi_co[counter1, counter2] = gmi_co.reshape(sos, 
                            len(sampling_hours)) * 1e6
    return comm_castnet, comm_gmi_o3, comm_gmi_no, comm_gmi_no2, comm_gmi_co
# # # # # # # # # # # # #
def commensurate_castnet_gmi_locations(castnet_sites_fr, sampling_months, 
                                       sampling_hours):
    """for mapping purposes, this function finds the latitudes and longitudes
    of CASTNet stations with observations for model configuration 'HindcastMR2'
    for the year 2008 (arbitrary) and their co-located (or nearly co-located)
    GMI sites.
    
    Parameters
    ----------    
    castnet_sites_fr : list
        CASTNET site names in focus region
    sampling_months : list 
        Months of interest
    sampling_hours : list 
        Hours during which trace gas concentrations are returned
        
    Returns
    ----------
    castnet_lats : list
        Latitudes of CASTNet stations with observations during 2008 in the 
        focus region, [stations in 'castnet_sites_fr',]
    castnet_lons : list
        Longitudes of CASTNet stations with observations during 2008 in the 
        focus region, [stations in 'castnet_sites_fr',]
    gmi_lats : list 
        Latitudes of GMI sites which are co-located (or nearly co-located with 
        CASTNet observations during 2008, [stations in 'castnet_sites_fr',]
    gmi_lons : list
        Longitudes of GMI sites which are co-located (or nearly co-located with 
        CASTNet observations during 2008, [stations in 'castnet_sites_fr',]    
    """
    import numpy as np
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # GMI/CASTNET site information Table S1 in Strode et al. [2015]
    gmi_castnet_matchup = np.genfromtxt(pollutants_constants.PATH_GMISTATIONINFO + 
        'gmi_castnet_matchup.csv', dtype = 'str', delimiter = ',', 
        skip_header = 1)
    gmi_sites = list(gmi_castnet_matchup[:, 0])     
    castnet_sites = list(gmi_castnet_matchup[:, 1])
    # find lat/lon of CASTNet sites in focus region
    csites = open_castnetsiting()
    # open GMI output model configuration 'HindcastMR2' for 2008 
    o3, co, no, no2, gmi_sites_i, gmi_lat, gmi_lon, times_ty = \
    commensurability.open_gmi_singyear('HindcastMR2', 2008, sampling_months,
                                       sampling_hours)
    # ensure input data matches reference table
    if set(gmi_sites) == set(gmi_sites_i):
        del gmi_sites_i
    else: 
        'Input file sites do not match Strode list!'
    # load CASTNet output
    castnet = commensurability.open_castnet_singyear(2008, sampling_months, 
                                                     sampling_hours)
    # find GMI sites commensurate with CASTNet sites
    gmi_sites_fr = []
    for castnet_site in castnet_sites_fr:
        idx = np.where(np.array(castnet_sites) == castnet_site)
        gmi_sites_fr.append(np.array(gmi_sites)[idx][0]) 
    # list creation; lists will be filled with latitude and longitudes of
    # CASTnet stations with observations in year of interest and their 
    # corresponding, co-located GMI sites
    gmi_lats, gmi_lons = [], []
    castnet_lats, castnet_lons = [], []    
    # loop through GMI and CASTNet sites and a counter in focus region 
    for gmi_site, castnet_site, counter2 in zip(gmi_sites_fr, castnet_sites_fr, 
                                               np.arange(0, len(castnet_sites_fr), 1)):
        # CASTNet observations at site of interest
        castnet_atsite = castnet.loc[castnet['SITE_ID'].isin([castnet_site])]
        castnet_atsite = castnet_atsite.sort_values(by = 'DATE_TIME')
        # find lat/lon of CASTNet station if CASTNet record has observations for 
        # the month(s) and year of interest
        if castnet_atsite.shape[0] > 0: 
            # find longitude/latitude of CASTNet station 
            csites_atsite = csites.loc[csites['SITE_ID'].isin([castnet_site])]
            lon_atsite = csites_atsite['LONGITUDE'].values[0]
            lat_atsite = csites_atsite['LATITUDE'].values[0]
            castnet_lats.append(lat_atsite) 
            castnet_lons.append(lon_atsite)
            # find longtiude/latitude of corresponding GMI site
            sitepos = np.where(np.array(gmi_sites) == gmi_site)[0]
            gmi_lats.append(gmi_lat[sitepos][0])
            gmi_lons.append(gmi_lon[sitepos][0])
    return castnet_lats, castnet_lons, gmi_lats, gmi_lons
# # # # # # # # # # # # #
def commensurate_castnet_gmi_ra(castnet_sites_fr, years, sampling_months,
                                sampling_hours):
    """function opens CASTNet and GMI trace gas concentrations for all four 
    GMI CTM model configurations (HindcastMR2, HindcastMR2-CCMI, 
    HindcastFFIgac2, HindcastFFIgac2-HighRes) using function 
    'commensurate_castnet_gmi' which produces a daily average over the hours 
    contained in variable 'sampling_hours.' Function then averages over the
    CASTNet sites to produce regionally-averaged CASTNet O3 and GMI O3, NO, 
    NO2, and CO. 
    
    Parameters
    ----------    
    castnet_sites_fr : list
        CASTNET site names in focus region
    years : list
        Years of interest
    sampling_months : list 
        Months of interest
    sampling_hours : list
        Hours of interest    
        
    Returns
    ----------
    result : dict
        Ordered dictionary containing regional mean (over CASTNet sites in 
        variable 'castnet_sites_fr') and time-averaged (over hours in variable 
        'sampling_hours') CASTNet O3 and GMI O3, NO, NO2, and CO 
    """
    import numpy as np
    # open different model configurations, averaged over 11-16 LST
    castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
    commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2', years, 
                             sampling_months, sampling_hours)
    temp, ccmi_o3, ccmi_no, ccmi_no2, ccmi_co, ccmi_gmi_sites_fr = \
    commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2-CCMI', years, 
                             sampling_months, sampling_hours)
    del temp
    temp, ffigac2_o3, ffigac2_no, ffigac2_no2, ffigac2_co, ffigac2_gmi_sites_fr = \
    commensurate_castnet_gmi(castnet_sites_fr, 'HindcastFFIgac2', years, 
                             sampling_months, sampling_hours)
    del temp
    temp, ffigac2hr_o3, ffigac2hr_no, ffigac2hr_no2, ffigac2hr_co, ffigac2hr_gmi_sites_fr = \
    commensurate_castnet_gmi(castnet_sites_fr, 'HindcastFFIgac2-HighRes', years, 
                             sampling_months, sampling_hours)
    del temp
    # average over all sites in region 
    # for O3
    castnet = np.nanmean(castnet, axis = 1)
    mr2_o3 = np.nanmean(mr2_o3, axis = 1)
    ccmi_o3 = np.nanmean(ccmi_o3, axis = 1)
    ffigac2_o3 = np.nanmean(ffigac2_o3, axis = 1)                      
    ffigac2hr_o3 = np.nanmean(ffigac2hr_o3, axis = 1)
    # for NO
    mr2_no = np.nanmean(mr2_no, axis = 1)
    ccmi_no = np.nanmean(ccmi_no, axis = 1)
    ffigac2_no = np.nanmean(ffigac2_no, axis = 1)                      
    ffigac2hr_no = np.nanmean(ffigac2hr_no, axis = 1)    
    # for NO2
    mr2_no2 = np.nanmean(mr2_no2, axis = 1)
    ccmi_no2 = np.nanmean(ccmi_no2, axis = 1)
    ffigac2_no2 = np.nanmean(ffigac2_no2, axis = 1)                      
    ffigac2hr_no2 = np.nanmean(ffigac2hr_no2, axis = 1)        
    # for CO 
    mr2_co = np.nanmean(mr2_co, axis = 1)
    ccmi_co = np.nanmean(ccmi_co, axis = 1)
    ffigac2_co = np.nanmean(ffigac2_co, axis = 1)                      
    ffigac2hr_co = np.nanmean(ffigac2hr_co, axis = 1)      
    return {'CASTNet':castnet, 'MR2 O3':mr2_o3,
            'MR2-CCMI O3':ccmi_o3, 'FFIgac2 O3':ffigac2_o3, 
            'FFIgac2-HighRes O3':ffigac2hr_o3, 'MR2 NO':mr2_no,
            'MR2-CCMI NO':ccmi_no, 'FFIgac2 NO':ffigac2_no, 
            'FFIgac2-HighRes NO':ffigac2hr_no, 'MR2 NO2':mr2_no2,
            'MR2-CCMI NO2':ccmi_no2, 'FFIgac2 NO2':ffigac2_no2,
            'FFIgac2-HighRes NO2':ffigac2hr_no2, 'MR2 CO':mr2_co,
            'MR2-CCMI CO':ccmi_co, 'FFIgac2 CO':ffigac2_co,
            'FFIgac2-HighRes CO':ffigac2hr_co}
# # # # # # # # # # # # #
def commensurate_castnet_gmi_diurnal_ra(castnet_sites_fr, years, 
                                        sampling_months):
    """function opens CASTNet and GMI trace gas concentrations for all four 
    GMI CTM model configurations (HindcastMR2, HindcastMR2-CCMI, 
    HindcastFFIgac2, HindcastFFIgac2-HighRes) and averages over the years, 
    days, and sites to produce average diurnal cycles (i.e. 0 - 23 UTC) for 
    CASTNet O3 and GMI O3, NO, NO2, and CO. 
    
    Parameters
    ----------    
    castnet_sites_fr : list
        CASTNET site names in focus region
    years : list
        Years of interest
    sampling_months : list 
        Months of interest
    
    Returns
    ----------
    result : dict
        Ordered dictionary containing regionally- and time-averaged diurnal 
        cycles of CASTNet O3 and GMI O3, NO, NO2, and CO over the CASTNet sites 
        and months contained in variables 'castnet_sites_fr' and 
        'sampling_months.'
    """
    import numpy as np
    # open different model configurations, diurnal cycle
    castnet_d, mr2_o3_d, mr2_no_d, mr2_no2_d, mr2_co_d = \
    commensurate_castnet_gmi_diurnal(castnet_sites_fr, 'HindcastMR2', 
                                                      years, sampling_months)
    temp, ccmi_o3_d, ccmi_no_d, ccmi_no2_d, ccmi_co_d = \
    commensurate_castnet_gmi_diurnal(castnet_sites_fr, 'HindcastMR2-CCMI', 
                                                      years, sampling_months)
    del temp
    temp, ffigac2_o3_d, ffigac2_no_d, ffigac2_no2_d, ffigac2_co_d = \
    commensurate_castnet_gmi_diurnal(castnet_sites_fr, 'HindcastFFIgac2', 
                                                      years, sampling_months)
    del temp
    temp, ffigac2hr_o3_d, ffigac2hr_no_d, ffigac2hr_no2_d, ffigac2hr_co_d = \
    commensurate_castnet_gmi_diurnal(castnet_sites_fr, 'HindcastFFIgac2-HighRes', 
                                                      years, sampling_months)
    del temp
    # average over years, CASTNet sites, and days in variable 'sampling_months'
    # in region to yield a single diurnal curve
    # for O3
    castnet_d = np.nanmean(castnet_d, axis = tuple(range(0, 3)))
    mr2_o3_d = np.nanmean(mr2_o3_d, axis = tuple(range(0, 3)))
    ccmi_o3_d = np.nanmean(ccmi_o3_d, axis = tuple(range(0, 3)))
    ffigac2_o3_d = np.nanmean(ffigac2_o3_d, axis = tuple(range(0, 3)))
    ffigac2hr_o3_d = np.nanmean(ffigac2hr_o3_d, axis = tuple(range(0, 3)))
    # for NO
    mr2_no_d = np.nanmean(mr2_no_d, axis = tuple(range(0, 3)))
    ccmi_no_d = np.nanmean(ccmi_no_d, axis = tuple(range(0, 3)))
    ffigac2_no_d = np.nanmean(ffigac2_no_d, axis = tuple(range(0, 3)))
    ffigac2hr_no_d = np.nanmean(ffigac2hr_no_d, axis = tuple(range(0, 3)))    
    # for NO2
    mr2_no2_d = np.nanmean(mr2_no2_d, axis = tuple(range(0, 3)))
    ccmi_no2_d = np.nanmean(ccmi_no2_d, axis = tuple(range(0, 3)))
    ffigac2_no2_d = np.nanmean(ffigac2_no2_d, axis = tuple(range(0, 3)))
    ffigac2hr_no2_d = np.nanmean(ffigac2hr_no2_d, axis = tuple(range(0, 3)))        
    # for CO
    mr2_co_d = np.nanmean(mr2_co_d, axis = tuple(range(0, 3)))
    ccmi_co_d = np.nanmean(ccmi_co_d, axis = tuple(range(0, 3)))
    ffigac2_co_d = np.nanmean(ffigac2_co_d, axis = tuple(range(0, 3)))
    ffigac2hr_co_d = np.nanmean(ffigac2hr_co_d, axis = tuple(range(0, 3)))   
    # return values using a dictionary
    return {'CASTNet':castnet_d, 'MR2 O3':mr2_o3_d,
            'MR2-CCMI O3':ccmi_o3_d, 'FFIgac2 O3':ffigac2_o3_d, 
            'FFIgac2-HighRes O3':ffigac2hr_o3_d, 'MR2 NO':mr2_no_d,
            'MR2-CCMI NO':ccmi_no_d, 'FFIgac2 NO':ffigac2_no_d, 
            'FFIgac2-HighRes NO':ffigac2hr_no_d, 'MR2 NO2':mr2_no2_d,
            'MR2-CCMI NO2':ccmi_no2_d, 'FFIgac2 NO2':ffigac2_no2_d,
            'FFIgac2-HighRes NO2':ffigac2hr_no2_d, 'MR2 CO':mr2_co_d,
            'MR2-CCMI CO':ccmi_co_d, 'FFIgac2 CO':ffigac2_co_d,
            'FFIgac2-HighRes CO':ffigac2hr_co_d}
# # # # # # # # # # # # #
def open_perturbed_emissions():
    """function opens perturbed emissions inventories used in GMI CTM 
    sensitivity simulations for the measuring period JJA 2008-2010 along with 
    latitude and longitude coordinates for the emissions inventory. The 
    emissions inventories retrieved are only for NO from the fossil fuel sector
    and vary on a daily basis, scaling with the observed temperature-industrial
    NOx relationship derived from MERRA-1/CEMS. These inventories were used as 
    inputs to the CTM hindcast simulation 'EGU_T.'

    Parameters
    ----------    
    None 

    Returns
    ----------      
    NO_ff_all : list
        JJA 2008-2010 NO from fossil fuels emissions inventory; emissions
        scale with temperature. Every element in list (numpy.ndarray) is 
        a month's worth of daily-varying NO emissions.
    lon : list
        GMI longitude grid, 1.25˚ resolution, units of degrees east, [lon,]
    lat : list
        GMI latitude grid, 1˚ resolution, units of degrees north, [lat,]  
    """
    import numpy as np
    from netCDF4 import Dataset    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    import glob
    # month abbreviations
    months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 
              'SEP', 'OCT', 'NOV', 'DEC']
    years = ['2008', '2009', '2010']
    # path to data holdings
    path = pollutants_constants.PATH_GMIEMISSIONS
    mfstring = 'hindcast_IAVanthQFED.daily.1x1.25_*_*_NOscaling.nc'
    infiles = glob.glob(path + mfstring)
    # sort file strings by both month and year; from 
    # stackoverflow.com/questions/33767899/sort-list-according-given-text-of-months
    infiles = sorted(infiles, key=lambda x: (years.index(x.split('_')[2]), 
                          months.index(x.split('_')[3])))
    # list creation; lists will be filled on subsequent iterations in loop
    lat, lon = [], []
    NO_ff_all = []
    # open perturbed NOx emissions inventories for each summer month in 
    # measuring period 2008-2010
    i = 0 
    for infile in infiles: 
        infile = Dataset(infile, 'r')
        # only fetch lat/lon on first iteration 
        if i == 0:
            lat.append(infile.variables['latitude_dim'][:])
            lon.append(infile.variables['longitude_dim'][:])
        # load chemical species names    
        species_raw = infile.variables['species'][:]
        species = []
        for row in np.arange(0, np.shape(species_raw)[0]):
            temp = []
            for element in species_raw[row]:
                temp.append(element.decode('UTF-8'))
            si = ''.join(temp)
            si = si.rstrip()
            species.append(si[:])      
        del species_raw
        # NO from fossil fuels field
        NO_ff = infile['emiss_2d'][:, np.where(np.array(species) == 'NO_ff')[0][0], :, :]
        NO_ff_all.append(NO_ff)
        i = i + 1
    return NO_ff_all, lon, lat
# # # # # # # # # # # # #    
def open_unperturbed_emissions(prefix):
    """same as function 'open_perturbed_emissions' but for unperturbed 
    emissions inventory used in HindcastMR2 (control). 
    
    Parameters
    ----------    
    prefix : str
        File prefix of emissions inventory (i.e. for MR2-EGU run prefix is
        IAVanthGFED4; for FFIgac2 run prefix is 'IAVanthGFED3gcEF')

    Returns
    ----------      
    NO_ff_all : list
        JJA 2008-2010 NO from fossil fuels emissions inventory; emissions
        scale with temperature. Every element in list (numpy.ndarray) 
        represents monthly mean emissions
    lon : list
        GMI longitude grid, 1.25˚ resolution, units of degrees east, [lon,]
    lat : list
        GMI latitude grid, 1˚ resolution, units of degrees north, [lat,]  
    """
    import numpy as np
    from netCDF4 import Dataset    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    import glob
    # path to data holdings
    path = pollutants_constants.PATH_GMIEMISSIONS
    mfstring = 'emisc_*_m_%s.nc' %(prefix)
    infiles = glob.glob(path + mfstring)
    infiles.sort()
    # list creation; lists will be filled on subsequent iterations in loop
    lat, lon = [], []
    NO_ff_all = []
    # open emissions inventories for measuring period 2008-2010
    i = 0 
    for infile in infiles: 
        infile = Dataset(infile, 'r')
        # only fetch lat/lon on first iteration 
        if i == 0:
            lat.append(infile.variables['latitude_dim'][:])
            lon.append(infile.variables['longitude_dim'][:])
        # load chemical species names    
        species_raw = infile.variables['species'][:]
        species = []
        for row in np.arange(0, np.shape(species_raw)[0]):
            temp = []
            for element in species_raw[row]:
                temp.append(element.decode('UTF-8'))
            si = ''.join(temp)
            si = si.rstrip()
            species.append(si[:])      
        del species_raw
        # NO from fossil fuels field
        if prefix == '2x2.5_IAVanthGFED3gcEF':
            NO_ff = infile['emiss'][5:8, 
                          np.where(np.array(species) == 'NO_ff')[0][0], 0, :, :]
        if prefix == '1x1.25_IAVanthGFED4':
            NO_ff = infile['emiss_2d'][5:8, 
                          np.where(np.array(species) == 'NO_ff')[0][0], :, :]
        NO_ff_all.append(NO_ff)    
        i = i + 1
    return NO_ff_all, lon, lat  
# # # # # # # # # # # # #
def commensurate_emiss_perturbed(comm_castnet, castnet_sites_fr, years, 
                                 sampling_months):
    """using function 'open_perturbed_emissions' this function loads GMI CTM
    perturbed emissions inventories and finds grid cells commensurate with 
    CASTNet sites with data and retrieves them. 

    Parameters
    ----------    
    comm_castnet : numpy.ndarray
        Existing CASTNet O3 observations at each station contained in variable 
        'castnet_sites_fr' for years in measuring period, units of ppbv, 
        [years in measuring period, stations in 'castnet_sites_fr', 
        days in months in 'sampling_months']
    castnet_sites_fr : list
        CASTNET site names in focus region
    years : list
        Years of interest
    sampling_months : list 
        Months of interest

    Returns
    ----------
    comm_emiss : numpy.ndarray
        GMI CTM NO from fossil fuel sector emissions inventory at grid cells
        co-located (or nearly colocated) with CASTNet stations, [years in 
        measuring period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']
    emiss_lats_fr : list  
        Longitudes of CTM GMI emissions inventory grid cells which are 
        co-located (or nearly co-located with CASTNet observations during years 
        of interest
    emiss_lons_fr : list
        Latitudes of CTM GMI emissions inventory grid cells which are 
        co-located (or nearly co-located with CASTNet observations during years 
        of interest     
    """
    import numpy as np
    from calendar import monthrange
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    # load perturbed emissions
    NO_ff, lon, lat = open_perturbed_emissions()
    # convert longitude from (0 to 360) to (-180 - 180)    
    lon = np.mod(lon[0] - 180.0, 360.0) - 180.0
    lat = lat[0]
    # load CASTNet stations' locations
    csites = open_castnetsiting()    
    # find days in summer
    sos = []
    sos += [(monthrange(years[0], x)[1]) for x in sampling_months]    
    sos = np.sum(sos)
    # stack emissions array and change dimensions such that first dimension 
    # is the year, second is day of year; n.b. this reshaping is based on a 92 
    # day summer
    NO_ff = np.vstack(NO_ff)
    NO_ff = np.reshape(NO_ff, [len(years), sos, len(lat), len(lon)])
    # create empty array to be filled with NO emissions commensurate with 
    # CASTNet sites; n.b. maybe explore aggregating emissions from upwind sites 
    # into a single value? 
    comm_emiss = np.empty([len(years), len(castnet_sites_fr), sos], 
                           dtype = float)
    comm_emiss[:] = np.nan
    # create lists to be filled with the locations of emissions' grid cells
    emiss_lats_fr, emiss_lons_fr = [], []        
    # loop through years in measuring period 
    for year, counter1 in zip(years, np.arange(0, len(years), 1)): 
        # in a given year, loop through CASTNet sites in focus region and 
        # determine if data exists
        for site, castnet_data, counter2 in zip(castnet_sites_fr, 
                                                comm_castnet[counter1], 
                                                np.arange(0, len(castnet_sites_fr), 1)):
            # if data exists, then find MERRA T2m commensurate with CASTNet site
            if np.where(np.isnan(castnet_data) == True)[0].shape[0] != sos:
                # find longitude/latitude of CASTNET station 
                csites_atsite = csites.loc[csites['SITE_ID'].isin([site])]
                lon_atsite = csites_atsite['LONGITUDE'].values[0]
                lat_atsite = csites_atsite['LATITUDE'].values[0]
                # find closest grid cell to CASTNet site in emissions inventory
                lat_idx = geo_idx(lat_atsite, lat)
                lon_idx = geo_idx(lon_atsite, lon)
                # NO emissions co-located with CASTNet site
                emiss_atsite = NO_ff[counter1, :, lat_idx, lon_idx]
                comm_emiss[counter1, counter2] = emiss_atsite
                emiss_lats_fr.append(lat[lat_idx]) 
                emiss_lons_fr.append(lon[lon_idx])
    print('CTM emission data loaded!')        
    return comm_emiss, emiss_lats_fr, emiss_lons_fr
# # # # # # # # # # # # #
def commensurate_emiss_unperturbed(comm_castnet, castnet_sites_fr, years, 
                                   sampling_months, prefix):
    """using function 'open_unperturbed_emissions' this function loads GMI CTM
    unperturbed emissions inventories and finds grid cells commensurate with 
    CASTNet sites with data and retrieves them. Since emissions in unperturbed 
    emissions inventory only vary on a monthly basis (monthly mean emissions), 
    function expands these mean emissions by the number of days in each 
    summer month.

    Parameters
    ----------    
    comm_castnet : numpy.ndarray
        Existing CASTNet O3 observations at each station contained in variable 
        'castnet_sites_fr' for years in measuring period, units of ppbv, 
        [years in measuring period, stations in 'castnet_sites_fr', 
        days in months in 'sampling_months']
    castnet_sites_fr : list
        CASTNET site names in focus region
    years : list
        Years of interest
    sampling_months : list 
        Months of interest
    prefix : str
        File prefix of emissions inventory (i.e. for MR2-EGU run prefix is
        IAVanthGFED4; for FFIgac2 run prefix is 'IAVanthGFED3gcEF')

    Returns
    ----------
    comm_emiss : numpy.ndarray
        GMI CTM NO from fossil fuel sector emissions inventory at grid cells
        co-located (or nearly colocated) with CASTNet stations, [years in 
        measuring period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']    
    """
    import numpy as np
    from calendar import monthrange
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    # open unperturbed emissions
    NO_ff, lon, lat = open_unperturbed_emissions(prefix)
    # convert longitude from (0 to 360) to (-180 - 180)    
    lon = np.mod(lon[0] - 180.0, 360.0) - 180.0
    lat = lat[0]
    # load CASTNet stations' locations
    csites = open_castnetsiting()    
    # stack emissions array and reshape (i.e. [number of days, number of months, 
    # lat, lon])
    NO_ff = np.vstack(NO_ff)
    NO_ff = np.reshape(NO_ff, [len(years), len(sampling_months), len(lat),
                               len(lon)])
    # find days in summer
    sos = []
    sos += [(monthrange(years[0], x)[1]) for x in sampling_months]    
    sos = np.sum(sos)
    # create empty array to be filled with NO emissions commensurate with 
    # CASTNet sites
    comm_emiss = np.empty([len(years), len(castnet_sites_fr), sos], 
                           dtype = float)
    comm_emiss[:] = np.nan
    # loop through years in measuring period 
    for year, counter1 in zip(years, np.arange(0, len(years), 1)): 
        # in a given year, loop through CASTNet sites in focus region and 
        # determine if data exists
        for site, castnet_data, counter2 in zip(castnet_sites_fr, 
                                                comm_castnet[counter1], 
                                                np.arange(0, len(castnet_sites_fr), 1)):
            # if data exists, then find MERRA T2m commensurate with CASTNet site
            if np.where(np.isnan(castnet_data) == True)[0].shape[0] != sos:
                # find longitude/latitude of CASTNET station 
                csites_atsite = csites.loc[csites['SITE_ID'].isin([site])]
                lon_atsite = csites_atsite['LONGITUDE'].values[0]
                lat_atsite = csites_atsite['LATITUDE'].values[0]
                # find closest grid cell to CASTNet site in emissions inventory
                lat_idx = geo_idx(lat_atsite, lat)
                lon_idx = geo_idx(lon_atsite, lon)
                # NO emissions co-located with CASTNet site
                emiss_atsite = NO_ff[counter1, :, lat_idx, lon_idx]
                # expand monthly mean NO emissions at site by the number of 
                # days in each summer month; n.b. this assumes JJA
                emiss_all = []
                for emiss, md in zip(emiss_atsite, [30, 31, 31]):
                    emiss_all.append(np.repeat(emiss, md))
                comm_emiss[counter1, counter2] = np.hstack(emiss_all)            
    print('CTM emission data loaded!')        
    return comm_emiss
# # # # # # # # # # # # #
def commensurate_geos_gmi(castnet_sites_fr, sampling_hours):
    """function finds trace gas concentration from 2 hourly GEOS-Chem output 
    from Katie Travis' simulation at CASTNet sites. Concentrations are 
    averaged over hours in variable 'sampling_hours' to produce a daily 
    average. n.b. function produces a daily average of hourly 
    values assuming that len(sampling_hours) = 6. Change code if this isn't the
    the case. 
    
    Parameters
    ----------    
    castnet_sites_fr : list
        CASTNET site names in focus region
    sampling_hours : list 
        Hours during which trace gas concentrations are returned

    Returns
    ----------      
    comm_castnet : numpy.ndarray
        Existing CASTNet O3 observations at each station contained in variable 
        'castnet_sites_fr' for years in measuring period, units of ppbv, 
        [stations in 'castnet_sites_fr', days in months in 'sampling_months']
    comm_geos_o3 : numpy.ndarray
        GEOS-Chem O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per billion volume, 
        [stations in 'castnet_sites_fr', days in months in 'sampling_months']
    comm_geos_no : numpy.ndarray
        GEOS-Chem NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per billion volume, 
        [stations in 'castnet_sites_fr', days in months in 'sampling_months']
    comm_geos_no2 : numpy.ndarray
        GEOS-Chem NO2 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per billion volume, 
        [stations in 'castnet_sites_fr', days in months in 'sampling_months']
    comm_geos_co : numpy.ndarray     
        GEOS-Chem CO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of parts per billion volume, 
        [stations in 'castnet_sites_fr', days in months in 'sampling_months']
    geos_lats_atsite : list
        Latitude of GEOS-Chem grid cells co-located with CASTNet sites
    geos_lons_atsite : list
        Longitude of GEOS-Chem grid cells co-located with CASTNet sites    
    """
    import numpy as np
    import pandas as pd
    from calendar import monthrange
    from netCDF4 import Dataset, num2date
    import warnings    
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    year = 2013
    sampling_months = [8, 9]
    # list creation
    geos_lats_atsite = []
    geos_lons_atsite = []
    # find number of hours in variable 'sampling_months', n.b. this assumes that
    # the number of hours in 'sampling_months' is constant from year to year,
    # so it wouldn't work for a leap year 
    sos = []
    sos += [(monthrange(year, x)[1]) for x in sampling_months]    
    sos = np.sum(sos)
    # create empty arrays to be filled with commensurate CASTNet and GMI O3 
    # concentrations for multi-year period
    comm_castnet = np.empty([len(castnet_sites_fr), sos], dtype = float)
    comm_castnet[:] = np.nan
    comm_geos_o3 = np.empty([len(castnet_sites_fr), sos], dtype = float)
    comm_geos_o3[:] = np.nan
    comm_geos_no = np.empty([len(castnet_sites_fr), sos], dtype = float)
    comm_geos_no[:] = np.nan
    comm_geos_no2 = np.empty([len(castnet_sites_fr), sos], dtype = float)
    comm_geos_no2[:] = np.nan
    comm_geos_co = np.empty([len(castnet_sites_fr), sos], dtype = float)
    comm_geos_co[:] = np.nan    
    # load GEOS-Chem output for July - August 2013
    path_geos = '/Users/ghkerr/phd/GMI/data/'
    geos = Dataset(path_geos + 'ts2013_concat.bpch.nc', 'r')
    geos_o3 = geos.variables['o3'][:]
    geos_no = geos.variables['no'][:]
    geos_no2 = geos.variables['no2'][:]
    geos_co = geos.variables['co'][:]
    geos_lat = geos.variables['latitude_dim'][:]
    geos_lon = geos.variables['longitude_dim'][:]
    # convert GEOS-Chem time to datetime
    geos_times = geos.variables['time_dim'][:]
    time_unit = 'hours since 1985-01-01 00:00:00 UTC'
    time_cal = u"gregorian" # or standard
    timevar = []
    timevar.append(num2date(geos_times, units = time_unit, 
                            calendar = time_cal))
    # find hours of GEOS-Chem timestamps 
    geos_hours = []
    for it in timevar[0]: 
        geos_hours.append(it.hour)
    # find indices of GEOS-Chem 2-hourly fields in sampling_hours
    geos_idx = []    
    for hour in sampling_hours:
        geos_idx.append(np.where(np.array(geos_hours) == hour))
    geos_idx = np.sort(np.hstack(geos_idx))
    # load CASTNet output
    castnet = open_castnet_singyear(year, sampling_months, sampling_hours)
    # load CASTNet siting information 
    csites = open_castnetsiting()
    # times in year of interest with hourly timestep, retrieve only 
    # hours in variable 'sampling_hours' in ET
    dates = pd.date_range('%s-01-%s' %(sampling_months[0], year), 
        '%s-01-%s' %(sampling_months[-1] + 1, year), freq = '1H')[:-1]                
    dates = dates[np.where(dates.hour.isin([x-4 for x in sampling_hours]))[0]]
    # loop through CASTNet sites in region 
    for castnet_site, counter in zip(castnet_sites_fr, 
                                     np.arange(0, len(castnet_sites_fr), 1)):
        castnet_atsite = castnet.loc[castnet['SITE_ID'].isin([castnet_site])]
        castnet_atsite = castnet_atsite.sort_values(by = 'DATE_TIME')
        # only consider CASTNet sites with observations 
        if castnet_atsite.shape[0] > 0: 
            # fill in potentially missing data with NaN
            castnet_atsite.index = pd.DatetimeIndex(castnet_atsite['DATE_TIME'])
            castnet_atsite = castnet_atsite.reindex(dates, fill_value = np.nan)
            # find lat/lon at CASTNet site 
            castnet_siteinfo = csites.loc[csites['SITE_ID'].isin([castnet_site])]
            castnet_lat = castnet_siteinfo['LATITUDE'].values[0]
            castnet_lon = castnet_siteinfo['LONGITUDE'].values[0]
            # find lat/lon index of closest GEOS-Chem grid cell
            geos_lat_idx = geo_idx(castnet_lat, geos_lat)
            geos_lon_idx = geo_idx(castnet_lon, geos_lon)
            # save off lat/lon of GEOS-Chem grid cells
            geos_lats_atsite.append(geos_lat[geos_lat_idx])
            geos_lons_atsite.append(geos_lon[geos_lon_idx])
            # GEOS-Chem output at grid cell closest to CASTNet site at time of 
            # sampling hours
            geos_o3_atsite = geos_o3[geos_idx, geos_lat_idx, geos_lon_idx][0]
            geos_no_atsite = geos_no[geos_idx, geos_lat_idx, geos_lon_idx][0]
            geos_no2_atsite = geos_no2[geos_idx, geos_lat_idx, geos_lon_idx][0]
            geos_co_atsite = geos_co[geos_idx, geos_lat_idx, geos_lon_idx][0]
            # to produce daily average with hourly or 2 hourly values, average 
            # every 6 values for CASTNet (n.b. this only works when 
            # len(sampling_hours) = 6); average every 3 values for GEOS-Chem 
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Mean of empty slice")
                comm_castnet[counter] = np.nanmean(np.array(
                    castnet_atsite['OZONE'].values).reshape(-1, 6), axis = 1)
                comm_geos_o3[counter] = np.nanmean(np.array(
                    geos_o3_atsite).reshape(-1, 3), axis = 1)      
                comm_geos_no[counter] = np.nanmean(np.array(
                    geos_no_atsite).reshape(-1, 3), axis = 1)   
                comm_geos_no2[counter] = np.nanmean(np.array(
                    geos_no2_atsite).reshape(-1, 3), axis = 1) 
                comm_geos_co[counter] = np.nanmean(np.array(
                    geos_co_atsite).reshape(-1, 3), axis = 1)   
    print('GEOS-Chem data loaded!')
    return (comm_castnet, comm_geos_o3, comm_geos_no, comm_geos_no2, 
            comm_geos_co, geos_lats_atsite, geos_lons_atsite)
# # # # # # # # # # # # #        
def commensurate_geos_gmi_diurnal(castnet_sites_fr):
    """function opens diurnal curves O3 from CASTNet sites contained in 
    'castnet_sites_fr' and finds co-located (or nearly co-located) trace 
    gas concentrations from GEOS-Chem and outputs their diurnal curves; n.b. 
    GEOS-Chem output is 2-hourly, so for a particular site, only 12 
    concentrations comprise the daily diurnal curve. 

    Parameters
    ----------    
    castnet_sites_fr : list
        CASTNET site names in focus region

    Returns
    ----------      
    comm_castnet : numpy.ndarray
        Diurnal measurements of existing CASTNet O3 observations at each 
        station contained in variable 'castnet_sites_fr' for Aug - Sep 2013, 
        units of ppbv, [stations in 'castnet_sites_fr', 61, 24]    
    comm_geos_o3 : numpy.ndarray
        Diurnal GEOS-Chem O3 concentrations co-located (or nearly colocated) 
        with corresponding CASTNet stations, units of ppbv, [stations in 
        'castnet_sites_fr', 61, 12]    
    comm_geos_no : numpy.ndarray
        Diurnal GEOS-Chem NO concentrations co-located (or nearly colocated) 
        with corresponding CASTNet stations, units of ppbv, [stations in 
        'castnet_sites_fr', 61, 12]   
    comm_geos_no2 : numpy.ndarray
        Diurnal GEOS-Chem NO2 concentrations co-located (or nearly colocated) 
        with corresponding CASTNet stations, units of ppbv, [stations in 
        'castnet_sites_fr', 61, 12]  
    comm_geos_co : numpy.ndarray
        Diurnal GEOS-Chem CO concentrations co-located (or nearly colocated) 
        with corresponding CASTNet stations, units of ppbv, [stations in 
        'castnet_sites_fr', 61, 12]  
    """
    import numpy as np
    import pandas as pd
    from calendar import monthrange
    from netCDF4 import Dataset
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    year = 2013
    sampling_months = [8, 9]
    sampling_hours = list(np.arange(0, 24, 1))
    # find number of hours in variable 'sampling_months', n.b. this assumes that
    # the number of hours in 'sampling_months' is constant from year to year,
    # so it wouldn't work for a leap year 
    sos = []
    sos += [(monthrange(year, x)[1]) for x in sampling_months]    
    sos = np.sum(sos)
    # create empty arrays to be filled with commensurate CASTNet and GEOS-Chem
    # trace gas concentrations (n.b. GEOS-Chem output is 2-hourly)
    comm_castnet = np.empty([len(castnet_sites_fr), sos, 24], dtype = float)
    comm_castnet[:] = np.nan
    comm_geos_o3 = np.empty([len(castnet_sites_fr), sos, 12], dtype = float)
    comm_geos_o3[:] = np.nan
    comm_geos_no = np.empty([len(castnet_sites_fr), sos, 12], dtype = float)
    comm_geos_no[:] = np.nan
    comm_geos_no2 = np.empty([len(castnet_sites_fr), sos, 12], dtype = float)
    comm_geos_no2[:] = np.nan
    comm_geos_co = np.empty([len(castnet_sites_fr), sos, 12], dtype = float)
    comm_geos_co[:] = np.nan    
    # load GEOS-Chem output for July - August 2013
    path_geos = '/Users/ghkerr/phd/GMI/data/'
    geos = Dataset(path_geos + 'ts2013_concat.bpch.nc', 'r')
    geos_o3 = geos.variables['o3'][:]
    geos_no = geos.variables['no'][:]
    geos_no2 = geos.variables['no2'][:]
    geos_co = geos.variables['co'][:]
    geos_lat = geos.variables['latitude_dim'][:]
    geos_lon = geos.variables['longitude_dim'][:]
    # load CASTNet output
    castnet = open_castnet_singyear(year, sampling_months, sampling_hours)
    # load CASTNet siting information 
    csites = open_castnetsiting()
    # times in year of interest with hourly timestep, retrieve only 
    # hours in variable 'sampling_hours' in ET
    dates = pd.date_range('%s-01-%s' %(sampling_months[0], year), 
        '%s-01-%s' %(sampling_months[-1] + 1, year), freq = '1H')[:-1]                
    # loop through CASTNet sites in region 
    for castnet_site, counter in zip(castnet_sites_fr, 
                                     np.arange(0, len(castnet_sites_fr), 1)):
        castnet_atsite = castnet.loc[castnet['SITE_ID'].isin([castnet_site])]
        castnet_atsite = castnet_atsite.sort_values(by = 'DATE_TIME')
        # only consider CASTNet sites with observations 
        if castnet_atsite.shape[0] > 0: 
            # fill in potentially missing data with NaN
            castnet_atsite.index = pd.DatetimeIndex(castnet_atsite['DATE_TIME'])
            castnet_atsite = castnet_atsite.reindex(dates, fill_value = np.nan)
            # find lat/lon at CASTNet site 
            castnet_siteinfo = csites.loc[csites['SITE_ID'].isin([castnet_site])]
            castnet_lat = castnet_siteinfo['LATITUDE'].values[0]
            castnet_lon = castnet_siteinfo['LONGITUDE'].values[0]
            # find lat/lon index of closest GEOS-Chem grid cell
            geos_lat_idx = geo_idx(castnet_lat, geos_lat)
            geos_lon_idx = geo_idx(castnet_lon, geos_lon)
            # GEOS-Chem output at grid cell closest to CASTNet site at time of 
            # sampling hours
            geos_o3_atsite = geos_o3[:, geos_lat_idx, geos_lon_idx]
            geos_no_atsite = geos_no[:, geos_lat_idx, geos_lon_idx]
            geos_no2_atsite = geos_no2[:, geos_lat_idx, geos_lon_idx]
            geos_co_atsite = geos_co[:, geos_lat_idx, geos_lon_idx]
            # add station's/model's diurnal cycle of trace gases to array
            comm_castnet[counter] = castnet_atsite['OZONE'].values.reshape(sos, 
                    len(sampling_hours)) 
            comm_geos_o3[counter] = geos_o3_atsite.reshape(sos, 12)
            comm_geos_no[counter] = geos_no_atsite.reshape(sos, 12)
            comm_geos_no2[counter] = geos_no2_atsite.reshape(sos, 12)
            comm_geos_co[counter] = geos_co_atsite.reshape(sos, 12)
    return comm_castnet, comm_geos_o3, comm_geos_no, comm_geos_no2, comm_geos_co                   
# # # # # # # # # # # # #
def commensurate_aqstracegas_siting(castnet_sites_fr, years, sampling_months,
                                    sampling_hours):
    """same as 'commensurate_aqstracegas' but function separates AQS CO, NO2, 
    and O3 based on their siting environment (i.e. rural, suburban, or 
    rural) for all AQS stations within the radius from CASTNet stations 
    defined in 'searchrad.'
    
    Parameters
    ----------    
    castnet_sites_fr : list
        CASTNET site names in focus region
    years : list
        Years of interest
    sampling_months : list 
        Months of interest
    sampling_hours : list
        Hours of interest
        
    Returns
    ----------
    comm_co_su : numpy.ndarray
        AQS CO observations from suburban locations located within a bounding 
        box defined by variable 'searchrad' of the location of CASTNet station 
        averaged over the hours in variable 'sampling_hours,' units of parts 
        per million, [years in measuring period, stations in 
        'castnet_sites_fr', days in months in 'sampling_months']    
    comm_co_r : numpy.ndarray
        AQS CO observations from rural locations located within a bounding 
        box defined by variable 'searchrad' of the location of CASTNet station 
        averaged over the hours in variable 'sampling_hours,' units of parts 
        per million, [years in measuring period, stations in 
        'castnet_sites_fr', days in months in 'sampling_months']       
    comm_co_u : numpy.ndarray
        AQS CO observations from urban locations located within a bounding 
        box defined by variable 'searchrad' of the location of CASTNet station 
        averaged over the hours in variable 'sampling_hours,' units of parts 
        per million, [years in measuring period, stations in 
        'castnet_sites_fr', days in months in 'sampling_months']            
    coords_co_su : list
        Coordinates of AQS stations (may not be unique) which measure CO in 
        suburban locations within bounding boxes defined by CASTNet stations      
    coords_co_r : list
        Coordinates of AQS stations (may not be unique) which measure CO in 
        rural locations within bounding boxes defined by CASTNet stations   
    coords_co_u : list
        Coordinates of AQS stations (may not be unique) which measure CO in 
        urban locations within bounding boxes defined by CASTNet stations       
    comm_no2_su : numpy.ndarray
        Same as 'comm_co_su' but for suburban NO2 measurements, units of 
        parts per billion
    comm_no2_r : numpy.ndarray
        Same as 'comm_co_r' but for rural NO2 measurements, units of parts per 
        billion    
    comm_no2_u : numpy.ndarray
        Same as 'comm_co_u' but for urban NO2 measurements, units of parts per 
        billion           
    coords_no2_su : list
        Same as 'coords_co_su' but for suburban NO2 stations    
    coords_no2_r : list
        Same as 'coords_co_su' but for rural NO2 stations        
    coords_no2_u : list
        Same as 'coords_co_su' but for urban NO2 stations        
    comm_o3_su : numpy.ndarray
        Same as 'comm_co_su' but for suburban O3 measurements, units of parts 
        per million       
    comm_o3_r : numpy.ndarray
        Same as 'comm_co_r' but for rural O3 measurements, units of parts per 
        million         
    comm_o3_u : numpy.ndarray
        Same as 'comm_co_u' but for urban O3 measurements, units of parts per 
        million                
    coords_o3_su : list
        Same as 'coords_co_su' but for suburban O3 stations
    coords_o3_r : list
        Same as 'coords_co_r' but for rural O3 stations
    coords_o3_u : list
        Same as 'coords_co_u' but for urban O3 stations
    """
    import numpy as np
    from calendar import monthrange
    import pandas as pd
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI')
    import commensurability
    # search radius 
    searchrad = 1.0    
    # open AQS NO2/CO observations 
    COdf, NO2df, O3df = commensurability.load_aqshourly(years)
    COdf['Date GMT'] = pd.to_datetime(COdf['Date GMT'])
    NO2df['Date GMT'] = pd.to_datetime(NO2df['Date GMT'])
    O3df['Date GMT'] = pd.to_datetime(O3df['Date GMT'])      
    # select summer and afternoon observations
    CO = COdf.loc[COdf['Date GMT'].dt.month.isin(sampling_months)]
    CO = CO.loc[CO['Time GMT'].isin(['%d:00' %(x) for x in sampling_hours])] 
    NO2 = NO2df.loc[NO2df['Date GMT'].dt.month.isin(sampling_months)]
    NO2 = NO2.loc[NO2['Time GMT'].isin(['%d:00' %(x) for x in sampling_hours])] 
    O3 = O3df.loc[O3df['Date GMT'].dt.month.isin(sampling_months)]
    O3 = O3.loc[O3['Time GMT'].isin(['%d:00' %(x) for x in sampling_hours])] 
    # open CASTNet observations for particular model configuration (in this 
    # case 'HindcastMR2') to see which stations have observations for the 
    # measuring period
    comm_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
    commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2', 
                                              years, sampling_months, 
                                              sampling_hours)  
    # find lat/lon of CASTNet sites in focus region
    csites = commensurability.open_castnetsiting()
    # read AQS site file containing information about the siting of AQS
    # sites
    datapath = '/Users/ghkerr/phd/aqs_station_siting/data/'
    aqs_sites = pd.read_csv((datapath + 'aqs_sites.csv'), header = 0)
    site_cols = ['State Code', 'County Code',	 'Site Number', 'Latitude',
                 'Longitude', 'Datum', 'Elevation', 'Land Use',	
                 'Location Setting', 'Site Established Date', 'Site Closed Date',
                 'Met Site State Code', 'Met Site County Code',
                 'Met Site Site Number', 'Met Site Type', 'Met Site Distance',
                 'Met Site Direction', 'GMT Offset', 'Owning Agency',
                 'Local Site Name', 'Address', 'Zip Code', 'State Name',
                 'County Name', 'City Name', 'CBSA Name', 'Tribe Name',	
                 'Extraction Date']
    aqs_sites = pd.DataFrame(aqs_sites, columns = site_cols) 
    # total number of days in each year
    sos = []
    sos += [(monthrange(years[0], x)[1]) for x in sampling_months]    
    sos = np.sum(sos)
    # create empty arrays to be filled with CO, NO2 observations commensurate 
    # to CASTNet observations in urban, rural, and suburban settings
    # for CO
    comm_co_u = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_co_u[:] = np.nan
    comm_co_r = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_co_r[:] = np.nan
    comm_co_su = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_co_su[:] = np.nan
    # for NO2
    comm_no2_u = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_no2_u[:] = np.nan
    comm_no2_r = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_no2_r[:] = np.nan
    comm_no2_su = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_no2_su[:] = np.nan
    # for O3
    comm_o3_u = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_o3_u[:] = np.nan
    comm_o3_r = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_o3_r[:] = np.nan
    comm_o3_su = np.empty([len(years), len(castnet_sites_fr), sos], 
                             dtype = float)
    comm_o3_su[:] = np.nan
    # create empty lists to be filled with latitude, longitude coordinates of 
    # AQS stations in urban, rural, and suburban environments
    # for CO
    coords_co_u = []
    coords_co_su = []
    coords_co_r = [] 
    # for NO2
    coords_no2_u = []
    coords_no2_su = []
    coords_no2_r = []
    # for O3
    coords_o3_u = []
    coords_o3_su = []
    coords_o3_r = []
    # loop through years in measuring period and select corresponding 
    # observations
    for year, counter1 in zip(years, np.arange(0, len(years), 1)): 
        # dates in year/months of interest
        idx = pd.date_range('%s-01-%s' %(sampling_months[0], year),                                         
                            '%s-31-%s' %(sampling_months[-1], year))
        CO_year = CO.loc[CO['Date GMT'].dt.year.isin([year])]      
        NO2_year = NO2.loc[NO2['Date GMT'].dt.year.isin([year])]
        O3_year = O3.loc[O3['Date GMT'].dt.year.isin([year])]
        # in a given year, loop through CASTNet sites in focus region and
        # determine if data exists  
        for site, castnet_data, counter2 in zip(castnet_sites_fr, 
                                                comm_castnet[counter1], 
                                                np.arange(0, len(castnet_sites_fr), 1)):
            # if data exists, the find AQS with the distance from the CASTNet site
            # defined by the search radius
            if np.where(np.isnan(castnet_data) == True)[0].shape[0] != sos:
                # find longitude/latitude of CASTNET station 
                csites_atsite = csites.loc[csites['SITE_ID'].isin([site])]
                lon_atsite = csites_atsite['LONGITUDE'].values[0]
                lat_atsite = csites_atsite['LATITUDE'].values[0]
                # define bounding box            
                top = lat_atsite + searchrad
                left = lon_atsite - searchrad
                bottom = lat_atsite - searchrad
                right = lon_atsite + searchrad
                # CO, NO2, and O3 observations in bounding box
                CO_atsite = CO_year[(CO_year['Latitude'] > bottom) & 
                                    (CO_year['Latitude'] <= top) &
                                    (CO_year['Longitude'] > left) &
                                    (CO_year['Longitude'] <= right)]  
                NO2_atsite = NO2_year[(NO2_year['Latitude'] > bottom) & 
                                      (NO2_year['Latitude'] <= top) &
                                      (NO2_year['Longitude'] > left) &
                                      (NO2_year['Longitude'] <= right)]    
                O3_atsite = O3_year[(O3_year['Latitude'] > bottom) & 
                                    (O3_year['Latitude'] <= top) &
                                    (O3_year['Longitude'] > left) &
                                    (O3_year['Longitude'] <= right)]                
                # loop through unique latitude and longitudes and find whether 
                # they are urban, rural, or suburban 
                # for CO                 
                uniqueCO = CO_atsite.groupby(['Latitude', 'Longitude'
                                              ]).size().reset_index(name = 'Freq')
                for index, row in uniqueCO.iterrows():
                    # find index in AQS site record corresponding to station's 
                    # observations
                    siting_lat = np.where(np.around(
                            aqs_sites['Latitude'].values, 4) == np.around(row['Latitude'], 4))
                    siting_lon = np.where(np.around(
                            aqs_sites['Longitude'].values, 4) == np.around(row['Longitude'], 4))
                    siting = np.intersect1d(siting_lat, siting_lon)
                    siting = aqs_sites.iloc[siting[0]]['Location Setting']
                    # CO at suburban locations                
                    if siting == 'SUBURBAN':
                        CO_su = CO_atsite.loc[CO_atsite['Latitude'].isin(
                            [row['Latitude']]) & CO_atsite['Longitude'].isin(
                            [row['Longitude']])]
                        CO_su = CO_su.groupby(['Date GMT']).mean()
                        CO_su.index = pd.DatetimeIndex(CO_su.index)
                        CO_su = CO_su.reindex(idx, 
                            fill_value = np.nan)['Sample Measurement']
                        comm_co_su[counter1, counter2] = CO_su
                        coords_co_su.append(tuple(row.values[:-1]))                
                    # CO at rural locations
                    if siting == 'RURAL':
                        CO_r = CO_atsite.loc[CO_atsite['Latitude'].isin(
                            [row['Latitude']]) & CO_atsite['Longitude'].isin(
                            [row['Longitude']])]   
                        CO_r = CO_r.groupby(['Date GMT']).mean()
                        CO_r.index = pd.DatetimeIndex(CO_r.index)
                        CO_r = CO_r.reindex(idx, 
                            fill_value = np.nan)['Sample Measurement']
                        comm_co_r[counter1, counter2] = CO_r
                        coords_co_r.append(tuple(row.values[:-1]))                    
                    # CO at urban locations
                    if siting == 'URBAN AND CENTER CITY':
                        CO_u = CO_atsite.loc[CO_atsite['Latitude'].isin(
                            [row['Latitude']]) & CO_atsite['Longitude'].isin(
                            [row['Longitude']])]                    
                        CO_u = CO_u.groupby(['Date GMT']).mean()
                        CO_u.index = pd.DatetimeIndex(CO_u.index)
                        CO_u = CO_u.reindex(idx, 
                            fill_value = np.nan)['Sample Measurement']
                        comm_co_u[counter1, counter2] = CO_u
                        coords_co_u.append(tuple(row.values[:-1]))                    
                # for NO2 
                uniqueNO2 = NO2_atsite.groupby(['Latitude', 'Longitude'
                                              ]).size().reset_index(name = 'Freq')
                for index, row in uniqueNO2.iterrows():
                    # find index in AQS site record corresponding to station's 
                    # observations
                    siting_lat = np.where(np.around(
                            aqs_sites['Latitude'].values, 4) == np.around(row['Latitude'], 4))
                    siting_lon = np.where(np.around(
                            aqs_sites['Longitude'].values, 4) == np.around(row['Longitude'], 4))
                    siting = np.intersect1d(siting_lat, siting_lon)
                    siting = aqs_sites.iloc[siting[0]]['Location Setting']
                    # NO2 at suburban locations                
                    if siting == 'SUBURBAN':
                        NO2_su = NO2_atsite.loc[NO2_atsite['Latitude'].isin(
                            [row['Latitude']]) & NO2_atsite['Longitude'].isin(
                            [row['Longitude']])]
                        NO2_su = NO2_su.groupby(['Date GMT']).mean()
                        NO2_su.index = pd.DatetimeIndex(NO2_su.index)
                        NO2_su = NO2_su.reindex(idx, 
                            fill_value = np.nan)['Sample Measurement']
                        comm_no2_su[counter1, counter2] = NO2_su
                        coords_no2_su.append(tuple(row.values[:-1]))                    
                    # NO2 at rural locations
                    if siting == 'RURAL':
                        NO2_r = NO2_atsite.loc[NO2_atsite['Latitude'].isin(
                            [row['Latitude']]) & NO2_atsite['Longitude'].isin(
                            [row['Longitude']])]   
                        NO2_r = NO2_r.groupby(['Date GMT']).mean()
                        NO2_r.index = pd.DatetimeIndex(NO2_r.index)
                        NO2_r = NO2_r.reindex(idx, 
                            fill_value = np.nan)['Sample Measurement']
                        comm_no2_r[counter1, counter2] = NO2_r
                        coords_no2_r.append(tuple(row.values[:-1]))                    
                    # NO2 at urban locations
                    if siting == 'URBAN AND CENTER CITY':
                        NO2_u = NO2_atsite.loc[NO2_atsite['Latitude'].isin(
                            [row['Latitude']]) & NO2_atsite['Longitude'].isin(
                            [row['Longitude']])]                    
                        NO2_u = NO2_u.groupby(['Date GMT']).mean()
                        NO2_u.index = pd.DatetimeIndex(NO2_u.index)
                        NO2_u = NO2_u.reindex(idx, 
                            fill_value = np.nan)['Sample Measurement']
                        comm_no2_u[counter1, counter2] = NO2_u
                        coords_no2_u.append(tuple(row.values[:-1]))                    
                # for O3 
                uniqueO3 = O3_atsite.groupby(['Latitude', 'Longitude'
                                              ]).size().reset_index(name = 'Freq')
                for index, row in uniqueO3.iterrows():
                    # find index in AQS site record corresponding to station's 
                    # observations
                    siting_lat = np.where(np.around(
                            aqs_sites['Latitude'].values, 3) == np.around(row['Latitude'], 3))
                    siting_lon = np.where(np.around(
                            aqs_sites['Longitude'].values, 3) == np.around(row['Longitude'], 3))
                    siting = np.intersect1d(siting_lat, siting_lon)
                    siting = aqs_sites.iloc[siting[0]]['Location Setting']
                    # O3 at suburban locations                
                    if siting == 'SUBURBAN':
                        O3_su = O3_atsite.loc[O3_atsite['Latitude'].isin(
                            [row['Latitude']]) & O3_atsite['Longitude'].isin(
                            [row['Longitude']])]
                        O3_su = O3_su.groupby(['Date GMT']).mean()
                        O3_su.index = pd.DatetimeIndex(O3_su.index)
                        O3_su = O3_su.reindex(idx, 
                            fill_value = np.nan)['Sample Measurement']
                        comm_o3_su[counter1, counter2] = O3_su
                        coords_o3_su.append(tuple(row.values[:-1]))                    
                    # O3 at rural locations
                    if siting == 'RURAL':
                        O3_r = O3_atsite.loc[O3_atsite['Latitude'].isin(
                            [row['Latitude']]) & O3_atsite['Longitude'].isin(
                            [row['Longitude']])]   
                        O3_r = O3_r.groupby(['Date GMT']).mean()
                        O3_r.index = pd.DatetimeIndex(O3_r.index)
                        O3_r = O3_r.reindex(idx, 
                            fill_value = np.nan)['Sample Measurement']
                        comm_o3_r[counter1, counter2] = O3_r
                        coords_o3_r.append(tuple(row.values[:-1]))                    
                    # O3 at urban locations
                    if siting == 'URBAN AND CENTER CITY':
                        O3_u = O3_atsite.loc[O3_atsite['Latitude'].isin(
                            [row['Latitude']]) & O3_atsite['Longitude'].isin(
                            [row['Longitude']])]                    
                        O3_u = O3_u.groupby(['Date GMT']).mean()
                        O3_u.index = pd.DatetimeIndex(O3_u.index)
                        O3_u = O3_u.reindex(idx, 
                            fill_value = np.nan)['Sample Measurement']
                        comm_o3_u[counter1, counter2] = O3_u
                        coords_o3_u.append(tuple(row.values[:-1]))                    
    return (comm_co_su, comm_co_r, comm_co_u, coords_co_su, coords_co_r, 
            coords_co_u, comm_no2_su, comm_no2_r, comm_no2_u, coords_no2_su, 
            coords_no2_r, coords_no2_u, comm_o3_su, comm_o3_r, comm_o3_u, 
            coords_o3_su, coords_o3_r, coords_o3_u)
# # # # # # # # # # # # #    
def open_gridded_idailyCTM(case, years):
    """GMI CTM output from idaily (gridded daily 12Z output) over the focus 
    region is opened and aggregated over the summers whose years are defined 
    in variable 'years.' Output reduced for focus region was created with 
    'yearly_idaily_combine.py' from HindcastMR2 output. 
    
    Parameters
    ----------   
    case : str
        Hindcast family (i.e. MR2, MR2-CCMI, FFIgac2, etc.)    
    years : list 
        Years of interest (n.b., as of 5 July 2018 files for 2008-2010 are 
        saved locally)
    
    Returns
    ----------
    lat : numpy.ndarrayn
        Latitude coordinates of focus region, units of degrees north, [lat,]
    lon : numpy.ndarray
        Longitude coordinates of focus region, units of degrees east, [lon,]
    pressure : numpy.ndarray
        Pressure coordinates of CTM, units of hPa, [pressure,]
    times : numpy.ndarray
        datetime.datetime objects corresponding to 12Z daily CTM output, 
        [time,]
    co : numpy.ndarray
        Gridded CTM CO output, units of volume mixing ratio, [time, pressure, 
        lat, lon]
    no : numpy.ndarray
        Gridded CTM NO output, units of volume mixing ratio, [time, pressure, 
        lat, lon]
    no2 : numpy.ndarray
        Gridded CTM NO2 output, units of volume mixing ratio, [time, pressure, 
        lat, lon]
    o3 : numpy.ndarray
        Gridded CTM O3 output, units of volume mixing ratio, [time, pressure, 
        lat, lon]
    """
    import numpy as np
    from datetime import datetime, timedelta
    from netCDF4 import Dataset
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # lists will be filled with CTM output from every summer
    co, no, no2, o3, times = [], [], [], [], []
    # open JJA file for each year in measuring period 
    for year in years: 
        # the latitudinal/longitudinal bounds of file are hard-coded; would have
        # to change if examining other regions
        infile = Dataset(pollutants_constants.PATH_GMI + '%s/' %(case) + 
                         'gmic_%s_%s_35N_275E_50N_295E_15.idaily.nc' 
                         %(case, year), 'r')
        # extract dimensional information only on the first iteration of year
        # loop 
        if year == years[0]:
            lat = infile.variables['lat'][:]
            lon = infile.variables['lon'][:]
            pressure = infile.variables['pressure'][:]
        # append 
        co.append(infile.variables['CO'][:])
        no.append(infile.variables['NO'][:])
        no2.append(infile.variables['NO2'][:])
        o3.append(infile.variables['O3'][:])
        # generate time dimension for netCDF file 
        def perdelta(start, end, delta):
            curr = start
            while curr < end:
                yield curr
                curr += delta
        time = []
        for result in perdelta(datetime(year, 6, 1, 12), 
                               datetime(year, 9, 1, 12), 
                               timedelta(hours = 24)):
            time.append(result)     
        # append time to multi-year list
        times.append(time)
        print('12Z gridded CTM output for %s loaded!' %year)   
    # convert lists to arrays
    times = np.hstack(np.vstack(times))
    co = np.vstack(co)
    no = np.vstack(no)
    no2 = np.vstack(no2)
    o3 = np.vstack(o3)
    return (lat, lon, pressure, times, co, no, no2, o3)  
# # # # # # # # # # # # #
def open_inst6_3d_ana_Np(years, lev):
    """function opens MERRA-2 inst6_3d_ana_Np output for the summers of 
    interest. Input files created with file 'dailynetCDF_combine.py.' 
    Geopotential height, temperature, U-/V-wind, and coordinate information are
    retrieved. 
    
    Parameters
    ----------      
    years : list 
        Years of interest
    lev : int/list
        Index or indices corresponding to the pressure level(s) of interest. 
        If a single level is desired, either an integer with a length of 1 
        can be passed to function. If multiple levels are desired, pass a 
        list with comma-separated indices. n.b., lev = 0 or lev = [0] would 
        correspond to the 1000 hPa level, lev = 1 corresponds to 975 hPa, 
        lev = 2 corresponds to 950, etc.              

    Returns
    ----------
    H : numpy.ndarray
        Geopotential height, units of m, [time, lev, lat, lon]
    T : numpy.ndarray 
        Air temperature, units of K, [time, lev, lat, lon]
    U : numpy.ndarray 
        U wind component, units of m/s, [time, lev, lat, lon]
    V : numpy.ndarray
        V wind component, units of m/s, [time, lev, lat, lon]    
    time : numpy.ndarray
        Datetime objects corresponding to each 6-hourly timestep of MERRA-2,
        [time,]
    lat : numpy.ndarray
        Latitude coordinates of focus region, units of degrees north, [lat,]
    lon : numpy.ndarray
        Longitude coordinates of focus region, units of degrees east, [lon,]
    pressure : 
        Pressure coordinates for level(s) of interest, units of hPa, [lev,]
    """
    import numpy as np
    from netCDF4 import Dataset
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/Desktop/')
    from generate_times import generate_times
    H, T, U, V, times = [], [], [], [], []
    # for each year in measuring period, open MERRA-2 meteorology and generate
    # times of MERRA-2 output
    for year in years: 
        times.append(generate_times(year, 6, 8, 6))
        # n.b., if focus region is changed, following line will need to be 
        # rewritten
        infile = Dataset(pollutants_constants.PATH_METEOROLOGY + 
            'MERRA2_300.inst6_3d_ana_Np_%d_35N_275E_50N_295E.nc' %year, 'r')
        # extract meteorology at level(s) of interest
        H.append(infile.variables['H'][:, lev])
        T.append(infile.variables['T'][:, lev])
        U.append(infile.variables['U'][:, lev])
        V.append(infile.variables['V'][:, lev])   
        # on first iteration of loop, extract coordinate information
        if year == years[0]:
            lat = infile.variables['lat'][:]
            lon = infile.variables['lon'][:]
            pressure = infile.variables['lev'][lev]
        print('6-hourly MERRA-2 meteorology for %s loaded!' %year)
    # stack along time dimension 
    H = np.vstack(H)
    T = np.vstack(T)
    U = np.vstack(U)
    V = np.vstack(V)
    times = np.hstack(times)        
    return H, T, U, V, times, lat, lon, pressure
# # # # # # # # # # # # #
def open_profile_singmonth(case, month, year, species, stations, levels):
    """for a given month, function retrieves a desired constituent at desired 
    stations and model pressure levels. Output constituent array represents 
    the volume mixing ratios of desired species at hourly temporal resolution.
    
    Parameters
    ----------   
    case : str
        GMI Simulation name
    month : str 
        Three month abbreviations, all lowercase (i.e. 'jun', 'may', etc.); 
        n.b. only June - August are available
    year : int
        Year of interest
    species : str
        Either 'CO', 'NO2', 'NO', or 'O3'
    stations : list
        'all' for all stations or a list of stations at which column diagnostic 
        information is desired
    levels : str/int
        'all' for full column (surface - 0.015 hPa), 0 for surface (i.e. 
        surface is 992.52405 hPa)
    times : numpy.ndarray
        Datetime objects for each model timestep, [step * no. days in months,]        

    Returns
    ----------
    col_latitude : numpy.ndarray
        Station latitudes, units awof degrees north, [no. stations,]
    col_longitude : numpy.ndarray
        Station longitudes, units of degrees east, [no. stations,]
    station_labels : numpy.ndarray
        Station name abbreviations, [no. stations,]
    pressure : numpy.ndarray 
        Pressure levels (centered), units of hPa, [no. pressure levels,]
    const : numpy.ndarray
        Constituent, units of volume mixing ratio, [no. stations, 24 * no. 
        days in month, no. pressure levels,]
    times : numpy.ndarray 
        Datetime objects for each model timestep, [step * no. days in months,]        
    """
    import numpy as np
    from netCDF4 import Dataset    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from generate_times import generate_times    
    # open Gmimod column diagnostic file of interest
    infile = Dataset(pollutants_constants.PATH_GMI + '%s/' %case +
                     'gmic_%s_%d_%s.profile.nc' %(case, year, month), 'r')
    # extract species names; change format from dtype='|S1'
    const_labels_raw = infile.variables['const_labels'][:]
    const_labels_raw = np.vstack(const_labels_raw)
    const_labels = []
    for row in np.arange(0, np.shape(const_labels_raw)[0]):
        temp = []
        for element in const_labels_raw[row]:
            temp.append(element.decode('UTF-8'))
        si = ''.join(temp)
        si = si.rstrip()
        const_labels.append(si[:])      
    # extract station names; change format from dtype='|S1'
    station_labels_raw = infile.variables['station_labels'][:]
    station_labels_raw = np.vstack(station_labels_raw)
    station_labels = []
    for row in np.arange(0, np.shape(station_labels_raw)[0]):
        temp = []
        for element in station_labels_raw[row]:
            temp.append(element.decode('UTF-8'))
        si = ''.join(temp)
        si = si.rstrip()
        station_labels.append(si[:])      
    del station_labels_raw, const_labels_raw    
    where_species = np.where(np.array(const_labels) == species)[0][0]
    # location position(s) of station(s) of interest
    if stations == 'all':
        where_stations = np.arange(0, len(station_labels), 1)
    else: 
        where_stations = np.where(np.in1d(station_labels, stations) == True)[0]
    # extract coordinate, pressure, and constituent information; n.b. 
    # const_surf is the same as const[:, :, :, 0]
    station_labels = np.array(station_labels)[where_stations]
    col_latitude = infile.variables['col_latitude'][where_stations]
    col_longitude = infile.variables['col_longitude'][where_stations]
    if levels == 'all':
        pressure = infile.variables['pressure'][:]
        const = infile.variables['const'][where_stations, :, where_species, :]
    else: 
        pressure = infile.variables['pressure'][levels]    
        const = infile.variables['const'][where_stations, :, 
                                where_species, levels]
    # generate times for hourly CTM output 
    mon_dict = {'jun' : 6, 'jul' : 7, 'aug' : 8}
    times = generate_times(year, mon_dict[month], mon_dict[month], 1)
    return (np.hstack(col_latitude), np.hstack(col_longitude), station_labels, 
            pressure, const, np.hstack(times))
# # # # # # # # # # # # #
def open_profile_multimonth(case, months, years, species, stations, 
                                 levels):     
    """similar to function open_profile_singmonth but for multiple months.
    
    Parameters
    ----------   
    case : str
        GMI Simulation name
    months : list
        Contains three month abbreviations, all lowercase (i.e. 'jun', 'may', 
        etc.); n.b. items must be in chronilogical order and only May -
        September are available
    years : int
        Years of interest        
    species : str
        Either 'CO', 'NO2', 'NO', or 'O3'
    stations : numpy.ndarray 
        'all' for all stations or a list of stations at which column diagnostic 
        information is desired
    levels : str
        'all' for full column (surface - 0.015 hPa), 0 for surface (i.e. 
        surface is 992.52405 hPa)

    Returns
    ----------
    col_latitude : numpy.ndarray
        Station latitudes, units of degrees north, [no. stations,]
    col_longitude: numpy.ndarray
        Station longitudes, units of degrees east, [no. stations,]
    station_labels : numpy.ndarray
        Station name abbreviations, [no. stations,]        
    pressure : numpy.ndarray 
        Pressure levels (centered), units of hPa, [no. pressure levels,]
    const_multiyear : numpy.ndarray
        Constituent, units of volume mixing ration, [no. stations, step * 
        no. days in months * no. years, no. pressure levels,]
    times_multiyear : numpy.ndarray 
        Datetime objects for each model timestep, 
        [step * no. days in months * no. years,]          
    """
    import numpy as np
    const_multiyear = []
    times_multiyear = []
    # loop through months
    for year in years:
        for month in months: 
            (col_latitude, col_longitude, station_labels, pressure, 
             const, times) = open_profile_singmonth(case, month, year, 
                species, stations, levels)
            # append constituent, times to list (all other information remains 
            # the same through iterations of the loop)
            const_multiyear.append(const)
            times_multiyear.append(times)            
    # stack along time dimension
    const_multiyear = np.hstack(const_multiyear)
    times_multiyear = np.hstack(times_multiyear)
    return (col_latitude, col_longitude, station_labels, pressure, 
            const_multiyear, times_multiyear)
# # # # # # # # # # # # #