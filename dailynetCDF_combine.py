#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAME
    dailynetCDF_combine.py
PURPOSE
    Function extracts relevant output from GMI idaily (12Z) files and MERRA
    daily files, aggregates over a summer (JJA), and outputs to a single file 
PROGRAMMER
    Gaige Hunter Kerr
REVISION HISTORY
    05072018 -- initial version created
    08072018 -- function 'yearly_idaily_combine' modified to include CTM 
                simulation name as an argument.
    18072018 -- filename changed to 'dailynetCDF_combine.py' and function 
                'tavg1_2d_aer_Nx_combine' added. 
"""
# # # # # # # # # # # # #
def yearly_idaily_combine(ulat, llat, llon, rlon, presslevel, case, year):
    """using Gmimod freq1 species concentration files for summer (JJA) months, 
    function uses latitude/longitude inputs to find the closest values in 
    the CTM's grid. Output is reduced to this smaller grid and daily gridded 
    CO, NO, NO2, and O3 are extracted along with dimensional information 
    (i.e., latitude, longitude, pressure). Daily (12Z) JJA output is saved in 
    single file. 
    
    Parameters
    ---------- 
    ulat : int/float
        Upper latitude bound for focus region
    llat : int/float
        Lower latitude bound for focus region    
    llon : int/float
        Left longitude bound for focus region    
    rlon : int/float
        Right longitude bound for focus region            
    presslevel : int
        The number of pressure levels (started from the lowest model level, 
        992.52 hPa) that will be saved to output file.
    case : str
        CTM simulation name
    year : int 
        Year of interest

    Returns
    ----------
    None        
    """
    import numpy as np
    from netCDF4 import Dataset
    # path to files 
    path = '/mnt/scratch4/gkerr4/data/%s/%d/' %(case, year)
    # summer month abbreviations
    months = ['jun', 'jul', 'aug']
    # list creation; will be filled with CTM output in region for each month
    co, no, no2, o3 = [], [], [], []
    # open input file for each summer month in year 
    for month in months:
        infile = Dataset(path + 'gmic_%s_%d_%s.idaily.nc' %(case, year, month), 
                         'r')
        if month == months[0]:
            # find indices corresponding to region of interest
            llat = np.abs(infile.variables['latitude_dim'][:] - llat).argmin()
            ulat = np.abs(infile.variables['latitude_dim'][:] - ulat).argmin()
            llon = np.abs(infile.variables['longitude_dim'][:] - llon).argmin()
            rlon = np.abs(infile.variables['longitude_dim'][:] - rlon).argmin()
            # freq1 species name
            freq1_labels = infile.variables['freq1_labels'][:]
            # extract species names to more usable format
            cs_gmi = []
            for row in np.arange(0, np.shape(freq1_labels)[0]):
                # concatenate chemical species 
                si = freq1_labels[row, :].tostring()            
                si = si.decode('utf-8') 
                # remove whitespace after chemical species
                si = si.rstrip()
                # write chemical species to cs_gmi list
                cs_gmi.append(si[:])
            # find latitude and longitude of focus region
            lat = infile.variables['latitude_dim'][llat:ulat+1]
            lon = infile.variables['longitude_dim'][llon:rlon+1]
            # eta dimension of interest
            eta = infile.variables['eta_dim'][:presslevel]
        # species of interest include CO, NO, NO2, and O3
        co_mon = infile.variables['const_freq1'][
                :, np.where(np.array(cs_gmi) == 'CO')[0][0], 
                :presslevel, llat:ulat+1, llon:rlon+1]
        no_mon = infile.variables['const_freq1'][
                :, np.where(np.array(cs_gmi) == 'NO')[0][0], 
                :presslevel, llat:ulat+1, llon:rlon+1]
        no2_mon = infile.variables['const_freq1'][
                :, np.where(np.array(cs_gmi) == 'NO2')[0][0], 
                :presslevel, llat:ulat+1, llon:rlon+1]
        o3_mon = infile.variables['const_freq1'][
                :, np.where(np.array(cs_gmi) == 'O3')[0][0], 
                :presslevel, llat:ulat+1, llon:rlon+1]
        # append monthly arrays to list comprised of multi-month list
        co.append(co_mon)
        no.append(no_mon)
        no2.append(no2_mon)
        o3.append(o3_mon)
    # create output file; naming convention is the same as input file, 
    # with trailing values corresponding to: lower latitude, left longitude, 
    # upper latitude, right longitude, and the number of pressure levels 
    # sought
    outfile = Dataset(path + 'gmic_%s_%d_%dN_%dE_%dN_%dE_%d.idaily.nc'
                      %(case, year, lat[0], lon[0], lat[-1], lon[-1], presslevel), 'w')    
    # define set of dimensions
    lat_dim = outfile.createDimension('lat', len(lat))
    lon_dim = outfile.createDimension('lon', len(lon))
    pressure_dim = outfile.createDimension('pressure', len(eta))
    time_dim = outfile.createDimension('time', None)
    # create coordinate variables for variables
    co_dataset = outfile.createVariable('CO', np.float32, 
                                        ('time', 'pressure', 'lat', 'lon')) 
    no_dataset = outfile.createVariable('NO', np.float32, 
                                        ('time', 'pressure', 'lat', 'lon')) 
    no2_dataset = outfile.createVariable('NO2', np.float32, 
                                         ('time', 'pressure', 'lat', 'lon')) 
    o3_dataset = outfile.createVariable('O3', np.float32, 
                                        ('time', 'pressure', 'lat', 'lon')) 
    lat_dim = outfile.createVariable('lat', np.float32, ('lat'))
    lon_dim = outfile.createVariable('lon', np.float32, ('lon'))
    pressure_dim = outfile.createVariable('pressure', np.float32, ('pressure'))
    # add data and define their attributes
    # latitude
    lat_dim[:] = lat[:]
    lat_dim.standard_name = 'latitude'
    lat_dim.units = 'degrees_north'
    # longitude
    lon_dim[:] = lon
    lon_dim.standard_name = 'longitude'
    lon_dim.units = 'degree_east'
    # pressure 
    pressure_dim[:] = eta
    pressure_dim.standard_name = 'air_pressure'
    pressure_dim.units = 'hPa'
    # CO
    co_dataset[:] = np.vstack(co)
    co_dataset.standard_name = 'mass_fraction_of_CO_in_air'
    co_dataset.units = 'volume mixing ratio'
    # NO
    no_dataset[:] = np.vstack(no)
    no_dataset.standard_name = 'mass_fraction_of_NO_in_air'
    no_dataset.units = 'volume mixing ratio'
    # NO2
    no2_dataset[:] = np.vstack(no2)
    no2_dataset.standard_name = 'mass_fraction_of_NO2_in_air'
    no2_dataset.units = 'volume mixing ratio'
    # O3
    o3_dataset[:] = np.vstack(o3)
    o3_dataset.standard_name = 'mass_fraction_of_O3_in_air'
    o3_dataset.units = 'volume mixing ratio'
    # global attributes
    outfile.title = 'Gmimod freq1 species concentration file for region of ' + \
        'interest, adapted from monthly files into a single summer (JJA) files'
    outfile.author = 'Gaige Hunter Kerr'
    outfile.comment = 'GMI CTM output in this file is daily (12Z) NO, NO2, ' + \
        'CO, and O3 for the first %d pressure levels. ' %(presslevel) + \
        'Output latitudinally bound by %d-%dN ' %(lat[0], lat[-1]) +\
        'and longitudinally bound by %d-%dE is included.' %(lon[0], lon[-1])
    outfile.close()   
    return 
# # # # # # # # # # # # #
def tavg1_2d_aer_Nx_combine(ulat, llat, llon, rlon):
    """function opens all daily files of 1-hourly, time-averaged, single-level,
    assimilation, aerosol diagnostics, extracts output spatially bound by the 
    input parameters and creates a single output file for the summer (MJJAS) of
    2006 for comparison with daily-varying dust and aerosols GMI simulation.

    Parameters
    ---------- 
    ulat : int/float
        Upper latitude bound for focus region
    llat : int/float
        Lower latitude bound for focus region    
    llon : int/float
        Left longitude bound for focus region    
    rlon : int/float
        Right longitude bound for focus region            

    Returns
    ----------
    None   
    """
    import numpy as np
    import glob
    from netCDF4 import Dataset
    # path to files 
    path = '/mnt/scratch3/gkerr4/data/MERRA-2/tavg1_2d_aer_Nx/'
    # file names with wildcard character for date
    mfstring = 'MERRA2_300.tavg1_2d_aer_Nx.*.SUB.nc4'
    infiles = []
    for file in glob.glob(path + mfstring):
        infiles.append(file)
    # sort input files by date (YYYYMMDD format)
    infiles.sort()
    # create lists, will be filled with hourly values over the NEUS
    dus, sea, bla, org, so4 = [], [], [], [], []
    for file in infiles:     
        infile = Dataset(file, 'r')
        if file == infiles[0]:
            # convert longitude to (0 - 360)        
            lon = infile.variables['lon'][:] + 360.
            lat = infile.variables['lat'][:]
            # find indices corresponding to region of interest on first iteration
            # through loop          
            llat = np.abs(lat - llat).argmin()
            ulat = np.abs(lat - ulat).argmin()
            llon = np.abs(lon - llon).argmin()
            rlon = np.abs(lon - rlon).argmin()
            lon = lon[llon:rlon+1]
            lat = lat[llat:ulat+1]    
        # S in variable names corresponds to the surface rather than column
        # Dust Surface Mass Concentration, units of kg m-3    
        dus_day = infile.variables['DUSMASS'][:, llat:ulat+1, llon:rlon+1]
        # Sea Salt Surface Mass Concentration, units of kg m-3    
        sea_day = infile.variables['SSSMASS'][:, llat:ulat+1, llon:rlon+1]
        # Black Carbon Surface Mass Concentration, units of kg m-3    
        bla_day = infile.variables['BCSMASS'][:, llat:ulat+1, llon:rlon+1]
        # Organic Carbon Surface Mass Concentration __ENSEMBLE__, units of kg m-3    
        org_day = infile.variables['OCSMASS'][:, llat:ulat+1, llon:rlon+1]
        # SO4 Surface Mass Concentration __ENSEMBLE__, units of kg m-3
        so4_day = infile.variables['SO4SMASS'][:, llat:ulat+1, llon:rlon+1]
        # append hourly daily output to multi-day lists
        dus.append(dus_day)
        sea.append(sea_day) 
        bla.append(bla_day) 
        org.append(org_day) 
        so4.append(so4_day)
    # create output file; naming convention is the same as input file, 
    # with trailing values corresponding to: lower latitude, left longitude, 
    # upper latitude, right longitude
    outfile = Dataset(path + 
                      'MERRA2_300.tavg1_2d_aer_Nx_2006_%dN_%dE_%dN_%dE.nc'
                      %(lat[0], lon[0], lat[-1], lon[-1]), 'w')    
    # define set of dimensions
    lat_dim = outfile.createDimension('lat', len(lat))
    lon_dim = outfile.createDimension('lon', len(lon))
    time_dim = outfile.createDimension('time', None)
    # create coordinate variables for variables
    dus_dataset = outfile.createVariable('DUST', np.float32, 
                                        ('time', 'lat', 'lon')) 
    sea_dataset = outfile.createVariable('SS', np.float32, 
                                        ('time', 'lat', 'lon')) 
    bla_dataset = outfile.createVariable('BC', np.float32, 
                                         ('time', 'lat', 'lon')) 
    org_dataset = outfile.createVariable('OC', np.float32, 
                                         ('time', 'lat', 'lon')) 
    so4_dataset = outfile.createVariable('SO4', np.float32, 
                                         ('time', 'lat', 'lon')) 
    lat_dim = outfile.createVariable('lat', np.float32, ('lat'))
    lon_dim = outfile.createVariable('lon', np.float32, ('lon'))
    # add data and define their attributes
    # latitude
    lat_dim[:] = lat[:]
    lat_dim.standard_name = 'latitude'
    lat_dim.units = 'degrees_north'
    # longitude
    lon_dim[:] = lon
    lon_dim.standard_name = 'longitude'
    lon_dim.units = 'degree_east'
    # Dust
    dus_dataset[:] = np.vstack(dus)
    dus_dataset.standard_name = 'Dust Surface Mass Concentration'
    dus_dataset.long_name = 'Dust Surface Mass Concentration'
    dus_dataset.units = 'kg m-3'
    # Sea salt
    sea_dataset[:] = np.vstack(sea)
    sea_dataset.standard_name = 'Sea Salt Surface Mass Concentration'
    sea_dataset.long_name = 'Sea Salt Surface Mass Concentration'
    sea_dataset.units = 'kg m-3'
    # Black carbon
    bla_dataset[:] = np.vstack(bla)
    bla_dataset.standard_name = 'Black Carbon Surface Mass Concentration'
    bla_dataset.long_name = 'Black Carbon Surface Mass Concentration'
    bla_dataset.units = 'kg m-3'
    # Organic carbon
    org_dataset[:] = np.vstack(org)
    org_dataset.standard_name = 'Organic Carbon Surface Mass Concentration __ENSEMBLE__'
    org_dataset.long_name = 'Organic Carbon Surface Mass Concentration __ENSEMBLE__'
    org_dataset.units = 'kg m-3'
    # SO4
    so4_dataset[:] = np.vstack(so4)
    so4_dataset.standard_name = 'SO4 Surface Mass Concentration __ENSEMBLE__'
    so4_dataset.long_name = 'SO4 Surface Mass Concentration __ENSEMBLE__'
    so4_dataset.units = 'kg m-3'
    outfile.title = 'MERRA2 tavg1_2d_aer_Nx: 2d,1-Hourly,Time-averaged,' + \
    'Single-Level,Assimilation,Aerosol Diagnostics for summer (MJJAS) ' + \
    'of 2006, adapted from daily files.'
    outfile.author = 'Gaige Hunter Kerr'
    outfile.comment = 'Output latitudinally bound by %d-%dN ' %(lat[0], lat[-1]) +\
        'and longitudinally bound by %d-%dE is included.' %(lon[0], lon[-1])
    outfile.close()
    return 
# # # # # # # # # # # # #    
def inst6_3d_ana_Np_combine(year, ulat, llat, llon, rlon):
    """function opens all daily files of 6-Hourly, Instantaneous, 
    Pressure-Level, Analysis, Analyzed Meteorological Fields for a particular
    summer (JJA) and extracts grid cells bounded by the input latitude and 
    longitudes and creates a single output file for the summer. 

    Parameters
    ----------
    year : int
        Year of interest
    ulat : int/float
        Upper latitude bound for focus region
    llat : int/float
        Lower latitude bound for focus region    
    llon : int/float
        Left longitude bound for focus region    
    rlon : int/float
        Right longitude bound for focus region            

    Returns
    ----------
    None    
    """
    import numpy as np
    from netCDF4 import Dataset
    import glob
    # path to MERRA-2 data holdings
    path = '/mnt/scratch3/gkerr4/data/MERRA-2/inst6_3d_ana_Np/'
    # file names with wildcard character for date
    mfstring = 'MERRA2_300.inst6_3d_ana_Np.%s*.SUB.nc4' %year
    infiles = []
    for file in glob.glob(path + mfstring):
        infiles.append(file)
    # sort input files by date (YYYYMMDD format)
    infiles.sort()
    H, T, U, V = [], [], [], []
    for file in infiles:     
        infile = Dataset(file, 'r')
        if file == infiles[0]:
            # convert longitude to (0 - 360)        
            lon = infile.variables['lon'][:] + 360.
            lat = infile.variables['lat'][:]
            # find indices corresponding to region of interest on first iteration
            # through loop          
            llat = np.abs(lat - llat).argmin()
            ulat = np.abs(lat - ulat).argmin()
            llon = np.abs(lon - llon).argmin()
            rlon = np.abs(lon - rlon).argmin()
            lon = lon[llon:rlon+1]
            lat = lat[llat:ulat+1]    
            pressure = infile.variables['lev'][:]
        # extract relevant variables and append 4-hourly fields to JJA 
        # list
        # Geopotential height
        H.append(infile.variables['H'][:, :, llat:ulat+1, llon:rlon+1])
        # Air temperature
        T.append(infile.variables['T'][:, :, llat:ulat+1, llon:rlon+1])
        # Eastward wind component
        U.append(infile.variables['U'][:, :, llat:ulat+1, llon:rlon+1])
        # Northward wind component
        V.append(infile.variables['V'][:, :, llat:ulat+1, llon:rlon+1])
    # create output file; naming convention is the same as input file, 
    # with trailing values corresponding to: lower latitude, left longitude, 
    # upper latitude, right longitude
    outfile = Dataset(path + 
                      'MERRA2_300.inst6_3d_ana_Np_%d_%dN_%dE_%dN_%dE.nc'
                      %(year, lat[0], lon[0], lat[-1], lon[-1]), 'w')    
    # define set of dimensions
    lat_dim = outfile.createDimension('lat', len(lat))
    lon_dim = outfile.createDimension('lon', len(lon))
    pressure_dim = outfile.createDimension('lev', len(pressure))
    time_dim = outfile.createDimension('time', None)
    # create coordinate variables for variables
    H_dataset = outfile.createVariable('H', np.float32, 
                                       ('time', 'lev', 'lat', 'lon')) 
    T_dataset = outfile.createVariable('T', np.float32, 
                                        ('time', 'lev', 'lat', 'lon')) 
    U_dataset = outfile.createVariable('U', np.float32, 
                                         ('time', 'lev', 'lat', 'lon')) 
    V_dataset = outfile.createVariable('V', np.float32, 
                                         ('time', 'lev', 'lat', 'lon')) 
    lat_dim = outfile.createVariable('lat', np.float32, ('lat'))
    lon_dim = outfile.createVariable('lon', np.float32, ('lon'))
    pressure_dim = outfile.createVariable('lev', np.float32, ('lev'))
    # add data and define their attributes
    # latitude
    lat_dim[:] = lat[:]
    lat_dim.standard_name = 'latitude'
    lat_dim.units = 'degrees_north'
    # longitude
    lon_dim[:] = lon
    lon_dim.standard_name = 'longitude'
    lon_dim.units = 'degree_east'
    # pressure
    pressure_dim[:] = pressure
    pressure_dim.standard_name = 'air_pressure'
    pressure_dim.long_name = 'vertical level'
    pressure_dim.units = 'hPa'
    # H
    H_dataset[:] = np.vstack(H)
    H_dataset.standard_name = 'geopotential_height'
    H_dataset.long_name = 'Geopotential height'
    H_dataset.units = 'm'
    # T
    T_dataset[:] = np.vstack(T)
    T_dataset.standard_name = 'air_temperature'
    T_dataset.long_name = 'Air temperature'
    T_dataset.units = 'K'
    # U 
    U_dataset[:] = np.vstack(U)
    U_dataset.standard_name = 'eastward_wind'
    U_dataset.long_name = 'Eastward wind component'
    U_dataset.units = 'm/s'
    # V
    V_dataset[:] = np.vstack(V)
    V_dataset.standard_name = 'northward_wind'
    V_dataset.long_name = 'Northward wind component'
    V_dataset.units = 'm/s'
    # attributes
    outfile.title = 'MERRA2_300.inst6_3d_ana_Np, 3d, MERRA2 6-Hourly, ' + \
         'Instantaneous, Pressure-Level, Analysis, Analyzed Meteorological ' + \
         'Fields, for summer (JJA) %s.' %year
    outfile.author = 'Gaige Hunter Kerr'
    outfile.comment = 'Output latitudinally bound by %d-%dN ' %(lat[0], lat[-1]) + \
        'and longitudinally bound by %d-%dE is included.' %(lon[0], lon[-1])
    outfile.close()
    return 
# # # # # # # # # # # # # #    
# define function inputs 
ulat = 50.
llat = 35.
llon = 275.
rlon = 295.
presslevel = 15
# for HindcastMR2-DiurnalAvgT
yearly_idaily_combine(ulat, llat, llon, rlon, presslevel, 
    'HindcastMR2-DiurnalAvgT', 2008)
# for tavg1_2d_aer_Nx
tavg1_2d_aer_Nx_combine(ulat, llat, llon, rlon)
# for inst6_3d_ana_Np
year = 2008
inst6_3d_ana_Np_combine(year, ulat, llat, llon, rlon)