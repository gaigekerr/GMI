#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAME
    yearly_idaily_combinea.py
PURPOSE
    Function extracts relevant output from GMI idaily (12Z) files and 
    aggregates over a summer (JJA) to output a single file
PROGRAMMER
    Gaige Hunter Kerr
REVISION HISTORY
    05072018 -- initial version created
    08072018 -- function 'yearly_idaily_combine' modified to include CTM 
                simulation name as an argument.
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
    path = '/mnt/scratch4/gkerr4/data/HindcastMR2/%d/' %year
    # summer month abbreviations
    months = ['jun', 'jul', 'aug']
    # list creation; will be filled with CTM output in region for each month
    co, no, no2, o3 = [], [], [], []
    # open input file for each summer month in year 
    for month in months:
        infile = Dataset(path + 'gmic_HindcastMR2_%d_%s.idaily.nc' %(year, month), 'r')
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
    outfile = Dataset(path + 'gmic_HindcastMR2_%d_%dN_%dE_%dN_%dE_%d.idaily.nc'
                      %(year, lat[0], lon[0], lat[-1], lon[-1], presslevel), 'w')    
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
# define function inputs 
ulat = 50.
llat = 35.
llon = 275.
rlon = 285.
presslevel = 15
case = 'HindcastMR2-DiurnalAvgT'
year = 2008
yearly_idaily_combine(ulat, llat, llon, rlon, presslevel, case, year)
