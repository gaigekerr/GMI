#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAME
    separateFiles_singleyear_combine.py
PURPOSE
    Combine GMI CTM hindcast output into smaller files (restrict model leves 
    and spatial domain).
PROGRAMMER
    Gaige Hunter Kerr
REVISION HISTORY
    20173108 -- initial version created
    20170109 -- moved to eps-cassini     
    20172511 -- function 'separateFiles_singleyear_combine' modified such that 
                it not only extracts O3 from input files but also CO, NO, and 
                NO2 at GMI/CASTNET sites
    20180501 -- changed path names to accomodate perturbed emissions runs                
    20181809 -- function 'monthly_overpass2_combine' added to combine monthly 
                GMI overpass2 files
"""
def separateFiles_singleyear_combine(case, year, sites): 
    """function opens monthly column diagnostic files at GMI/CASTNET sites 
    specified in parameter 'sites' for year specified in parameter 'year'. 
    Column O3, NO, NO2, and CO is extracted along with site's lat/lon. 
    Information from all sites are used to create a yearly netCDF file for the 
    specified case.
    
    Parameters
    ----------    
    case : str
        Hindcast family (i.e. MR2, MR2-CCMI, FFIgac2, etc.)
    year : int
        unfiltered EPA AQS Data Mart O3 observations
    sites : list
        GMI site codes
    
    Returns
    ----------    
    None        
    """
    import numpy as np
    import glob
    from netCDF4 import MFDataset, Dataset
    # path to data on EPS-Cassini or local
    path = '/Volumes/GAIGEKERR/GMI/egu_t_sensitivity/%d/' %(year)
    #path = '/home/local/WIN/gkerr4/data/%s/%d/' %(case, year)
    # relevant parameters for GMI/CASTNET sites
    site_multisites = []
    o3_multisites, co_multisites, no_multisites, no2_multisites = [], [], [], []
    lat_multisites = []
    lon_multisites = []
    pressure_multisites = []
    # month abbreviations used for file naming convenctions
    months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 
              'sep', 'oct', 'nov', 'dec']
    # loop through GMI/CASTNET sites
    i = 1 
    for site in sites:
        print('extracting data from %s!' %(site))
        # for sensitivity simulations
        mfstring = 'gmic_GHKerr-DailyEmiss%s_%d*_%s.profile.nc' %(str(year)[2:], year, site)
        # for Hindcast simulations
        #mfstring = 'gmic_%s_%d_*_%s.profile.nc' %(case, year, site)
        # open Gmimod column diagnostic file for a particular site 
        infiles = []
        for file in glob.glob(path + mfstring):
            infiles.append(file)
        # sort list according given text of months; for more information see
        # https://stackoverflow.com/questions/33767899/sort-list-according-given-text-of-months
        infiles.sort()
        #infiles.sort(key = lambda x: months.index(x.split('_')[3]))
        # MFDataset aggregates along the slowest varying dimension, opens files 
        # and creates an instance of the netCDF4 class
        infiles = MFDataset(infiles, 'r')
        # extract station information (lat, lon converted from (0 to 360) to (-180 
        # to 180), and station name)
        lat = infiles.variables['col_latitude'][:][0]
        lon = (((infiles.variables['col_longitude'][:] + 180.) % 360.) - 180.)[0]
        # fetch Species name
        const_labels = infiles.variables['const_labels'][:]
        # convert chemical species listed in freq1_labels into more user-friendly
        # format contained in list cs_gmi
        cs_gmi = []
        for row in np.arange(0, np.shape(const_labels)[0]):
            # concatenate chemical species 
            si = const_labels[row, :].tostring()            
            si = si.decode('utf-8') 
            # remove whitespace after chemical species
            si = si.rstrip()
            # write chemical species to cs_gmi array
            cs_gmi.append(si[:])
        # position of chemical species
        o3 = np.where(np.array(cs_gmi) == 'O3')[0][0]
        co = np.where(np.array(cs_gmi) == 'CO')[0][0]
        no = np.where(np.array(cs_gmi) == 'NO')[0][0]
        no2 = np.where(np.array(cs_gmi) == 'NO2')[0][0]
        # fetch O3 concentrations at station from infile variable 
        # const (Constituent at const [volume mixing ratio])
        # dimensions are (hours in year, species, eta)
        o3 = infiles.variables['const'][:, o3]
        co = infiles.variables['const'][:, co]
        no = infiles.variables['const'][:, no]
        no2 = infiles.variables['const'][:, no2]
        # append information from individual GMI (CASTNET) sites to lists for
        # multiple sites
        site_multisites.append(site)
        lat_multisites.append(lat)
        lon_multisites.append(lon)
        o3_multisites.append(o3)
        co_multisites.append(co)
        no_multisites.append(no)
        no2_multisites.append(no2)
        # only fetch pressure at first iteration 
        if i == 1: 
            pressure_multisites.append(infiles.variables['pressure'][:])
            i = i + 1
    # ensure that data was correctly opened and parsed
    if sites != site_multisites:
        print('site order does not match!')
    # dimensions of final ozone data extracted for collection of sites should 
    # be [number of sites, hours in year, eta levels]
    outfile = Dataset(path + '%s_%d.nc' %(case, year), 'w', 
                      format = 'NETCDF4')
    # define set of dimensions
    sites_dataset = outfile.createDimension('sites', len(site_multisites))
    pressure_dataset = outfile.createDimension('pressure', len(pressure_multisites[0]))
    time_dataset = outfile.createDimension('time', None)
    # create coordinate variables for 3-dimensions
    sites = outfile.createVariable('sites', np.str, ('sites',)) 
    pressure = outfile.createVariable('pressure', np.float32, ('pressure',)) 
    # create the data variables
    o3_dataset = outfile.createVariable('o3', np.float32, 
                                        ('sites', 'time', 'pressure'))
    co_dataset = outfile.createVariable('co', np.float32, 
                                        ('sites', 'time', 'pressure'))
    no_dataset = outfile.createVariable('no', np.float32, 
                                        ('sites', 'time', 'pressure'))
    no2_dataset = outfile.createVariable('no2', np.float32, 
                                        ('sites', 'time', 'pressure'))    
    lat_dataset = outfile.createVariable('lat', np.float32, ('sites'))
    lon_dataset = outfile.createVariable('lon', np.float32, ('sites'))
    # grow data along site, unlimited, and pressure dimensions
    o3_dataset[:] = np.stack(o3_multisites)
    co_dataset[:] = np.stack(co_multisites)
    no_dataset[:] = np.stack(no_multisites)
    no2_dataset[:] = np.stack(no2_multisites)
    lat_dataset[:] = np.stack(lat_multisites)
    lon_dataset[:] = np.stack(lon_multisites)
    sites[:] = np.stack(site_multisites)
    pressure[:] = np.stack(pressure_multisites)
    # variable attributes 
    o3_dataset.standard_name = 'mass_fraction_of_ozone_in_air'
    o3_dataset.units = '1'
    co_dataset.standard_name = 'mass_fraction_of_CO_in_air'
    co_dataset.units = '1'
    no_dataset.standard_name = 'mass_fraction_of_NO_in_air'
    no_dataset.units = '1'
    no2_dataset.standard_name = 'mass_fraction_of_NO2_in_air'
    no2_dataset.units = '1'
    pressure.standard_name = 'air_pressure'
    pressure.units = 'hPa'
    lat_dataset.standard_name = 'latitude'
    lat_dataset.units = 'degrees_north'
    lon_dataset.standard_name = 'longitude'
    lon_dataset.units = 'degree_east'
    # global attributes
    #title — a succinct description of the data set
    #insititute — the organisation where the data were produced
    #source — how the data were produced, e.g. model type, run number and circumstances
    #history — an audit trail of data set processing
    #references — a list of references that describe the data or the methodology used
    #comment — other useful information not covered elsewhere that adds value
    #author — the person(s) who generated the data
    outfile.title = 'Gmimod column diagnostic file for all GMI sites' 
    outfile.author = 'Gaige Hunter Kerr'
    outfile.close()
    return
# # # # # # # # # # # # #
def monthly_overpass2_combine(year, case, ulat, llat, llon, rlon):
    """function opens monthly GMI overpass2 files, finds surface-level 
    model output over the specified latitude/longitude bounds, and creates 
    a single output file for JJA of a particular year

    Parameters
    ----------
    year : int
        Year of interest
    case : str
        Hindcast family (i.e. MR2, MR2-CCMI, FFIgac2, etc.)
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
    # month abbreviations
    months = ['jan.overpass2.nc', 'feb.overpass2.nc', 'mar.overpass2.nc', 
              'apr.overpass2.nc', 'may.overpass2.nc', 'jun.overpass2.nc', 
              'jul.overpass2.nc', 'aug.overpass2.nc', 'sep.overpass2.nc', 
              'oct.overpass2.nc', 'nov.overpass2.nc', 'dec.overpass2.nc']
    # path to MERRA-2 data holdings
    path = '/mnt/scratch1/gaige/data/GMI/overpass2/%s/%s/' %(case, year)
    # file names with wildcard character for date
    mfstring = 'gmic_%s_%s_*.overpass2.nc' %(case, year)
    infiles = []
    for file in glob.glob(path + mfstring):
        infiles.append(file)
    # sort input files by month (n.b. from 
    # https://stackoverflow.com/questions/33767899/sort-list-according-given-text-of-months)
    infiles.sort(key = lambda x: months.index(x.split('_')[3]))
    gridBoxHeight_overpass, cloudFraction_overpass = [], []
    o3_all, no_all, no2_all, co_all = [], [], [], []
    for file in infiles:     
        infile = Dataset(file, 'r')
        if file == infiles[0]:
            # convert longitude to (0 - 360)        
            lon = infile.variables['longitude_dim'][:]
            lat = infile.variables['latitude_dim'][:]
            begTime_overpass2 = infile.variables['begTime_overpass2'][:]
            endTime_overpass2 = infile.variables['endTime_overpass2'][:]
            # find indices corresponding to region of interest on first iteration
            # through loop          
            llat = np.abs(lat - llat).argmin()
            ulat = np.abs(lat - ulat).argmin()
            llon = np.abs(lon - llon).argmin()
            rlon = np.abs(lon - rlon).argmin()
            lon = lon[llon:rlon+1]
            lat = lat[llat:ulat+1]    
        # convert chemical species listed in freq1_labels into more user-friendly
        # format contained in list cs_gmi
        overpass_labels = infile.variables['overpass_labels'][:]
        cs_gmi = []
        for row in np.arange(0, np.shape(overpass_labels)[0]):
            # concatenate chemical species 
            si = overpass_labels[row, :].tostring()            
            si = si.decode('utf-8') 
            # remove whitespace after chemical species
            si = si.rstrip()
            # write chemical species to cs_gmi array
            cs_gmi.append(si[:])
        # position of chemical species
        o3 = np.where(np.array(cs_gmi) == 'O3')[0][0]
        co = np.where(np.array(cs_gmi) == 'CO')[0][0]
        no = np.where(np.array(cs_gmi) == 'NO')[0][0]
        no2 = np.where(np.array(cs_gmi) == 'NO2')[0][0]
        # extract relevant variables and append daily fields to JJA list
        gridBoxHeight_overpass.append(infile.variables[
            'gridBoxHeight_overpass'][:, 0, llat:ulat+1, llon:rlon+1])
        cloudFraction_overpass.append(infile.variables[
            'cloudFraction_overpass'][:, 0, llat:ulat+1, llon:rlon+1])
        o3_all.append(infile.variables['const_overpass']
            [:, o3, 0, llat:ulat+1, llon:rlon+1])
        no_all.append(infile.variables['const_overpass']
            [:, no, 0, llat:ulat+1, llon:rlon+1])
        no2_all.append(infile.variables['const_overpass']
            [:, no2, 0, llat:ulat+1, llon:rlon+1])
        co_all.append(infile.variables['const_overpass']
            [:, co, 0, llat:ulat+1, llon:rlon+1])
    # create JJA output file
    outfile = Dataset(path + 'gmic_%s_%s_jja_%dN_%dE_%dN_%dE.overpass2.nc'
                      %(case, year, lat[0], lon[0], lat[-1], lon[-1]), 'w', 
                      format = 'NETCDF4')
    # define set of dimensions
    longitude_dim = outfile.createDimension('longitude_dim', len(lon))
    latitude_dim = outfile.createDimension('latitude_dim', len(lat))
    scalar_dim = outfile.createDimension('scalar_dim', len(endTime_overpass2))
    rec_dim = outfile.createDimension('rec_dim', None)
    # create data/coordinate variables
    lon_dataset = outfile.createVariable('longitude_dim', np.float32, 
                                         ('longitude_dim',)) 
    lat_dataset = outfile.createVariable('latitude_dim', np.float32, 
                                         ('latitude_dim',)) 
    begTime_dataset = outfile.createVariable('begTime_overpass2', np.float32, 
                                         ('scalar_dim',)) 
    endTime_dataset = outfile.createVariable('endTime_overpass2', np.float32, 
                                         ('scalar_dim',)) 
    o3_dataset = outfile.createVariable('o3', np.float32, ('rec_dim', 
                                        'latitude_dim', 'longitude_dim'))
    co_dataset = outfile.createVariable('co', np.float32, ('rec_dim', 
                                        'latitude_dim', 'longitude_dim'))
    no_dataset = outfile.createVariable('no', np.float32, ('rec_dim', 
                                        'latitude_dim', 'longitude_dim'))
    no2_dataset = outfile.createVariable('no2', np.float32, ('rec_dim', 
                                        'latitude_dim', 'longitude_dim'))
    cloud_dataset = outfile.createVariable('cloudFraction_overpass', np.float32, 
                                           ('rec_dim', 'latitude_dim', 
                                            'longitude_dim'))
    gridbox_dataset = outfile.createVariable('gridBoxHeight_overpass', np.float32, 
                                             ('rec_dim', 'latitude_dim', 
                                              'longitude_dim'))
    # grow data and define attributes
    # latitude
    lat_dataset[:] = np.stack(lat)
    lat_dataset.long_name = 'latitude'
    lat_dataset.units = 'degrees_north'
    # longitude
    lon_dataset[:] = np.stack(lon)
    lon_dataset.long_name = 'longitude'
    lon_dataset.units = 'degree_east'    
    # Beginning Time for overpass
    begTime_dataset[:] = begTime_overpass2
    begTime_dataset.long_name = 'Beginning Time for overpass'
    # End Time for overpass
    endTime_dataset[:] = endTime_overpass2
    endTime_dataset.long_name = 'End Time for overpass'
    # O3
    o3_dataset[:] = np.vstack(o3_all)
    o3_dataset.long_name = 'mass_fraction_of_ozone_in_air'
    o3_dataset.units = '1'
    # CO
    co_dataset[:] = np.vstack(co_all)
    co_dataset.long_name = 'mass_fraction_of_CO_in_air'
    co_dataset.units = '1'
    # NO
    no_dataset[:] = np.vstack(no_all)
    no_dataset.long_name = 'mass_fraction_of_NO_in_air'
    no_dataset.units = '1'
    # NO2
    no2_dataset[:] = np.vstack(no2_all)
    no2_dataset.long_name = 'mass_fraction_of_NO2_in_air'
    no2_dataset.units = '1'
    # Total cloud fraction
    cloud_dataset[:] = np.vstack(cloudFraction_overpass)
    cloud_dataset.long_name = 'Total cloud fraction'
    cloud_dataset.units = 'unitless'
    # Grid Box Height
    gridbox_dataset[:] = np.vstack(gridBoxHeight_overpass)
    gridbox_dataset.long_name = 'Grid Box Height'
    gridbox_dataset.units = 'm'
    # global attributes
    outfile.title = 'Gmimod overpass variable file for summer (JJA) %s.' %year
    outfile.author = 'Gaige Hunter Kerr'
    outfile.comment = 'Output latitudinally bound by %d-%dN ' %(lat[0], lat[-1]) + \
        'and longitudinally bound by %d-%dE is included.' %(lon[0], lon[-1])
    outfile.close()
    return 
# # # # # # # # # # # # #
#import numpy as np
#path = '/Users/ghkerr/phd/GMI/docs/model_info/'
#gmi_cooper_siteinfo = np.genfromtxt(path + 'GMIstats_forcomp_alt_sep.txt', 
#                                    dtype = 'str', skip_header = 2)
## columns for site information are 
## [:, 0] - GMI code (for GMI file naming convention)
## [:, 1] - Cooper (used in Cooper et al. [2012])
## [:, 2] - Timezone (string, e.g. EST, MST) 
## [:, 3] - Sitenum 
## [:, 4] - GMIaltbox (where to sample idaily, etc GMI grid wrt pressure)
#sites = list(gmi_cooper_siteinfo[:, 0])
#case  = 'HindcastMR2'
#year = 2009
#separateFiles_singleyear_combine(case, year, sites)
case = 'HindcastMR2'
year = 2008
ulat = 50.
llat = 25.
llon = 230.
rlon = 300.
monthly_overpass2_combine(year, case, ulat, llat, llon, rlon)



    
