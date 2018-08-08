#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  3 11:53:48 2018

@author: ghkerr
"""
#    GMT_offset : int
#        The difference in hours and minutes from Coordinated Universal Time 
#        (UTC); i.e., for the Northeastern United States this will be -4 during
#        JJA

#



import numpy as np
import sys
sys.path.append('/Users/ghkerr/phd/GMI/')
import commensurability, gmi_visualizations
years = [2008, 2009, 2010]            
sampling_months = [6, 7, 8]
sampling_hours = [15, 16, 17, 18, 19, 20]
# for the NEUS s
castnet_sites_neus = ['ASH', 'WFM', 'WST', 'APT', 'SAR', 'CTH', 'WSP',
                    'ARE', 'BEL', 'PSU', 'LRL', 'PAR', 'PED', 'SHN', 
                    'CDR', 'VPI', 'MKG', 'KEF']







#sampling_hours_neus = [15, 16, 17, 18, 19, 20]
#sampling_hours_south = [15, 16, 17, 18, 19, 20]
#sampling_hours_midwest = [16, 17, 18, 19, 20, 21]
#sampling_hours_west = [17, 18, 19, 20, 21, 22]


#castnet_sites_south = ['ESP', 'SPD', 'GRS', 'BLR', 'COW', 'PNF', 'GRS', 'SND',
#                       'GAS', 'CVL', 'CAD', 'SUM']
#castnet_sites_west = ['LAV', 'PIN', 'SEK', 'JOT', 'GLR', 'CHA', 'GRC', 'YEL', 
#                      'PND', 'CNT', 'ROM', 'GTH']
#castnet_sites_midwest = ['PRK', 'BVL', 'VIN', 'ALH', 'OXF', 'SAL', 'LYK', 
#                         'DCP', 'UVL', 'ANA']
#
regions = ['northeast']
castnet_sites = [castnet_sites_neus]



for castnet_sites_fr, region in zip(castnet_sites, regions):
#    # # # # open daily mean (11-16 LST average) CASTNet O3 and GMI trace gases
    meandaily2 = commensurability.commensurate_castnet_gmi_ra(castnet_sites_fr, 
        years, sampling_months, sampling_hours)
#    # # # # open MERRA-2 2-meter temperatures which are commensurate with 
#    # CASTNet sites, n.b. MR2 CASTNet data (with no regional average) needs 
#    # to be loaded
#    mr2_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
#    commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2', 
#        years, sampling_months, sampling_hours)
#    (comm_t2m, comm_t10m, comm_u2m, comm_u10m, comm_v2m, comm_v10m, 
#     merra_lats_fr, merra_lons_fr) = commensurability.commensurate_MERRA2(
#        mr2_castnet, castnet_sites_fr, years, sampling_months, sampling_hours)        
#    ## Old: load MERRA-1 2-meter temperature
#    #comm_t2m, merra_lats_fr, merra_lons_fr = \
#    #commensurability.commensurate_t2m(mr2_castnet, castnet_sites_fr, years, 
#    #    sampling_months)
#
#    # open MERRA-2 meteorology at surface and 500 hPa
#    H0, T0, U0, V0, times, lat_merra, lon_merra, pressure_merra = \
#    commensurability.open_inst6_3d_ana_Np(years, 0)
#    H500, T500, U500, V500, times, lat_merra, lon_merra, pressure_merra = \
#    commensurability.open_inst6_3d_ana_Np(years, 16)

    
    
#    del mr2_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr
#    # # # # open AQS NO2 and CO
#    aqs_co, aqs_no2, aqs_o3, aqs_co_coords, aqs_no2_coords, aqs_o3_coords = \
#    commensurability.commensurate_aqstracegas(castnet_sites_fr, years, sampling_months, 
#                                           sampling_hours)
#    # average AQS NO2, CO, and O3 over CASTNet sites
#    aqs_co = np.nanmean(aqs_co, axis = 1)
#    aqs_no2 = np.nanmean(aqs_no2, axis = 1)
#    aqs_o3 = np.nanmean(aqs_o3, axis = 1)
#    # # # # open time- and regionally-averaged diurnal cycles of CASTNet O3 and
#    # GMI trace gases from 
#    diurnal = commensurability.commensurate_castnet_gmi_diurnal_ra(castnet_sites_fr, 
#                                                                   years, 
#                                                                   sampling_months)
#    # # # # open open time- and regionally-averaged diurnal cycles of AQS NO2, 
#    # CO, and O3 (n.b. function 'commensurate_aqstracegas_diurnal' required 
#    # unaveraged values)
#    castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
#    commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2', 
#                                              years, sampling_months, 
#                                              sampling_hours)
#    del mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr
#    aqs_co_d, aqs_no2_d, aqs_o3_d, aqs_co_coords_d, aqs_no2_coords_d, \
#    aqs_o3_coords_d = commensurability.commensurate_aqstracegas_diurnal(
#            castnet, castnet_sites_fr, years, sampling_months)
#    # average AQS CO over years, days and CASTNet sites
#    aqs_co_d = np.nanmean(aqs_co_d, axis = tuple(range(0, 3)))
#    aqs_no2_d = np.nanmean(aqs_no2_d, axis = tuple(range(0, 3)))
#    aqs_o3_d = np.nanmean(aqs_o3_d, axis = tuple(range(0, 3)))
#    
#    
#    
#    
#    
#    
    # # # # visualization functions for daily-averaged values
    # data location, n.b. only works for NEUS 
#    coords = [-84., 36., -65., 49.]    
#    gmi_visualizations.map_gmiaqscastnet(meandaily['CASTNet'], castnet_sites_fr, 
#        aqs_co_coords, aqs_no2_coords, coords[0], coords[1], coords[2], coords[3], 
#        years, sampling_months, sampling_hours, region)
#    # O3 boxplot
#    gmi_visualizations.boxplot_tracegas(meandaily['CASTNet'], meandaily['MR2 O3'], 
#                           meandaily['MR2-CCMI O3'], meandaily['FFIgac2 O3'],
#                           meandaily['FFIgac2-HighRes O3'], years, region, 
#                           'O$_{3}$ [ppbv]', 'CASTNet', 'castnetgmio3')
#    # NO2 boxplot
#    gmi_visualizations.boxplot_tracegas(aqs_no2, meandaily['MR2 NO2'], meandaily['MR2-CCMI NO2'], 
#                           meandaily['FFIgac2 NO2'], 
#                           meandaily['FFIgac2-HighRes NO2'], years, region,
#                           'NO$_{2}$ [ppbv]', 'AQS', 'aqsgmino2')
#    # CO boxplot
#    gmi_visualizations.boxplot_tracegas(aqs_co, meandaily['MR2 CO'], meandaily['MR2-CCMI CO'], 
#                     meandaily['FFIgac2 CO'], meandaily['FFIgac2-HighRes CO'], 
#                     years, region, 'CO [ppmv]', 'AQS', 'aqsgmico')
    # O3 CDF
#    gmi_visualizations.cdf_castneto3gmio3(meandaily['CASTNet'], meandaily['MR2 O3'], 
#                       meandaily['MR2-CCMI O3'], meandaily['FFIgac2 O3'],
#                       meandaily['FFIgac2-HighRes O3'], years, region)
#    # O3 PDF
#    gmi_visualizations.pdf_castneto3gmio3(meandaily['CASTNet'], meandaily['MR2 O3'], 
#                       meandaily['MR2-CCMI O3'], meandaily['FFIgac2 O3'],
#                       meandaily['FFIgac2-HighRes O3'], years, region)
#    # NO2 PDF
#    gmi_visualizations.pdf_aqsno2gmino2(aqs_no2, meandaily['MR2 NO2'], meandaily['MR2-CCMI NO2'], 
#                     meandaily['FFIgac2 NO2'], meandaily['FFIgac2-HighRes NO2'], 
#                     years, region)
#    # CO PDF
#    gmi_visualizations.pdf_aqscogmico(aqs_co, meandaily['MR2 CO'], meandaily['MR2-CCMI CO'], 
#                   meandaily['FFIgac2 CO'], meandaily['FFIgac2-HighRes CO'], 
#                   years, region)
#    # O3 scatterplot 
#    gmi_visualizations.scatter_castneto3_gmio3(meandaily['CASTNet'], meandaily['MR2 O3'], 
#                            meandaily['MR2-CCMI O3'], meandaily['FFIgac2 O3'],
#                            meandaily['FFIgac2-HighRes O3'], years, region)
#    # O3 timeseries 
#    gmi_visualizations.timeseries_castneto3gmio3(meandaily['CASTNet'], meandaily['MR2 O3'], 
#                              meandaily['MR2-CCMI O3'], meandaily['FFIgac2 O3'],
#                              meandaily['FFIgac2-HighRes O3'], years, region)
#    ## O3/temperature relationship  
#    #scatter_t2m_castneto3(comm_t2m, meandaily['CASTNet'], years, region)
#    # diurnal O3
#    gmi_visualizations.diurnal_castneto3gmio3(diurnal['CASTNet'], diurnal['MR2 O3'], 
#                           diurnal['MR2-CCMI O3'], diurnal['FFIgac2 O3'],
#                           diurnal['FFIgac2-HighRes O3'], years, region)
#    # diurnal CO
#    gmi_visualizations.diurnal_aqscogmico(aqs_co_d, diurnal['MR2 CO'], diurnal['MR2-CCMI CO'], 
#                       diurnal['FFIgac2 CO'], diurnal['FFIgac2-HighRes CO'], 
#                       years, region)
#    # diurnal NO2 
#    gmi_visualizations.diurnal_aqsno2gmino2(aqs_no2_d, diurnal['MR2 NO2'], diurnal['MR2-CCMI NO2'], 
#                         diurnal['FFIgac2 NO2'], diurnal['FFIgac2-HighRes NO2'], 
#                         years, region)
#    # PDFs of O3 for all simulations
#    gmi_visualizations.pdf_allgmio3(meandaily['CASTNet'], meandaily['MR2 O3'], 
#        meandaily['MR2-CCMI O3'], meandaily['EGU_T O3'], meandaily['DAT O3'], 
#        meandaily['MERRA O3'], region)
#    # plot scatterplot of O3 vs T2m for all simulations
#    gmi_visualizations.scatter_allgmio3(meandaily['CASTNet'], 
#        meandaily['MR2 O3'], meandaily['MR2-CCMI O3'], meandaily['EGU_T O3'], 
#        meandaily['DAT O3'], meandaily['MERRA O3'], t2m_ra, region)
#    # plot difference in O3 on hot/cold days for EGU and MR2 simulations
#    gmi_visualizations.scatter_t2m_dato3eguo3_highlowdays(meandaily['MR2 O3'], 
#        meandaily['DAT O3'], meandaily['EGU_T O3'], 
#        np.nanmean(comm_t2m, axis = 1), region)
#    # plot scatter matrix of meteorology and MR2 O3
#    gmi_visualizations.scatter_matrix_mr2o3meteorology(times, 
#        meandaily['MR2 O3'], H500, U500, V500, comm_t2m, comm_u10m, comm_v10m, 
#        sampling_hours, region)
#    # plot timeseries of O3 and threshold for events from MR2 and DAT 
#    # simulations
#    gmi_visualizations.timeseries_mr2o3dato3(meandaily['MR2 O3'], 
#        meandaily['DAT O3'], region)
#    # plot historgram of O3 and threshold for events from MR2 and DAT 
#    # simulations
#    gmi_visualizations.histogram_mr2o3dato3(meandaily['MR2 O3'], 
#        meandaily['DAT O3'], region)
#    # plot O3-temperature relationship from MR2 and DAT
#    gmi_visualizations.scatter_t2m_mr2o3dato3_deltat2m_delta_o3(
#        meandaily['MR2 O3'], meandaily['DAT O3'], comm_t2m, region)
    
    
    
    
#    
#    
#    
#    
#    
#    
#        
#    
#merra_time = np.copy(times)
#
#
#t2m = comm_t2m
#mr2 = meandaily['MR2 O3']
#egu = meandaily['EGU_T O3']
#dat = meandaily['DAT O3']
#
#
#import pandas as pd
## find indices of MERRA-2 meteorological fields in 'sampling_hours;' n.b.
## this will not work if 'sampling_hours' is not 15-20 local time
#times = pd.to_datetime(times)
#times_intersect = np.in1d(times.hour, sampling_hours)
#H500_intersect = H500[times_intersect]
#U500_intersect = U500[times_intersect]
#V500_intersect = V500[times_intersect]     
#    
#total_wind_intersect = np.hypot(U500_intersect, V500_intersect)    
#
## find regional averages of fields 
#H500_ra = np.nanmean(H500_intersect, axis = tuple((1, 2)))
#U500_ra = np.nanmean(U500_intersect, axis = tuple((1, 2)))
#V500_ra = np.nanmean(V500_intersect, axis = tuple((1, 2)))
#total_wind_ra = np.nanmean(total_wind_intersect, axis = tuple((1, 2)))
#t2m_ra = np.hstack(np.nanmean(t2m, axis = 1))
#mr2_ra = np.hstack(mr2)
#egu_ra = np.hstack(egu)
#dat_ra = np.hstack(dat)
## find high days 
#H500_high = np.where(H500_ra > np.percentile(H500_ra, 80))[0]
#t2m_high = np.where(t2m_ra > np.percentile(t2m_ra, 80))[0]
#U500_high = np.where(U500_ra > np.percentile(U500_ra, 80))[0]
#V500_high = np.where(V500_ra > np.percentile(V500_ra, 80))[0]
#total_wind_high = np.where(total_wind_ra > np.percentile(total_wind_ra, 80))[0]
#
#
#
#    
#wind_abs = np.sqrt(U500_ra**2 + V500_ra**2)
#wind_dir_trig_to = np.arctan2(U500_ra/wind_abs, V500_ra/wind_abs) 
#wind_dir_trig_to_degrees = wind_dir_trig_to * 180/np.pi ## -111.6 degrees    
#    
#    
#
#
#wind090 = np.where((wind_dir_trig_to_degrees > 0) & (wind_dir_trig_to_degrees < 90))[0]
#wind90180 = np.where((wind_dir_trig_to_degrees > 90) & (wind_dir_trig_to_degrees < 180))[0]
#windn0n90 = np.where((wind_dir_trig_to_degrees > -90) & (wind_dir_trig_to_degrees < 0))[0]
#windn90n180= np.where((wind_dir_trig_to_degrees > -180) & (wind_dir_trig_to_degrees < -90))[0]
#
#
#
#






#""" times, sampling_hours
#
#"""
#import pandas as pd
## find indices of MERRA-2 meteorological fields in 'sampling_hours;' n.b.
## this will not work if 'sampling_hours' is not 15-20 local time
#times = pd.to_datetime(times)
#times_intersect = np.in1d(times.hour, sampling_hours)
#H500_intersect = H500[times_intersect]
#U500_intersect = U500[times_intersect]
#V500_intersect = V500[times_intersect]
#
#
## stack yearly time series into one continuous time series
#mr2 = np.hstack(mr2)
#dat = np.hstack(dat)
#mr2_80 = np.percentile(mr2, 80)
## find percentage of extreme events 
#where_high_mr2 = np.where(mr2 > mr2_80)[0]
#where_high_dat = np.where(dat > mr2_80)[0]
## find and plot high extremes in DAT that are not in MR2
#high_only_dat = np.where(np.in1d(where_high_dat, where_high_mr2) == 
#                         False)[0]
#high_only_dat = where_high_dat[high_only_dat]
#
#
## 
#H500_only_dat = H500_intersect[high_only_dat]
#H500_only_dat = np.nanmean(H500_only_dat, axis = 0)
#
#U500_only_dat = U500_intersect[high_only_dat]
#U500_only_dat = np.nanmean(U500_only_dat, axis = 0)
#
#V500_only_dat = V500_intersect[high_only_dat]
#V500_only_dat = np.nanmean(V500_only_dat, axis = 0)
#
#
#
#H500_clim = np.nanmean(H500, axis = 0)
#U500_clim = np.nanmean(U500, axis = 0)
#V500_clim = np.nanmean(V500, axis = 0)
#
#
#plt.contourf(lon_merra, lat_merra, H500_clim - H500_only_dat)
#plt.colorbar()
#plt.show()
#
#
#
#plt.contourf(lon_merra, lat_merra, U500_clim - U500_only_dat)
#plt.colorbar()
#plt.show()
#
#plt.contourf(lon_merra, lat_merra, V500_clim - V500_only_dat)
#plt.colorbar()
#plt.show()
    
    
    
