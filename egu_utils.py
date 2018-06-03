#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAME
    egu_utils.py
PURPOSE
    Script compares control CTM simulation (i.e. HindcastMR2) with sensitivity
    simulation using perturbed NO emissions which scale with 2 meter 
    temperature
PROGRAMMER
    Gaige Hunter Kerr (gaige.kerr@jhu.edu)
REVISION HISTORY
    17042018 -- initial version created
    17042018 -- added functions 'scatter_controlno2eguno2_bysite' and 
                'scatter_controlnoeguno_bysite'
    23042018 -- function 'scatter_inventory_ctm' added 
    26042018 -- function 'scatter_tracegases_runs' and added 
    27042018 -- functions 'scatter_deltanoxdeltao3' and 
                'scatter_deltanoxdeltao3_percent' added
    02052018 -- function 'scatter_deltanoxdeltao3_byyear' added
    21052018 -- function 'timeseries_noemiss_t2m' added 
    29052018 -- function 'scatter_dEmissdNOxdO3' modified to output plots of 
                dNOx/dEmiss and dO3/dNOx
"""
# # # # # # # # # # # # #
def scatter_controlnoeguno_bysite(mr2, egu, castnet_sites_fr, years, region):
    """function plots the percentage change in NO between perturbed NOx 
    emissions and control GMI CTMs runs versus NO concentrations from the 
    control run for every site in variable 'castnet_sites_fr.' Also plotted are 
    the lines of best fit from linear regression and the fitted slope/
    intercept. 
  
    Parameters
    ----------    
    mr2 : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from 'control' model run (i.e. 
        HindcastMR2), units of parts per billion volume, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']    
    egu : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from 'control' model run (i.e. 
        HindcastMR2), units of parts per billion volume, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']    
    castnet_sites_fr : list
        CASTNET site names in focus region
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
    # initialize figure, axes
    fig, axs = plt.subplots(3, 5, figsize = (15, 6))
    axs = axs.ravel()
    # index for axis
    j = 0
    # loop through CASTNet sites in region 
    for site, i in zip(castnet_sites_fr, np.arange(0, len(castnet_sites_fr), 1)):
        # fetch trace gas at site
        mr2_atsite = mr2[:, i]
        egu_atsite = egu[:, i]
        # normalize control by itself (always 1) and perturbed emissions run 
        # by the control run (this will allow for greater values on days where 
        # the perturbed run has higher trace gas concentrations and vice versa)
        egu_atsite = egu_atsite/mr2_atsite * 100.
        # only plot if CASTNet data at site exists
        if np.unique(np.isnan(mr2_atsite))[0] == False:
            axs[j].scatter(mr2_atsite, egu_atsite, s = 5, c = 'k')
            fit = np.polyfit(np.hstack(mr2_atsite), np.hstack(egu_atsite), 
                             deg = 1)
            x = np.linspace(mr2_atsite.min() - 0.01, mr2_atsite.max() + 0.01, 100)
            axs[j].plot(x, fit[0] * x + fit[1], color = 'darkred')
            # print site name abbreviation and slope/intercept of linear fit
            axs[j].set_title(site + ', m = %.1f, b = %.1f' %(fit[0], fit[1]))        
            # axis limits and ticks
            axs[j].set_ylim([75, 125])
            axs[j].set_yticks([75, 100, 125])        
            axs[j].set_ylabel('Change [%]')
            # set maximum number of tick labels on subplots' x axes 
            axs[j].locator_params(axis = 'x', nbins = 5)
            # remove unwanted y-axis tick labels and axis labels
            if (j == 1) or (j == 2) or (j == 3) or (j == 4) or (j == 6) or \
            (j == 7) or (j == 8) or (j == 9) or (j == 11) or (j == 12) or \
            (j == 13) or (j == 14) or (j == 15):
                axs[j].set_yticklabels([''])        
                axs[j].set_ylabel('')
            if (j == 10) or (j == 11) or (j == 12) or (j == 13) or (j == 14):
                axs[j].set_xlabel('NO [ppbv]')
            j = j + 1
    plt.tight_layout()        
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_controlnoeguno_bysite_%d-%d_%s.eps' %(years[0], 
                years[-1], region), dpi = 300)    
    return 
# # # # # # # # # # # # #
def scatter_controlno2eguno2_bysite(mr2, egu, castnet_sites_fr, years, region):
    """function plots the percentage change in NO2 between perturbed NOx 
    emissions and control GMI CTMs runs versus NO2 concentrations from the 
    control run. Also plotted are the lines of best fit from linear regression
    and the fitted slope/intercept. 
  
    Parameters
    ----------    
    mr2 : numpy.ndarray
        GMI CTM NO2 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from 'control' model run (i.e. 
        HindcastMR2), units of parts per billion volume, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']    
    egu : numpy.ndarray
        GMI CTM NO2 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from 'control' model run (i.e. 
        HindcastMR2), units of parts per billion volume, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']    
    castnet_sites_fr : list
        CASTNET site names in focus region
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
    # initialize figure, axes
    fig, axs = plt.subplots(3, 5, figsize = (15, 6))
    axs = axs.ravel()
    # index for axis
    j = 0
    # loop through CASTNet sites in region 
    for site, i in zip(castnet_sites_fr, np.arange(0, len(castnet_sites_fr), 1)):
        # fetch trace gas at site
        mr2_atsite = mr2[:, i]
        egu_atsite = egu[:, i]
        # normalize control by itself (always 1) and perturbed emissions run 
        # by the control run (this will allow for greater values on days where 
        # the perturbed run has higher trace gas concentrations and vice versa)
        egu_atsite = egu_atsite/mr2_atsite * 100.
        # only plot if CASTNet data at site exists
        if np.unique(np.isnan(mr2_atsite))[0] == False:
            axs[j].scatter(mr2_atsite, egu_atsite, s = 5, c = 'k')
            fit = np.polyfit(np.hstack(mr2_atsite), np.hstack(egu_atsite), 
                             deg = 1)
            x = np.linspace(mr2_atsite.min() - 0.01, mr2_atsite.max() + 0.01, 100)
            axs[j].plot(x, fit[0] * x + fit[1], color = 'darkred')
            # print site name abbreviation and slope/intercept of linear fit
            axs[j].set_title(site + ', m = %.1f, b = %.1f' %(fit[0], fit[1]))        
            # axis limits and ticks
            axs[j].set_ylim([75, 125])
            axs[j].set_yticks([75, 100, 125])        
            axs[j].set_ylabel('Change [%]')
            # set maximum number of tick labels on subplots' x axes 
            axs[j].locator_params(axis = 'x', nbins = 5)
            # remove unwanted y-axis tick labels and axis labels
            if (j == 1) or (j == 2) or (j == 3) or (j == 4) or (j == 6) or \
            (j == 7) or (j == 8) or (j == 9) or (j == 11) or (j == 12) or \
            (j == 13) or (j == 14) or (j == 15):
                axs[j].set_yticklabels([''])        
                axs[j].set_ylabel('')
            if (j == 10) or (j == 11) or (j == 12) or (j == 13) or (j == 14):
                axs[j].set_xlabel('NO$_{2}$ [ppbv]')
            j = j + 1
    plt.tight_layout()        
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_controlno2eguno2_bysite_%d-%d_%s.eps' %(years[0], 
                years[-1], region), dpi = 300)    
    return 
# # # # # # # # # # # # #
def scatter_controlo3eguno3_bysite(mr2, egu, castnet_sites_fr, years, region):
    """function plots the percentage change in O3 between perturbed NOx 
    emissions and control GMI CTMs runs versus O3 concentrations from the 
    control run. Also plotted are the lines of best fit from linear regression
    and the fitted slope/intercept. 
  
    Parameters
    ----------    
    mr2 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from 'control' model run (i.e. 
        HindcastMR2), units of parts per billion volume, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']    
    egu : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from 'control' model run (i.e. 
        HindcastMR2), units of parts per billion volume, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']    
    castnet_sites_fr : list
        CASTNET site names in focus region
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
    # initialize figure, axes
    fig, axs = plt.subplots(3, 5, figsize = (15, 6))
    axs = axs.ravel()
    # index for axis
    j = 0
    # loop through CASTNet sites in region 
    for site, i in zip(castnet_sites_fr, np.arange(0, len(castnet_sites_fr), 1)):
        # fetch trace gas at site
        mr2_atsite = mr2[:, i]
        egu_atsite = egu[:, i]
        # normalize control by itself (always 1) and perturbed emissions run 
        # by the control run (this will allow for greater values on days where 
        # the perturbed run has higher trace gas concentrations and vice versa)
        egu_atsite = egu_atsite/mr2_atsite * 100.
        # only plot if CASTNet data at site exists
        if np.unique(np.isnan(mr2_atsite))[0] == False:
            axs[j].scatter(mr2_atsite, egu_atsite, s = 5, c = 'k')
            fit = np.polyfit(np.hstack(mr2_atsite), np.hstack(egu_atsite), 
                             deg = 1)
            x = np.linspace(mr2_atsite.min() - 0.01, mr2_atsite.max() + 0.01, 100)
            axs[j].plot(x, fit[0] * x + fit[1], color = 'darkred')
            # print site name abbreviation and slope/intercept of linear fit
            axs[j].set_title(site + ', m = %.1f, b = %.1f' %(fit[0], fit[1]))        
            # axis limits and ticks
            axs[j].set_ylim([95, 105])
            axs[j].set_yticks([95, 100, 105])        
            axs[j].set_ylabel('Change [%]')
            # set maximum number of tick labels on subplots' x axes 
            axs[j].locator_params(axis = 'x', nbins = 5)
            # remove unwanted y-axis tick labels and axis labels
            if (j == 1) or (j == 2) or (j == 3) or (j == 4) or (j == 6) or \
            (j == 7) or (j == 8) or (j == 9) or (j == 11) or (j == 12) or \
            (j == 13) or (j == 14) or (j == 15):
                axs[j].set_yticklabels([''])        
                axs[j].set_ylabel('')
            if (j == 10) or (j == 11) or (j == 12) or (j == 13) or (j == 14):
                axs[j].set_xlabel('O$_{3}$ [ppbv]')
            j = j + 1
    plt.tight_layout()        
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_controlno3eguno3_bysite_%d-%d_%s.eps' %(years[0], 
                years[-1], region), dpi = 300)    
    return          
# # # # # # # # # # # # #    
def diurnal_controlegunono2nox(mr2_no_d, mr2_no2_d, egu_no_d, egu_no2_d, 
                               castnet_sites_fr, years, region):
    """function finds the difference between the diurnal curves of NO, NO2, and 
    NOx from the control (i.e. HindcastMR2) and perturbed NOx emissions 
    simulations. Plotted values are Control - Variable Emissions, so positive
    values imply that greater values for a given trace gas are higher in 
    control run. 
    
    Parameters
    ----------   
    mr2_no_d : numpy.ndarray
        Regionally- and time-averaged NO diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2, units of ppmv, [24,]        
    mr2_no2_d : numpy.ndarray
        Regionally- and time-averaged NO2 diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2, units of ppmv, [24,]      
    egu_no_d : numpy.ndarray
        Regionally- and time-averaged NO diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2 with added time-varying NOx 
        emissions, units of ppmv, [24,]       
    egu_no2_d : numpy.ndarray
        Regionally- and time-averaged NO2 diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2 with added time-varying NOx 
        emissions, units of ppmv, [24,]      
    castnet_sites_fr : list
        CASTNET site names in focus region
    years : list
        Years of interest    
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function

    Returns
    ----------          
    None    
    """
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
    path = pollutants_constants.FONTPATH_BOLD
    prop = font_manager.FontProperties(fname = path)
    mpl.rcParams['mathtext.bf'] = prop.get_name()
    # initialize figure, axis
    fig = plt.figure(figsize = (9, 4))
    ax = plt.subplot2grid((1, 2), (0, 0), colspan = 2)
    ax.set_title('Control - Variable Emissions')
    # for NO
    ax.plot(mr2_no_d - egu_no_d, lw = 1.5, color = '#66c2a5', 
            label = 'NO')
    # for NO2
    ax.plot(mr2_no2_d - egu_no2_d, lw = 1.5, color = '#fc8d62', 
            label = 'NO$_{2}$')
    # for NOx = NO + NO2
    ax.plot((mr2_no2_d + mr2_no_d) - (egu_no_d + egu_no2_d), 
            lw = 1.5, color = '#8da0cb', label = 'NO$_{x}$')
    ax.legend(loc = 'upper left', frameon = False, fontsize = 10)
    ax.set_xlim([0, 23])
    ax.set_xlabel('Hour [GMT]')
    ax.set_ylabel('Difference [ppbv]')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'diurnal_controlegunono2nox_%d-%d_%s.eps' %(years[0], years[-1],
                                                           region), dpi = 300)
    return
# # # # # # # # # # # # #    
def scatter_inventory_ctm(no_emiss_perturbed_ra, no_emiss_unperturbed_ra, 
                          egu_no_ra, mr2_no_ra, egu_no2_ra, mr2_no2_ra, 
                          egu_o3_ra, mr2_o3_ra, years, region): 
    """function finds the percentage difference between regionally-averaged
    NO from perturbed NOx emissions and control emissions inventory and 
    the percentage difference in modeled NO, NO2, and O3 between both 
    simulations. These relationships are plotted along with regression lines 
    and their slopes.
   
    Parameters
    ----------    
    no_emiss_perturbed_ra : numpy.ndarray
        Regionally-averaged NO emissions at grid cells co-located (or nearly
        co-located) with CASTNet sites in region from perturbed NO emissions 
        simulation which vary on a daily basis, units of mol/mol, [years in 
        measuring period, days in months in 'sampling_months']    
    no_emiss_unperturbed_ra : numpy.ndarray
        Regionally-averaged NO emissions at grid cells co-located (or nearly
        co-located) with CASTNet sites in region from perturbed NO emissions 
        simulation which vary on a monthly basis, units of mol/mol, [years in 
        measuring period, days in months in 'sampling_months']        
   egu_no_ra : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from perturbed NO emissions simulation, 
        units of parts per billion volume, [years in measuring period, stations 
        in 'castnet_sites_fr', days in months in 'sampling_months']           
    mr2_no_ra : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from 'control' model run (i.e. 
        HindcastMR2), units of parts per billion volume, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']        
    egu_no2_ra : numpy.ndarray
        Same as 'egu_no_ra' but for NO2
    mr2_no2_ra : numpy.ndarray
        Same as 'mr2_no_ra' but for NO2
    egu_o3_ra : numpy.ndarray
        Same as 'egu_no_ra' but for O3    
    mr2_o3_ra : numpy.ndarray
        Same as 'mr2_no_ra' but for O3
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
    # initialize figure, axis
    fig = plt.figure(figsize = (7, 7))
    ax = plt.subplot2grid((1, 1), (0, 0))
    # plot percentage change in trace gas between control and perturbed NOx
    # emissions simulation vs. percentage change of NO emissions between control 
    # and perturbed NOx emissions emissions inventory
    ax.plot((np.hstack(no_emiss_perturbed_ra)/np.hstack(no_emiss_unperturbed_ra)), 
            (np.hstack(egu_no_ra)/np.hstack(mr2_no_ra)), 'o', markersize = 4, 
            color = '#7570b3', label = 'NO')
    m_no, b_no = np.polyfit((np.hstack(no_emiss_perturbed_ra)/np.hstack(no_emiss_unperturbed_ra)), 
                            (np.hstack(egu_no_ra)/np.hstack(mr2_no_ra)), 1)
    ax.plot((np.hstack(no_emiss_perturbed_ra)/np.hstack(no_emiss_unperturbed_ra)), 
            (np.hstack(egu_no2_ra)/np.hstack(mr2_no2_ra)), 'o', markersize = 4, 
            color = '#d95f02', label = 'NO$_{2}$')
    m_no2, b_no2 = np.polyfit((np.hstack(no_emiss_perturbed_ra)/np.hstack(no_emiss_unperturbed_ra)), 
                              (np.hstack(egu_no2_ra)/np.hstack(mr2_no2_ra)), 1)
    ax.plot((np.hstack(no_emiss_perturbed_ra)/np.hstack(no_emiss_unperturbed_ra)), 
            (np.hstack(egu_o3_ra)/np.hstack(mr2_o3_ra)), 'o', markersize = 4, 
            color = '#1b9e77', label = 'O$_{3}$')
    m_o3, b_o3 = np.polyfit((np.hstack(no_emiss_perturbed_ra)/np.hstack(no_emiss_unperturbed_ra)), 
                            (np.hstack(egu_o3_ra)/np.hstack(mr2_o3_ra)), 1)
    # overlay regression lines
    x = np.arange(0.8, 1.15, 0.001)
    ax.plot(x, m_no*x + b_no, '-', lw = 2., color = '#7570b3')
    ax.plot(x, m_no2*x + b_no2, '-', lw = 2., color = '#d95f02')
    ax.plot(x, m_o3*x + b_o3, '-', lw = 2., color = '#1b9e77')
    # add slope of regression lines 
    ax.text(0.05, 0.9, '%.2f' %m_no, color = '#7570b3',
            transform = ax.transAxes, fontsize = 16)
    ax.text(0.05, 0.84, '%.2f' %m_no2, color = '#d95f02',
            transform = ax.transAxes, fontsize = 16)
    ax.text(0.05, 0.78, '%.2f' %m_o3, color = '#1b9e77',
            transform = ax.transAxes, fontsize = 16)
    ax.set_xlim([0.75, 1.2])
    ax.set_ylim([0.75, 1.2])
    plt.tick_params(axis = 'both', which = 'major', labelsize = 14)
    ax.set_xlabel('NO$_{\mathregular{ff, Daily Varying}}$ : NO$_{\mathregular{ff, HindcastMR2}}$', fontsize = 16)
    ax.set_ylabel('$\chi_{\mathregular{Daily Varying}}$ : $\chi_{\mathregular{HindcastMR2}}$', fontsize = 16)
    leg = ax.legend(loc = 'lower right', fontsize = 16)
    leg.get_frame().set_linewidth(0.0)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_inventory_ctm_%d-%d_%s.eps' %(years[0], years[-1], region), 
                dpi = 300)
# # # # # # # # # # # # #    
def scatter_tracegases_runs(no_emiss_perturbed, no_emiss_unperturbed, egu_no, 
                            mr2_no, egu_no2, mr2_no2, egu_o3, mr2_o3,
                            t2m, years, region): 
    """He et al. [2012] gives the following equations: 
    d[O3]/dT = dEmissions[NOx]/dT x d[NOx]/dEmissions[NOx] x OPE
    This function plots the first three terms (all besides OPE) for the control
    and perturbed NO emissions simulations and finds the slope of the linear
    regression between all the quantities. 
    
    Parameters
    ----------    
    no_emiss_perturbed : numpy.ndarray
        Regionally-averaged NO emissions at grid cells co-located (or nearly
        co-located) with CASTNet sites in region from perturbed NO emissions 
        simulation which vary on a daily basis, units of mol/mol, [years in 
        measuring period, days in months in 'sampling_months']    
    no_emiss_unperturbed : numpy.ndarray
        Regionally-averaged NO emissions at grid cells co-located (or nearly
        co-located) with CASTNet sites in region from perturbed NO emissions 
        simulation which vary on a monthly basis, units of mol/mol, [years in 
        measuring period, days in months in 'sampling_months']        
    egu_no : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from perturbed NO emissions simulation, 
        units of parts per billion volume, [years in measuring period, stations 
        in 'castnet_sites_fr', days in months in 'sampling_months']           
    mr2_no : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from 'control' model run (i.e. 
        HindcastMR2), units of parts per billion volume, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']        
    egu_no2 : numpy.ndarray
        Same as 'egu_no_ra' but for NO2
    mr2_no2 : numpy.ndarray
        Same as 'mr2_no_ra' but for NO2
    egu_o3 : numpy.ndarray
        Same as 'egu_no_ra' but for O3    
    mr2_o3 : numpy.ndarray
        Same as 'mr2_no_ra' but for O3
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
    import matplotlib as mpl
    import matplotlib.font_manager as font_manager    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # set custom font
    path = pollutants_constants.FONTPATH_LIGHT
    prop = font_manager.FontProperties(fname=path)
    mpl.rcParams['font.family'] = prop.get_name()
    path = pollutants_constants.FONTPATH_LIGHT
    mpl.rcParams['mathtext.bf'] = prop.get_name()
    # calculate NOx
    egu_nox = egu_no + egu_no2
    mr2_nox = mr2_no + mr2_no2    
    # initialize figure, axes
    fig = plt.figure(figsize = (10, 8))
    ax1 = plt.subplot2grid((2, 3), (0, 0))
    ax2 = plt.subplot2grid((2, 3), (0, 1))
    ax3 = plt.subplot2grid((2, 3), (0, 2))
    ax4 = plt.subplot2grid((2, 3), (1, 0))
    ax5 = plt.subplot2grid((2, 3), (1, 1))
    ax6 = plt.subplot2grid((2, 3), (1, 2))
    # plot NO emissions vs. 2 meter temperature for control simulation
    ax1.plot(np.hstack(t2m), np.hstack(no_emiss_unperturbed), 'ko', markersize = 4)
    m_emiss_t = np.polyfit(np.hstack(t2m), np.hstack(no_emiss_perturbed), 1)[0]
    ax1.text(0.05, 0.85, 'm = %.4f' %m_emiss_t, color = 'k',
             transform = ax1.transAxes, fontsize = 12)
    ax1.set_xlabel('T$_{2 m}$ [K]')
    ax1.set_ylabel('Control\nNO Emissions [ppbv]')
    # plot NOx vs. NO emissions 
    ax2.plot(np.hstack(no_emiss_unperturbed), np.hstack(mr2_nox), 'ko', markersize = 4)
    m_nox_emiss = np.polyfit(np.hstack(no_emiss_unperturbed), np.hstack(mr2_nox), 1)[0]
    ax2.text(0.05, 0.85, 'm = %.4f' %m_nox_emiss, color = 'k',
             transform = ax2.transAxes, fontsize = 12)
    ax2.set_xlabel('NO Emissions [ppbv]')
    ax2.set_ylabel('NO$_{x}$ [ppbv]')
    # plot O3 vs 2 meter temperature
    ax3.plot(np.hstack(t2m), np.hstack(mr2_o3), 'ko', markersize = 4)
    m_o3_t = np.polyfit(np.hstack(t2m), np.hstack(mr2_o3), 1)[0]
    ax3.text(0.05, 0.85, 'm = %.4f' %m_o3_t, color = 'k',
             transform = ax3.transAxes, fontsize = 12)
    ax3.set_xlabel('T$_{2 m}$ [K]')
    ax3.set_ylabel('O$_{3}$ [ppbv]')
    # plot NO emissions vs. 2 meter temperature for perturbed NO simulation
    ax4.plot(np.hstack(t2m), np.hstack(no_emiss_perturbed), 'ko', markersize = 4)
    m_emiss_t = np.polyfit(np.hstack(t2m), np.hstack(no_emiss_perturbed), 1)[0]
    ax4.text(0.05, 0.85, 'm = %.4f' %m_emiss_t, color = 'k',
             transform = ax4.transAxes, fontsize = 12)
    ax4.set_xlabel('T$_{2 m}$ [K]')
    ax4.set_ylabel('Perturbed\nNO Emissions [ppbv]')
    # plot NOx vs. NO emissions 
    ax5.plot(np.hstack(no_emiss_perturbed), np.hstack(egu_nox), 'ko', markersize = 4)
    m_nox_emiss = np.polyfit(np.hstack(no_emiss_perturbed), np.hstack(egu_nox), 1)[0]
    ax5.text(0.05, 0.85, 'm = %.4f' %m_nox_emiss, color = 'k',
             transform = ax5.transAxes, fontsize = 12)
    ax5.set_xlabel('NO Emissions [ppbv]')
    ax5.set_ylabel('NO$_{x}$ [ppbv]')
    # plot O3 vs 2 meter temperature
    ax6.plot(np.hstack(t2m), np.hstack(egu_o3), 'ko', markersize = 4)
    m_o3_t = np.polyfit(np.hstack(t2m), np.hstack(egu_o3), 1)[0]
    ax6.text(0.05, 0.85, 'm = %.4f' %m_o3_t, color = 'k',
             transform = ax6.transAxes, fontsize = 12)
    ax6.set_xlabel('T$_{2 m}$ [K]')
    ax6.set_ylabel('O$_{3}$ [ppbv]')
    plt.tight_layout()
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_tracegases_runs_%d-%d_%s.eps' %(years[0], years[-1], 
                                                         region), dpi = 300)
    return
# # # # # # # # # # # # #
#def scatter_deltanoxdeltao3(egu_no, mr2_no, ffigac2_no, emfix_no, egu_no2, 
#                            mr2_no2, ffigac2_no2, emfix_no2, egu_o3, mr2_o3,
#                            ffigac2_o3, emfix_o3, years, region):
#    """function plots the difference in O3 versus the difference in NOx for 
#    two pairs of models: (left subplot) Fixed emissions - standard simulations
#    from Strode et al. [2015] and (right subplot) perturbed NO emissions - 
#    control simulation. For both pairs of model, the slope of the linear 
#    regression is plotted. 
#
#    Parameters
#    ----------         
#    egu_no : numpy.ndarray
#        GMI CTM NO concentrations co-located (or nearly colocated) with 
#        corresponding CASTNet stations from perturbed NO emissions simulation, 
#        units of parts per billion volume, [years in measuring period, stations 
#        in 'castnet_sites_fr', days in months in 'sampling_months']           
#    mr2_no : numpy.ndarray
#        GMI CTM NO concentrations co-located (or nearly colocated) with 
#        corresponding CASTNet stations from 'control' model run (i.e. 
#        HindcastMR2), units of parts per billion volume, [years in measuring 
#        period, stations in 'castnet_sites_fr', days in months in 
#        'sampling_months']   
#    ffigac2_no : numpy.ndarray
#        GMI CTM NO concentrations co-located (or nearly colocated) with 
#        corresponding CASTNet stations from Strode et al. [2015] 'Std.' run 
#        (i.e. HindcastFFIgac2), units of parts per billion volume, 
#        [years in measuring period, stations in 'castnet_sites_fr', days in 
#        months in 'sampling_months']   
#    emfix_no : numpy.ndarray
#        GMI CTM NO concentrations co-located (or nearly colocated) with 
#        corresponding CASTNet stations from Strode et al. [2015] 'EmFix' run 
#        (i.e. Hindcast3Igac2), units of parts per billion volume, 
#        [years in measuring period, stations in 'castnet_sites_fr', days in 
#        months in 'sampling_months']           
#    egu_no2 : numpy.ndarray
#        Same as 'egu_no' but for NO2
#    mr2_no2 : numpy.ndarray
#        Same as 'mr2_no' but for NO2       
#    ffigac2_no2 : numpy.ndarray 
#        Same as 'ffigac2_no' but for NO2   
#    emfix_no2 : numpy.ndarray
#        Same as 'emfix_no' but for NO       
#    egu_o3 : numpy.ndarray
#        Same as 'egu_no' but for O3    
#    mr2_o3 : numpy.ndarray
#        Same as 'mr2_no' but for O3
#    ffigac2_o3 : numpy.ndarray
#        Same as 'ffigac2_no' but for O3    
#    emfix_o3 : numpy.ndarray
#        Same as 'emfix_no' but for O3    
#    years : list
#        Years of interest    
#    region : str
#        Region over which regionally-averaged concentrations are supplied to 
#        function
#
#    Returns
#    ----------          
#    None             
#    """
#    import numpy as np
#    import matplotlib.pyplot as plt
#    # calculate NOx (= NO + NO2)
#    egu_nox = egu_no + egu_no2
#    mr2_nox = mr2_no + mr2_no2
#    emfix_nox = emfix_no + emfix_no2
#    ffigac2_nox = ffigac2_no + ffigac2_no2    
#    # initialize figure, axes
#    fig = plt.figure(figsize = (8, 4))
#    ax1 = plt.subplot2grid((1, 2), (0, 0))
#    ax2 = plt.subplot2grid((1, 2), (0, 1))
#    # difference in O3 vs NOx between EmFix and Std simulations (i.e. 
#    # Hindcast3Igac2 - HindcastFFIgac2)
#    ax1.plot(np.hstack(emfix_nox - ffigac2_nox), 
#             np.hstack(emfix_o3 - ffigac2_o3), 'ko', markersize = 4)
#    ax1.set_title('EmFix - Std')
#    m, b = np.polyfit(np.hstack(emfix_nox - ffigac2_nox), 
#                      np.hstack(emfix_o3 - ffigac2_o3), 1)
#    x = np.arange((emfix_nox - ffigac2_nox).min(), (emfix_nox - ffigac2_nox).max(), 0.01)
#    ax1.plot(x, x*m + b, '-k')
#    ax1.text(0.75, 0.9, '%.2f' %m, color = 'k',
#             transform = ax1.transAxes, fontsize = 14)
#    ax1.set_ylabel('ppbv $\Delta$ [O$_{3}$]')
#    ax1.set_xlabel('ppbv $\Delta$ [NO$_{x}]$')
#    # difference in O3 vs NOx between pertubed NO emissions and control simulations 
#    #(i.e. EGU_T - HindcastMR2)
#    ax2.plot(np.hstack(egu_nox - mr2_nox), 
#             np.hstack(egu_o3 - mr2_o3), 'ko', markersize = 4)
#    ax2.set_title('CEMS NO - control')
#    m, b = np.polyfit(np.hstack(egu_nox - mr2_nox), 
#                      np.hstack(egu_o3 - mr2_o3), 1)
#    x = np.arange((egu_nox - mr2_nox).min(), (egu_nox - mr2_nox).max(), 0.01)
#    ax2.plot(x, x*m + b, '-k')
#    ax2.text(0.05, 0.9, '%.2f' %m, color = 'k',
#             transform = ax2.transAxes, fontsize = 14)
#    ax2.set_xlabel('ppbv $\Delta$ [NO$_{x}]$')
#    # fit linear regression for delta NOx values greater than 0.0
#    g0 = np.where((egu_nox - mr2_nox) > 0.)
#    m_g0, b_g0 = np.polyfit(np.hstack((egu_nox - mr2_nox)[g0]), 
#                            np.hstack((egu_o3 - mr2_o3)[g0]), 1)
#    x_g0 = np.arange((egu_nox - mr2_nox)[g0].min(), 
#                  (egu_nox - mr2_nox)[g0].max(), 0.01)
#    ax2.plot(x_g0, x_g0*m_g0 + b_g0, '#1f78b4')
#    ax2.text(0.05, 0.84, '%.2f' %m_g0, color = '#1f78b4',
#            transform = ax2.transAxes, fontsize = 14)
#    # for values less than or equal to 0.0
#    l0 = np.where((egu_nox - mr2_nox) <= 0.)
#    m_l0, b_l0 = np.polyfit(np.hstack((egu_nox - mr2_nox)[l0]), 
#                            np.hstack((egu_o3 - mr2_o3)[l0]), 1)
#    x_l0 = np.arange((egu_nox - mr2_nox)[l0].min(), 
#                  (egu_nox - mr2_nox)[l0].max(), 0.01)
#    ax2.plot(x_l0, x_l0*m_l0 + b_l0, '#b2df8a')
#    ax2.text(0.05, 0.78, '%.2f' %m_l0, 
#             color = '#b2df8a', transform = ax2.transAxes, fontsize = 14)         
#    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
#                'scatter_deltanoxdeltao3_%d-%d_%s.eps' %(years[0], years[-1], 
#                                                         region), dpi = 300)         
#    return 
# # # # # # # # # # # # #
def scatter_deltanoxdeltao3_percent(egu_no, mr2_no, ffigac2_no, emfix_no, 
                                    egu_no2, mr2_no2, ffigac2_no2, emfix_no2, 
                                    egu_o3, mr2_o3, ffigac2_o3, emfix_o3, 
                                    years, region):
    """same as 'scatter_deltanoxdeltao3_percent,' but function finds the 
    percentage difference of NOx and O3 differences between simulations; i.e. 
    %\Delta O3 = (O3_{EmFix} - O3_{Std})/O3_{Std}

    Parameters
    ----------         
    egu_no : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from perturbed NO emissions simulation, 
        units of parts per billion volume, [years in measuring period, stations 
        in 'castnet_sites_fr', days in months in 'sampling_months']           
    mr2_no : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from 'control' model run (i.e. 
        HindcastMR2), units of parts per billion volume, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']   
    ffigac2_no : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from Strode et al. [2015] 'Std.' run 
        (i.e. HindcastFFIgac2), units of parts per billion volume, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months']   
    emfix_no : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from Strode et al. [2015] 'EmFix' run 
        (i.e. Hindcast3Igac2), units of parts per billion volume, 
        [years in measuring period, stations in 'castnet_sites_fr', days in 
        months in 'sampling_months']           
    egu_no2 : numpy.ndarray
        Same as 'egu_no' but for NO2
    mr2_no2 : numpy.ndarray
        Same as 'mr2_no' but for NO2       
    ffigac2_no2 : numpy.ndarray 
        Same as 'ffigac2_no' but for NO2   
    emfix_no2 : numpy.ndarray
        Same as 'emfix_no' but for NO       
    egu_o3 : numpy.ndarray
        Same as 'egu_no' but for O3    
    mr2_o3 : numpy.ndarray
        Same as 'mr2_no' but for O3
    ffigac2_o3 : numpy.ndarray
        Same as 'ffigac2_no' but for O3    
    emfix_o3 : numpy.ndarray
        Same as 'emfix_no' but for O3    
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
    # calculate NOx (= NO + NO2)
    egu_nox = egu_no + egu_no2
    mr2_nox = mr2_no + mr2_no2
    emfix_nox = emfix_no + emfix_no2
    ffigac2_nox = ffigac2_no + ffigac2_no2    
    # initialize figure, axes
    fig = plt.figure(figsize = (8, 4))
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    # difference in O3 vs NOx between EmFix and Std simulations (i.e. 
    # Hindcast3Igac2 - HindcastFFIgac2)
    ax1.plot(np.hstack((emfix_nox - ffigac2_nox)/ffigac2_nox) * 100., 
             np.hstack((emfix_o3 - ffigac2_o3)/ffigac2_o3) * 100., 'ko', markersize = 4)
    ax1.set_title('EmFix - Std')
    m, b = np.polyfit(np.hstack((emfix_nox - ffigac2_nox)/ffigac2_nox) * 100., 
                      np.hstack((emfix_o3 - ffigac2_o3)/ffigac2_o3) * 100., 1)
    x = np.arange(((emfix_nox - ffigac2_nox)/ffigac2_nox).min() * 100., 
                  ((emfix_nox - ffigac2_nox)/ffigac2_nox).max() * 100., 0.01)
    ax1.plot(x, x*m + b, '-k')
    ax1.text(0.75, 0.9, '%.2f' %m, color = 'k',
             transform = ax1.transAxes, fontsize = 14)
    ax1.set_ylabel('%$\Delta$ [O$_{3}$]')
    ax1.set_xlabel('%$\Delta$ [NO$_{x}]$')
    # percent difference in O3 vs NOx between pertubed NO emissions and control 
    # simulations (i.e. EGU_T - HindcastMR2)
    ax2.plot(np.hstack((egu_nox - mr2_nox)/mr2_nox) * 100., 
             np.hstack((egu_o3 - mr2_o3)/mr2_o3) * 100., 'ko', markersize = 4)
    ax2.set_title('CEMS NO - control')
    m, b = np.polyfit(np.hstack((egu_nox - mr2_nox)/mr2_nox) * 100., 
                      np.hstack((egu_o3 - mr2_o3)/mr2_o3) * 100., 1)
    x = np.arange(np.hstack((egu_nox - mr2_nox)/mr2_nox).min() * 100., 
                  np.hstack((egu_nox - mr2_nox)/mr2_nox).max() * 100., 0.01)
    ax2.plot(x, x*m + b, '-k')
    ax2.text(0.05, 0.9, '%.2f' %m, color = 'k',
             transform = ax2.transAxes, fontsize = 14)
    ax2.set_xlabel('%$\Delta$ [NO$_{x}]$')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_deltanoxdeltao3_percent_%d-%d_%s.eps' 
                %(years[0], years[-1], region), dpi = 300)
    return
# # # # # # # # # # # # #
def diurnal_controlegunono2o3(mr2_no_d_ra, egu_no_d_ra, mr2_no2_d_ra, 
                              egu_no2_d_ra, mr2_o3_d_ra, egu_o3_d_ra, t2m, 
                              years, region): 
    """function plots the regionally-averaged diurnal curves of NO, NO2, and 
    O3 on hot (2 meter temperature > 80th percentile) and cool (2 meter 
    temperature < 20th percentile) days. Curves are plotted for both control 
    and perturbed NO emissions simulations.
    
    Parameters
    ----------   
    mr2_no_d_ra : numpy.ndarray
        Regionally- NO diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2, units of ppmv, [years in 
        measuring period, days in months in 'sampling_months', 24,]     
    egu_no_d_ra : numpy.ndarray
        Regionally- NO diurnal cycle concentrations from GMI 
        CTM model configuration EGU_T, units of ppmv, [years in 
        measuring period, days in months in 'sampling_months', 24,]        
    mr2_no2_d : numpy.ndarray
        Regionally- NO2 diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2, units of ppmv, [years in 
        measuring period, days in months in 'sampling_months', 24,]   
    egu_no2_d_ra : numpy.ndarray
        Regionally- NO2 diurnal cycle concentrations from GMI 
        CTM model configuration EGU_T, units of ppmv, [years in 
        measuring period, days in months in 'sampling_months', 24,]              
    mr2_o3_d_ra : numpy.ndarray   
        Regionally- O3 diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2, units of ppmv, [years in 
        measuring period, days in months in 'sampling_months', 24,]       
    egu_o3_d_ra : numpy.ndarray
        Regionally- O3 diurnal cycle concentrations from GMI 
        CTM model configuration EGU_T, units of ppmv, [years in 
        measuring period, days in months in 'sampling_months', 24,] 
    t2m : numpy.ndarray    
        Regionally-averaged MERRA 2 meter temperatures co-located (or nearly 
        colocated) with corresponding CASTNet stations, units of K, [years in 
        measuring period, days in months in 'sampling_months']     
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
    import matplotlib.lines as mlines
    import matplotlib.patches as mpatches
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
    path = pollutants_constants.FONTPATH_BOLD
    prop = font_manager.FontProperties(fname = path)
    mpl.rcParams['mathtext.bf'] = prop.get_name()
    # find hot (top 20%-ile) and cool (bottom 20%-ile) days 
    hotday = np.where(t2m >= np.percentile(t2m, 80))
    coolday = np.where(t2m <= np.percentile(t2m, 20))
    # find NO, NO2, and O3 on hot and cool days
    mr2_no_hot = mr2_no_d_ra[hotday]
    mr2_no_cool = mr2_no_d_ra[coolday]
    mr2_no2_hot = mr2_no2_d_ra[hotday]
    mr2_no2_cool = mr2_no2_d_ra[coolday]
    mr2_o3_hot = mr2_o3_d_ra[hotday]
    mr2_o3_cool = mr2_o3_d_ra[coolday]
    egu_no_hot = egu_no_d_ra[hotday]
    egu_no_cool = egu_no_d_ra[coolday]
    egu_no2_hot = egu_no2_d_ra[hotday]
    egu_no2_cool = egu_no2_d_ra[coolday]
    egu_o3_hot = egu_o3_d_ra[hotday]
    egu_o3_cool = egu_o3_d_ra[coolday]
    # average over all hot/cool days
    mr2_no_hot = np.nanmean(mr2_no_hot, axis = 0)
    mr2_no_cool = np.nanmean(mr2_no_cool, axis = 0)
    mr2_no2_hot = np.nanmean(mr2_no2_hot, axis = 0)
    mr2_no2_cool = np.nanmean(mr2_no2_cool, axis = 0)
    mr2_o3_hot = np.nanmean(mr2_o3_hot, axis = 0)
    mr2_o3_cool = np.nanmean(mr2_o3_cool, axis = 0)
    egu_no_hot = np.nanmean(egu_no_hot, axis = 0)
    egu_no_cool = np.nanmean(egu_no_cool, axis = 0)
    egu_no2_hot = np.nanmean(egu_no2_hot, axis = 0)
    egu_no2_cool = np.nanmean(egu_no2_cool, axis = 0)
    egu_o3_hot = np.nanmean(egu_o3_hot, axis = 0)
    egu_o3_cool = np.nanmean(egu_o3_cool, axis = 0)
    # initialize figure, axes
    fig = plt.figure(figsize = (4, 8))
    ax1 = plt.subplot2grid((3, 1), (0, 0))
    ax2 = plt.subplot2grid((3, 1), (1, 0))
    ax3 = plt.subplot2grid((3, 1), (2, 0))
    # difference in modeled NO on hot/cool days between control and perturbed
    # emissions simulations
    ax1.plot(mr2_no_hot, '-r')
    ax1.plot(mr2_no_cool, '-b')
    ax1.plot(egu_no_hot, '-.r')
    ax1.plot(egu_no_cool, '-.b')
    # difference in NO2
    ax2.plot(mr2_no2_hot, '-r')
    ax2.plot(mr2_no2_cool, '-b')
    ax2.plot(egu_no2_hot, '-.r')
    ax2.plot(egu_no2_cool, '-.b')
    # difference in O3
    ax3.plot(mr2_o3_hot, '-r')
    ax3.plot(mr2_o3_cool, '-b')
    ax3.plot(egu_o3_hot, '-.r', label = '')
    ax3.plot(egu_o3_cool, '-.b', label = '')
    ax1.set_xlim([0, 23])
    ax1.set_xticks([0, 6, 12, 18])
    ax1.set_xticklabels([''])
    ax1.set_ylabel('NO [ppbv]')
    ax2.set_xlim([0, 23])
    ax2.set_xticks([0, 6, 12, 18])
    ax2.set_xticklabels([''])
    ax2.set_ylabel('NO$_{2}$ [ppbv]')
    ax3.set_xlim([0, 23])
    ax3.set_xticks([0, 6, 12, 18])
    ax3.set_xlabel('Hour [GMT]')
    ax3.set_ylabel('O$_{3}$ [ppbv]')
    # manually add items to legend
    mr2_line = mlines.Line2D([], [], color = 'k', linestyle = '-',  
                             label='Control')
    egu_line = mlines.Line2D([], [], color = 'k', linestyle = '-.', 
                             label = 'Perturbed NO')
    red_patch = mpatches.Patch(color = 'r', label = 'T$_{2 m}$ > 80%-ile')
    blue_patch = mpatches.Patch(color = 'b', label = 'T$_{2 m}$ < 20%-ile')
    plt.subplots_adjust(bottom = 0.2, left = 0.25)
    leg = ax3.legend(handles = [mr2_line, egu_line, red_patch, blue_patch], 
                     ncol = 2, bbox_to_anchor = (1.1, -0.35))
    leg.get_frame().set_linewidth(0.0)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'diurnal_controlegunono2o3_%d-%d_%s.eps' 
                %(years[0], years[-1], region), dpi = 300)
    return 
# # # # # # # # # # # # #    
def scatter_deltanoxdeltao3_byyear(ffigac2_no, emfix_no, ffigac2_no2, 
                                   emfix_no2, ffigac2_o3, emfix_o3, years, 
                                   region): 
    """for each year in the measuring period 2000 - 2010 function the slope of
    the linear regression line through data points represent the change in 
    O3 between EmFix and standard simulations (y: dO3) and the change in NOx
    between the EmFix and standard simulations (x: NOx).

    Parameters
    ----------   
    ffigac2_no : numpy.ndarray
        GMI CTM NO concentrations regionally-averaged over co-located (or nearly 
        colocated) with corresponding CASTNet stations from Strode et al. 
        [2015] 'Std.' run (i.e. HindcastFFIgac2), units of parts per billion 
        volume, [years in measuring period, days in months in 
        'sampling_months']   
    emfix_no : numpy.ndarray
        GMI CTM NO concentrations regionally-averaged over co-located (or nearly 
        colocated) with  corresponding CASTNet stations from Strode et al. 
        [2015] 'EmFix' run (i.e. Hindcast3Igac2), units of parts per billion 
        volume, [years in measuring period, days in months in 
        'sampling_months']        
    ffigac2_no2 : numpy.ndarray 
        Same as 'ffigac2_no' but for NO2   
    emfix_no2 : numpy.ndarray
        Same as 'emfix_no' but for NO       
    ffigac2_o3 : numpy.ndarray
        Same as 'ffigac2_no' but for O3    
    emfix_o3 : numpy.ndarray
        Same as 'emfix_no' but for O3  
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
    # initialize figure, axes
    fig, axs = plt.subplots(3, 4)#, figsize = (15, 6))
    axs = axs.ravel()
    # loop through years
    for year, counter1 in zip(years, np.arange(0, len(years), 1)):
        # std and EmFix O3 for year of interest
        ffigac2_o3_ty = ffigac2_o3[counter1]
        emfix_o3_ty = emfix_o3[counter1]
        # std and EmFix NOx for year of interest
        ffigac2_nox_ty = ffigac2_no[counter1] + ffigac2_no2[counter1]
        emfix_nox_ty = emfix_no[counter1] + emfix_no2[counter1]
        # find difference between NOx and O3 in standard and fixed emissions 
        # simulations
        dO3 = emfix_o3_ty - ffigac2_o3_ty
        dNOx = emfix_nox_ty - ffigac2_nox_ty
        # linear regression
        dO3dNOx = np.polyfit(dNOx, dO3, 1)
        # plot data points
        axs[counter1].scatter(np.hstack(emfix_nox_ty - ffigac2_nox_ty), 
                              np.hstack(emfix_o3_ty - ffigac2_o3_ty), 
                              s = 5, c = 'k')
        # plot linear regression
        x = np.arange(dNOx.min(), dNOx.max(), 0.01)    
        axs[counter1].plot(x, dO3dNOx[0] * x + dO3dNOx[1], '-r')    
        # plot slope on subplot 
        axs[counter1].text(0.50, 0.75, '%.1f' %dO3dNOx[0], color = 'r',
                transform = axs[counter1].transAxes, fontsize = 12)
        # axis labels    
        axs[counter1].set_title(years[counter1])
        if (counter1 == 8) or (counter1 == 9) or (counter1 == 10) or (counter1 == 11):
            axs[counter1].set_xlabel('ppbv $\Delta$NO$_{x}$')    
        if (counter1 == 0) or (counter1 == 4) or (counter1 == 8):
            axs[counter1].set_ylabel('ppbv $\Delta$O$_{3}$')        
    axs[-1].axis('off')
    plt.subplots_adjust(wspace = 0.4, hspace = 0.7)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_deltanoxdeltao3_byyear_%d-%d_%s.eps' 
                %(years[0], years[-1], region), dpi = 300)
# # # # # # # # # # # # #    
def timeseries_deltanoxdeltao3_percent(ffigac2_no, emfix_no, mr2_no, egu_no, 
                                       ffigac2_no2, emfix_no2, mr2_no2, 
                                       egu_no2, ffigac2_o3, emfix_o3, mr2_o3, 
                                       egu_o3, t2m, years, region): 
    """function plots a timeseries of dO3/dNOx, calculating each of these 
    differences as a percentage difference (i.e. (new - old)/old). This is
    done for the (1) 2000 - 2010 Std and EmFix simulations from Strode et al.
    [2015], (2) the 2008 - 2010 MR2 and MR2+CEMS simulations, and (3) the 
    2008 - 2010 MR2 and MR2+CEMS simulations on hot days. 

    Parameters
    ----------   
    ffigac2_no : numpy.ndarray
        GMI CTM NO concentrations regionally-averaged over co-located (or nearly 
        colocated) with  corresponding CASTNet stations from Strode et al. 
        [2015] 'Std.' run (i.e. HindcastFFIgac2), units of parts per billion 
        volume, [11, days in months in 
        'sampling_months']   
    emfix_no : numpy.ndarray
        GMI CTM NO concentrations regionally-averaged over co-located (or nearly 
        colocated) with  corresponding CASTNet stations from Strode et al. 
        [2015] 'EmFix' run (i.e. Hindcast3Igac2), units of parts per billion 
        volume, [11, days in months in 
        'sampling_months']        
    mr2_no : numpy.ndarray
        GMI CTM NO concentrations regionally-averaged over co-located (or nearly 
        colocated) with  corresponding CASTNet stations HindcastMR2 simulation, 
        units of parts per billion volume, [3, days in months in 
        'sampling_months']     
    egu_no : numpy.ndarray        
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations from perturbed NO emissions simulation, 
        units of parts per billion volume, [3, stations in 'castnet_sites_fr', 
        days in months in 'sampling_months']            
    ffigac2_no2 : numpy.ndarray 
        Same as 'ffigac2_no' but for NO2   
    emfix_no2 : numpy.ndarray
        Same as 'emfix_no' but for NO2 
    mr2_no2 : numpy.ndarray
        Same as 'mr2_no' but for NO2
    egu_no2 : numpy.ndarray
        same as 'egu_no' but for NO2        
    ffigac2_o3 : numpy.ndarray
        Same as 'ffigac2_no' but for O3    
    emfix_o3 : numpy.ndarray
        Same as 'emfix_no' but for O3  
    mr2_o3 : numpy.ndarray
        Same as 'mr2_no' but for O3    
    egu_o3 : numpy.ndarray
        Same as 'egu_no' but for O3     
    t2m : numpy.ndarray
        Regionally-averaged MERRA 2 meter temperatures co-located (or nearly 
        colocated) with corresponding CASTNet stations, units of K, [3, 
        days in months in 'sampling_months']    
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
    import matplotlib.lines as mlines
    # initialize figure, axis
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    # calculate NOx (= NO + NO2)
    egu_nox = egu_no + egu_no2
    mr2_nox = mr2_no + mr2_no2
    emfix_nox = emfix_no + emfix_no2
    ffigac2_nox = ffigac2_no + ffigac2_no2  
    # loop through years in measuring period of EmFix/Std simulations
    for year, counter1 in zip(years, np.arange(0, len(years), 1)):
        # std and EmFix O3 for year of interest
        ffigac2_o3_ty = ffigac2_o3[counter1]
        emfix_o3_ty = emfix_o3[counter1]
        # std and EmFix NOx for year of interest
        ffigac2_nox_ty = ffigac2_nox[counter1]
        emfix_nox_ty = emfix_nox[counter1]
        # find percentage change between NOx and O3 in standard and fixed 
        # emissions simulations; i.e. % change = (perturbed - control)/control
        dO3 = (emfix_o3_ty - ffigac2_o3_ty)/ffigac2_o3_ty
        dNOx = (emfix_nox_ty - ffigac2_nox_ty)/ffigac2_nox_ty
        # calculate dO3/dNOx for percentage change quantities
        dO3dNOx = np.polyfit(dNOx, dO3, 1)
        ax.plot(years[counter1], dO3dNOx[0], 'ko', 
                label = 'EmFix-Std')
    # loop through years in measuring period of MR2/CEMS simulations
    for counter2 in np.arange(0, len([2008, 2009, 2010]), 1):
        egu_o3_ty = egu_o3[counter2]
        mr2_o3_ty = mr2_o3[counter2]
        egu_nox_ty = egu_nox[counter2]
        mr2_nox_ty = mr2_nox[counter2]
        dO3 = (egu_o3_ty - mr2_o3_ty)/mr2_o3_ty
        dNOx = (egu_nox_ty - mr2_nox_ty)/mr2_nox_ty
        dO3dNOx = np.polyfit(dNOx, dO3, 1)
        ax.plot([2008, 2009, 2010][counter2], dO3dNOx[0], 'b*', zorder = 20)
    # identify hot days for MR2-CEMS simulation and find dO3/dNOx on those 
    # days 
    for counter3 in np.arange(0, len([2008, 2009, 2010]), 1):    
        t2m_ty = t2m[counter3]
        hotdays = np.where(t2m_ty > np.percentile(t2m_ty, 80))[0]
        dO3_hot = (egu_o3_ty[hotdays] - mr2_o3_ty[hotdays])/mr2_o3_ty[hotdays]
        dNOx_hot = (egu_nox_ty[hotdays] - mr2_nox_ty[hotdays])/mr2_nox_ty[hotdays]
        dO3dNOx_hot = np.polyfit(dNOx_hot, dO3_hot, 1)
        ax.plot([2008, 2009, 2010][counter3], dO3dNOx_hot[0], 'rx', markersize = 8)
    ax.set_xlabel('Year')
    ax.set_ylabel('%$\Delta$ [O$_{3}$] %$\Delta$ [NO$_{x}$]$^{-1}$')
    # manually add items to legend
    ffigac2_line = mlines.Line2D([], [], color = 'k', linestyle = '',  
                                 marker = 'o', label = 'EmFix - Std')
    mr2_line = mlines.Line2D([], [], color = 'b', linestyle = '', marker = '*',
                             label = 'CEMS - MR2')
    mr2_hot_line = mlines.Line2D([], [], color = 'r', linestyle = '', 
                                 marker = 'x', markersize = 8, 
                                 label = '(CEMS - MR2)$_{\mathregular{hot}}$')
    plt.subplots_adjust(bottom = 0.3)
    leg = ax.legend(handles = [ffigac2_line, mr2_line, mr2_hot_line],
                     ncol = 3, bbox_to_anchor = (1.0, -0.2))
    leg.get_frame().set_linewidth(0.0)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'timeseries_deltanoxdeltao3_percent_%s.eps' 
                %(region), dpi = 300)
    return 
# # # # # # # # # # # # # 
def timeseries_noemiss_t2m(egu_no_emiss, mr2_no_emiss, t2m, years, 
                           region):
    """this function is similar to function 'NO_inventory_atpoint' in script 
    't2mdependence.py.' For the measuring period JJA 2008-2010 function plots 
    the regionally-averaged MR2 (control) NO_ff emission inventory and the 
    perturbed ("EGU") NO_ff emission inventory alongside regionally-averaged
    MERRA 2 meter temperatures (top subplot). The percentage change of the two 
    NO_ff emission inventories is also plotted (bottom subplot).
    
    Parameters
    ----------   
    egu_no_emiss : numpy.ndarray
        Regionally-averaged NO_ff emissions inventory from HindcastMR2 
        simulation, units of mol/mol, (3, 92)
    emfix_no : numpy.ndarray
        Regionally-averaged NO_ff emissions inventory from "EGU" simulation
        in which NO_ff emissions are daily-varying and scale with surface 
        temperatures, units of mol/mol, (3, 92)
    t2m : numpy.ndarray
        Regionally-averaged MERRA 2 meter temperatures co-located (or nearly 
        colocated) with corresponding CASTNet stations, units of K, (3, 92)
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
    # calculate percentage change in emissions (i.e. the percentage difference 
    # between the perturbed vs. monthly-mean emission inventories)
    no_percent_change = ((np.hstack(egu_no_emiss) - np.hstack(mr2_no_emiss))/
                         np.hstack(mr2_no_emiss)) * 100.    
    # initialize axes
    ax1 = plt.subplot2grid((2, 2), (0, 0), rowspan = 1, colspan = 2)
    ax2 = plt.subplot2grid((2, 2), (1, 0), rowspan = 1, colspan = 2)
    # plot perturbed emissions 
    ax1.plot(np.hstack(egu_no_emiss), linewidth = 1.5, clip_on = True, 
             label = 'Perturbed', color = '#1b9e77', linestyle = ':')
    # plot unperturbed emissions (n.b. loop through indices of months)
    month_idx = [0, 30, 61, 92, 122, 153, 184, 214, 245, 276]
    for a, b in zip(month_idx[:-1], month_idx[1:]):
        ax1.plot(np.arange(a, b, 1), np.hstack(mr2_no_emiss)[a:b], 
                 linewidth = 2., clip_on = True,  color = '#1b9e77')    
    # twin axis, plot temperature 
    ax1t = ax1.twinx()
    ax1t.plot(np.hstack(t2m), linewidth = 1.5, color = '#d95f02', clip_on = True)          
    # plot percentage change of emissions         
    ax2.plot(no_percent_change, linewidth = 1.5, clip_on = True,  color = 'k', 
             linestyle = '-')    
    # set axis limits and ticks
    ax1.set_xlim([0, len(np.hstack(egu_no_emiss))])
    ax2.set_xlim([0, len(np.hstack(egu_no_emiss))])
    ax1.set_xticks(month_idx[:-1])
    ax1.set_xticklabels([])
    ax2.set_xticks(month_idx[:-1])
    ax2.set_xticklabels(['1 June 2008', '', '', '1 June 2009', '', '', 
                         '1 June 2010', '', ''])
    ax2.set_yticks([-20, -15, -10, -5, 0, 5, 10, 15, 20])
    ax2.set_yticklabels(['-20', '', '-10', '', '0', '', '10', '', '20'])    
        
    # add axis labels and color ticks labels
    ax1.set_ylabel('Emiss$_{\mathregular{NO}}$ [mol mol$^{-1}$]', 
                   color = '#1b9e77')           
    ax1t.set_ylabel('T$_{\mathregular{2m}}$ [K]', color = '#d95f02')
    for label in ax1.get_yticklabels():
        label.set_color('#1b9e77')
    for label in ax1t.get_yticklabels():
        label.set_color('#d95f02')
    ax2.set_ylabel('$\Delta$ Emiss$_{\mathregular{NO}}$ [%]', 
                   color = 'k')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'timeseries_noemiss_t2m_%d-%d_%s.eps' 
                %(years[0], years[-1], region), dpi = 300)
    return 
# # # # # # # # # # # # #
def scatter_dEmissdNOxdO3(ffigac2_no_emiss, mr2_no_emiss, egu_no_emiss, 
                          ffigac2_no, emfix_no, ffigac2_no2, emfix_no2, 
                          mr2_no, egu_no, mr2_no2, egu_no2, ffigac2_o3, 
                          emfix_o3, mr2_o3, egu_o3, years, years_strode, cems, 
                          ffigac2_castnet, t2m):
    """function finds the relative difference in NO and O3 between Std and 
    EmFix simulations discussed in Strode et al. (2016) and the relative 
    difference between NO and O3 from the MR2-CEMS simulations along with the
    relative variability of NO in these simulations. Relative differences are 
    calculated for each summer (JJA) season along with each summer month. 
    These differences are plotted in scatterplot along with the relative 
    difference in NOx from CEMS (i.e. CEMS for a particular year - CEMS from 
    2000) and CASTNet O3 (i.e. CASTNet for a particular year - 2000-2003 
    average CASTNet).
    
    Parameters
    ----------   
    ffigac2_no_emiss : numpy.ndarray
        Regionally-averaged monthly-varying NO from fossil fuels from the 
        emissions inventory from HindcastFFIgac2 simulation, units of mol/mol, 
        (years in 'years_strode', 92)
    mr2_no_emiss : numpy.ndarray
        Regionally-averaged monthly-varying NO from fossil fuels from the 
        emissions inventory from HindcastMR2 simulation, units of mol/mol, 
        (years in 'years', 92)
    egu_no_emiss : numpy.ndarray
        Regionally-averaged daily-varying NO from fossil fuels which scale with
        daily surface temperatures from the emissions inventory from 
        HindcastMR2 simulation, units of mol/mol, (years in 'years', 92)
    ffigac2_no : numpy.ndarray
        Regionally-averaged NO from HindcastFFIgac2, units of ppbv, 
        (years in 'years_strode', 92)    
    emfix_no : numpy.ndarray
        Regionally-averaged NO from Hindcast3Igac2, units of ppbv, 
        (years in 'years_strode', 92)    
    ffigac2_no2 : numpy.ndarrayh
        Regionally-averaged NO2 from HindcastFFIgac2, units of ppbv, 
        (years in 'years_strode', 92)    
    emfix_no2 : numpy.ndarray
        Regionally-averaged NO2 from Hindcast3Igac2, units of ppbv, 
        (years in 'years_strode', 92)        
    mr2_no : numpy.ndarray 
        Regionally-averaged NO from HindcastMR2, units of ppbv, 
        (years in 'years_strode', 92)        
    egu_no : numpy.ndarray
        Regionally-averaged NO from HindcastMR2 with daily-varying NO 
        emissions, units of ppbv, (years in 'years', 92)      
    mr2_no2 : numpy.ndarray
        Regionally-averaged NO2 from HindcastMR2, units of ppbv, 
        (years in 'years_strode', 92)            
    egu_no2 : numpy.ndarray
        Regionally-averaged NO2 from HindcastMR2 with daily-varying NO 
        emissions, units of ppbv, (years in 'years', 92)          
    ffigac2_o3 : numpy.ndarray
        Regionally-averaged O3 from HindcastFFIgac2, units of ppbv, 
        (years in 'years_strode', 92)
    emfix_o3 : numpy.ndarray
        Regionally-averaged O3 from Hindcast3Igac2, units of ppbv, 
        (years in 'years_strode', 92)    
    mr2_o3 : numpy.ndarray
        Regionally-averaged O3 from HindcastMR2, units of ppbv, (years in 
        'years', 92)        
    egu_o3 : numpy.ndarray
        Regionally-averaged O3 from HindcastMR2 with daily-varying NO 
        emissions, units of ppbv, (years in 'years', 92)            
    years : list
        Years in measuring period for MR2-CEMS simulations
    years_strode : list
        Years in measuring period for Strode et al. (2016) simulations
    cems : pandas.core.series.Series
        Summertime (JJA) CEMS NOx emissions from the 14 states in the 
        Northeastern United States for the measuring period 2000-2014, units of
        tons, (1380,)
    ffigac2_castnet : numpy.ndarray
        Regionally-averaged CASTNet O3, (years in 'years_strode', 92)
    t2m : numpy.ndarray
        Regionally-averaged MERRA 2 meter temperatures co-located (or nearly 
        colocated) with corresponding CASTNet stations, units of K, (3, 92)

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
    # create lists to be filled with relative changes of NO emissions, NOx, 
    # and O3
    no_emiss_monthly, no_emiss = [], []
    nox_monthly, nox = [], [] 
    o3_monthly, o3 = [], []
    no_emiss_strode_monthly, no_emiss_strode = [], []
    nox_strode_monthly, nox_strode = [], []
    o3_strode_monthly, o3_strode = [], []
    #cems_nox, cems_nox_mean = [], []
    #castnet_o3, castnet_o3_mean = [], []
    # first day of month indices
    month_idx = [0, 30, 61, 92]
    # for "EmFix" run in Strode et al. (2015), emissions were kept fixed at 
    # 2000 levels 
    emfix_no_emiss = ffigac2_no_emiss[0]
    #cems_2000 = cems.loc['2000-06-01':'2000-08-31'].values
    #castnet_2000 = ffigac2_castnet[0:3]
    # loop through 2000-2010 measuring period of Strode et al. (2015) 
    for year in years_strode:
        # find position of year in year list  
        ty = np.where(np.array(years_strode) == year)[0][0]
        # NO emissions from year of interest 
        ffigac2_no_emiss_ty = ffigac2_no_emiss[ty]
        # NOx from year of interest
        ffigac2_nox_ty = ffigac2_no[ty] + ffigac2_no2[ty]
        emfix_nox_ty = emfix_no[ty] + emfix_no2[ty]
        # O3 from year of interest
        ffigac2_o3_ty = ffigac2_o3[ty]
        emfix_o3_ty = emfix_o3[ty]    
        # relative changes in O3 and NO for each summer month
        for start, stop in zip(month_idx[:-1], month_idx[1:]):
            no_emiss_strode_monthly.append(np.mean(emfix_no_emiss[start:stop-1] - 
                ffigac2_no_emiss_ty[start:stop-1]))
            nox_strode_monthly.append(np.mean(emfix_nox_ty[start:stop-1] - 
                    ffigac2_nox_ty[start:stop-1]))
            o3_strode_monthly.append(np.mean(emfix_o3_ty[start:stop-1] - 
                                     ffigac2_o3_ty[start:stop-1]))
        # relative changes in O3 and NO for season
        no_emiss_strode.append((emfix_no_emiss - ffigac2_no_emiss_ty))    
        nox_strode.append((emfix_nox_ty - ffigac2_nox_ty))
        o3_strode.append((emfix_o3_ty - ffigac2_o3_ty))
        ## CEMS NOx emissions for year
        #cems_ty = cems.loc['%s-06-01'%year:'%s-08-31'%year].values
        ## CASTNet O3
        #castnet_ty = ffigac2_castnet[ty]
        ## relative changes in emissions from CEMS
        #cems_nox.append(cems_2000 - cems_ty)
        #cems_nox_mean.append(np.mean(cems_2000 - cems_ty))
        ## relative changes in CASTNet O3
        #castnet_o3.append(castnet_2000 - castnet_ty)
        #castnet_o3_mean.append(np.mean(castnet_2000 - castnet_ty))    
    # loop through 2008-2010 measuring period of MR2-CEMS simulations
    for year in years:
        # find position of year in list 'years'
        ty = np.where(np.array(years) == year)[0][0]
        # T2m for year of interest
        t2m_ty = t2m[ty]
        # NO emissions from MR2-CEMS simulations for year of interest
        egu_no_emiss_ty = egu_no_emiss[ty]
        mr2_no_emiss_ty = mr2_no_emiss[ty]
        # NOx from year of interest
        mr2_nox_ty = mr2_no[ty] + mr2_no2[ty]
        egu_nox_ty = egu_no[ty] + egu_no2[ty]
        # O3 from year of interest
        mr2_o3_ty = mr2_o3[ty]
        egu_o3_ty = egu_o3[ty]
        # find relative difference in NO emissions, NOx, and O3 for each month
        for start, stop in zip(month_idx[:-1], month_idx[1:]):
            # find hot and cold days in each month, i.e. top and lowest 20th %-ile 
            hot = np.where(t2m_ty[start:stop-1] >= 
                           np.percentile(t2m_ty[start:stop-1], 80))[0]
            cold = np.where(t2m_ty[start:stop-1] <= 
                            np.percentile(t2m_ty[start:stop-1], 20))[0]
            emiss_hot = (egu_no_emiss_ty[start:stop-1][hot] - 
                         mr2_no_emiss_ty[start:stop-1][hot])
            emiss_cold = (egu_no_emiss_ty[start:stop-1][cold] - 
                          mr2_no_emiss_ty[start:stop-1][cold])
            nox_hot = (egu_nox_ty[start:stop-1][hot] - 
                       mr2_nox_ty[start:stop-1][hot])
            nox_cold = (egu_nox_ty[start:stop-1][cold] - 
                        mr2_nox_ty[start:stop-1][cold])
            o3_hot = (egu_o3_ty[start:stop-1][hot] - 
                      mr2_o3_ty[start:stop-1][hot])
            o3_cold = (egu_o3_ty[start:stop-1][cold] - 
                       mr2_o3_ty[start:stop-1][cold])
            # append monthly relative differences to lists
            no_emiss_monthly.append(np.mean(emiss_hot - emiss_cold))
            nox_monthly.append(np.mean(nox_hot - nox_cold))
            o3_monthly.append(np.mean(o3_hot - o3_cold))
        # find hot and cold days for each season
        hot = np.where(t2m_ty >= np.percentile(t2m_ty, 80))[0]
        cold = np.where(t2m_ty <= np.percentile(t2m_ty, 20))[0]
        emiss_hot = egu_no_emiss_ty[hot] - mr2_no_emiss_ty[hot]
        emiss_cold = egu_no_emiss_ty[cold] - mr2_no_emiss_ty[cold]
        nox_hot = egu_nox_ty[hot] - mr2_nox_ty[hot]
        nox_cold = egu_nox_ty[cold] - mr2_nox_ty[cold]
        o3_hot = egu_o3_ty[hot] - mr2_o3_ty[hot]
        o3_cold = egu_o3_ty[cold] - mr2_o3_ty[cold]
        o3_hot_mean = (egu_o3_ty[hot]) - (mr2_o3_ty[hot])
        o3_cold_mean = (egu_o3_ty[cold]) - (mr2_o3_ty[cold])            
        # append seasonal relative differences to lists
        no_emiss.append((emiss_hot) - (emiss_cold))
        nox.append((nox_hot) - (nox_cold))
        o3.append(o3_hot_mean - o3_cold_mean)
    COLOR_FFIGAC2 = '#e41a1c'
    COLOR_MR2 = '#377eb8'
    # # # # for dO3/dEmiss
    ax = plt.subplot2grid((1, 1), (0, 0))
    # for Strode et al. (2015) simulation daily values 
    ax.plot(no_emiss_strode, o3_strode, 'o', markersize = 1, color = COLOR_MR2)        
    # for MR2-CEMS simulation daily values 
    ax.plot(no_emiss, o3, 'o', markersize = 1, color = COLOR_MR2)      
    # for Strode et al. (2015) simulation monthly values 
    ax.plot(no_emiss_strode_monthly, o3_strode_monthly, 's', markersize = 8, 
            markerfacecolor = 'w', markeredgecolor = COLOR_FFIGAC2)
    ax.plot(no_emiss_strode_monthly, o3_strode_monthly, 'o', markersize = 3, 
            color = COLOR_FFIGAC2, label = 'Strode et al. (2015)')    
    # for MR2-CEMS simulation 
    ax.plot(no_emiss_monthly, o3_monthly, 's', markersize = 8, markerfacecolor = 'w', 
            markeredgecolor = COLOR_MR2)  
    ax.plot(no_emiss_monthly, o3_monthly, 'o', markersize = 3, color = COLOR_MR2, 
            label = 'Daily-varying NO')        
    # fit linear regression through seasonal means
    x = np.hstack((no_emiss_monthly, no_emiss_strode_monthly))
    y = np.hstack((o3_monthly, o3_strode_monthly))
    fit = np.polyfit(x, y, 1)
    # plot linear regression
    fit_fn = np.poly1d(fit) 
    # fit_fn is now a function which takes in x and returns an estimate for y
    ax.plot(x, fit_fn(x), '--k', zorder = 1)
    ax.text(0.05, 0.9, 
            r'$\frac{\mathregular{\Delta O_{\mathregular{3}}}}' \
            '{\mathregular{\Delta} \mathregular{Emiss}_{\mathregular{NO}}}$ = %.2f' %fit[0], 
            color = 'k', transform = ax.transAxes, fontsize = 14)
    # axis labels 
    ax.set_xlabel('$\mathregular{\Delta}$ Emiss$_{\mathregular{NO}}$ [kg s$^{-1}$ grid cell$^{-1}$]')
    ax.set_ylabel('$\mathregular{\Delta}$ O$_{\mathregular{3}}$ [ppbv]')
    plt.legend(loc = 'lower right', frameon = False)
    ## twin y-axis and add dO3/dNOx from CASTNet/CEMS
    #axt = ax.twiny()
    #axt.plot(cems_nox, castnet_o3, 'ko', label = 'CEMS-CASTNet')
    #axt.set_xlabel('$\mathregular{\Delta}$ CEMS NO$_{x}$ [tons]')
    ## remove old legend, create new
    #ax.legend_.remove()
    #lines, labels = ax.get_legend_handles_labels()
    #lines2, labels2 = axt.get_legend_handles_labels()
    #plt.legend(lines + lines2, labels + labels2, loc = 'lower right', 
    #           frameon = False)
    ## save figure WITH CEMS NOx/CASTNet O3 
    #plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
    #            'scatter_deltanoxdeltao3_withCEMS.eps', dpi = 300)    
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'scatter_dEmissdO3.eps', dpi = 300)
    # # # # for dNOx/dEmiss
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.plot(no_emiss_strode, nox_strode, 'o', markersize = 1, 
            color = COLOR_FFIGAC2)    
    ax.plot(no_emiss, nox, 'o', markersize = 1, color = COLOR_MR2)  
    ax.plot(no_emiss_strode_monthly, nox_strode_monthly, 's', markersize = 8, 
            markerfacecolor = 'w', markeredgecolor = COLOR_FFIGAC2)
    ax.plot(no_emiss_strode_monthly, nox_strode_monthly, 'o', markersize = 3, 
            color = COLOR_FFIGAC2, label = 'Strode et al. (2015)')
    ax.plot(no_emiss_monthly, nox_monthly, 's', markersize = 8, markerfacecolor = 'w', 
            markeredgecolor = COLOR_MR2)  
    ax.plot(no_emiss_monthly, nox_monthly, 'o', markersize = 3, color = COLOR_MR2, 
            label = 'Daily-varying NO')    
    x = np.hstack((no_emiss_monthly, no_emiss_strode_monthly))
    y = np.hstack((nox_monthly, nox_strode_monthly))
    fit = np.polyfit(x, y, 1)
    fit_fn = np.poly1d(fit) 
    ax.plot(x, fit_fn(x), '--k', zorder = 1)
    ax.text(0.05, 0.9, 
            r'$\frac{\mathregular{\Delta NO_{x}}}' \
            '{\mathregular{\Delta} \mathregular{Emiss}_{\mathregular{NO}}}$ = %.2f' %fit[0], 
            color = 'k', transform = ax.transAxes, fontsize = 14)
    ax.set_xlabel('$\mathregular{\Delta}$ Emiss$_{\mathregular{NO}}$ [kg s$^{-1}$ grid cell$^{-1}$]')
    ax.set_ylabel('$\mathregular{\Delta}$ NO$_{x}$ [ppbv]')
    plt.legend(loc = 'lower right', frameon = False)        
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'scatter_dEmissdNOx.eps', dpi = 300)
    # # # # for dO3/dNOx
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.plot(nox_strode, o3_strode, 'o', markersize = 1, color = COLOR_FFIGAC2)
    ax.plot(nox, o3, 'o', markersize = 1, color = COLOR_MR2)    
    ax.plot(nox_strode_monthly, o3_strode_monthly, 's', markersize = 8, 
            markerfacecolor = 'w', markeredgecolor = COLOR_FFIGAC2)
    ax.plot(nox_strode_monthly, o3_strode_monthly, 'o', markersize = 3, 
            color = COLOR_FFIGAC2, label = 'Strode et al. (2015)')    
    ax.plot(nox_monthly, o3_monthly, 's', markersize = 8, markerfacecolor = 'w', 
            markeredgecolor = COLOR_MR2)  
    ax.plot(nox_monthly, o3_monthly, 'o', markersize = 3, color = COLOR_MR2, 
            label = 'Daily-varying NO')    
    x = np.hstack((nox_monthly, nox_strode_monthly))
    y = np.hstack((o3_monthly, o3_strode_monthly))
    fit = np.polyfit(x, y, 1)
    fit_fn = np.poly1d(fit) 
    ax.plot(x, fit_fn(x), '--k', zorder = 1)
    ax.text(0.05, 0.9, 
            r'$\frac{\mathregular{\Delta O_{\mathregular{3}}}}' \
            '{\mathregular{\Delta} \mathregular{NO}_{x}}$ = %.2f' %fit[0], 
            color = 'k', transform = ax.transAxes, fontsize = 14)
    ax.set_xlabel('$\mathregular{\Delta}$ NO$_{x}$ [ppbv]')
    ax.set_ylabel('$\mathregular{\Delta}$ O$_{\mathregular{3}}$ [ppbv]')    
    plt.legend(loc = 'lower right', frameon = False)    
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'scatter_dNOxdO3.eps', dpi = 300)
    return
# # # # # # # # # # # # #    
def pdf_controlo3eguo3(egu_o3, mr2_o3, mr2_castnet, years, region):
    """function plots probability density functions (PDFs) of regionally-
    averaged O3 from HindcastMR2 and daily-varying NO simulations alongside 
    CASTNet O3 observations. 
    
    Parameters
    ----------       
    egu_o3 : numpy.ndarray
        Regionally-averaged O3 from HindcastMR2 with daily-varying NO 
        emissions, units of ppbv, (years in 'years', 92) 
    mr2_castnet : numpy.ndarray
        Regionally-averaged O3 observations from CASTNet, units of ppbv, 
        (years in 'years', 92)
    mr2_o3 : numpy.ndarray
        Regionally-averaged O3 from HindcastMR2, units of ppbv, (years in 
        'years', 92)        
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
    import matplotlib as mpl
    import matplotlib.mlab as mlab     
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
    COLOR_MR2 = '#e41a1c'
    COLOR_NO = '#337eb8'
    # initialize figure, axis
    ax = plt.subplot2grid((1, 1), (0, 0), rowspan = 1, colspan = 1)
    # for Daily Varying NO simultion 
    n, bins, patches = ax.hist(np.hstack(egu_o3), 25, fc = 'r', ec = 'r', 
                               density = 1, alpha = 0.)
    pdf_egu = mlab.normpdf(bins, np.mean(egu_o3), np.std(egu_o3))
    ax.plot(bins, pdf_egu, lw = 1.5, color = COLOR_NO, ls = '-', label = 'Daily Varying NO')
    # for HindcastMR2 simulation
    n, bins, patches = ax.hist(np.hstack(mr2_o3), 25, fc = 'r', ec = 'r', 
                               density = 1, alpha = 0.)
    pdf_mr2 = mlab.normpdf(bins, np.mean(mr2_o3), np.std(mr2_o3))
    ax.plot(bins, pdf_mr2, lw = 1.5, color = COLOR_MR2, ls = '-', label = 'HindcastMR2')
    # for CASTNet        
    n, bins, patches = ax.hist(np.hstack(mr2_castnet), 25, fc = 'r', ec = 'r', 
                               density = 1, alpha = 0.)
    pdf_castnet = mlab.normpdf(bins, np.mean(mr2_castnet), np.std(mr2_castnet))
    ax.plot(bins, pdf_castnet, lw = 1.5, color = 'k', ls = '-', label = 'CASTNet')
    # add legend and titles
    leg = ax.legend(bbox_to_anchor = (1.02, 0.65))
    leg.get_frame().set_linewidth(0.0)
    ax.set_xlabel('Ozone [ppbv]')
    ax.set_ylabel('Probability Density')
    plt.subplots_adjust(right = 0.70, bottom = 0.15, left = 0.15)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'pdf_controlo3eguo3_%d-%d_%s.eps' %(years[0], years[-1], 
                                                    region), dpi = 300)
    return 
# # # # # # # # # # # # # #    
def timeseries_controlo3eguo3(egu_o3, mr2_o3, year, years, region):
    """function plots timeseries of regionally-averaged O3 from HindcastMR2 and
    HindcastMR2 with daily-varying NO emissions simulations. 
    
    Parameters
    ----------       
    egu_o3 : numpy.ndarray
        Regionally-averaged O3 from HindcastMR2 with daily-varying NO 
        emissions, units of ppbv, (years in 'years', 92) 
    mr2_o3 : numpy.ndarray
        Regionally-averaged O3 from HindcastMR2, units of ppbv, (years in 
        'years', 92)
    year : int
        Year for which time series will be plotted        
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
    COLOR_MR2 = '#e41a1c'
    COLOR_NO = '#337eb8'
    ty = np.where(np.array(years) == year)[0][0]
    # O3 from simulations for year of interest
    egu_o3_ty = egu_o3[ty]
    mr2_o3_ty = mr2_o3[ty]
    # initialize figure, axis
    ax = plt.subplot2grid((1, 1), (0, 0), rowspan = 1, colspan = 1)
    ax.plot(egu_o3_ty, lw = 1.5, ls = '-', color = COLOR_NO)
    ax.plot(mr2_o3_ty, lw = 1.5, ls = '-', color = COLOR_MR2)
    ax.set_xlim([0, len(egu_o3_ty)])
    ax.set_xticks([0, 30, 61])
    ax.set_xticklabels(['1 June %s'%year, '1 July %s'%year, '1 August %s'%year])
    ax.set_ylabel('O$_{3}$ [ppbv]')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' +
                'timeseries_controlo3eguo3_%s.eps' %(year), dpi = 300)    
    return 
# # # # # # # # # # # # #    
import numpy as np
import sys
sys.path.append('/Users/ghkerr/phd/GMI/')
import commensurability
from calculate_regional_average import calculate_regional_average
import sys
sys.path.append('/Users/ghkerr/phd/emissions/')
import AQSCEMSobs
years = [2008, 2009, 2010]            
years_strode = [2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010]
sampling_months = [6, 7, 8]
sampling_hours = [15, 16, 17, 18, 19, 20]  
castnet_sites_fr = ['ASH', 'WFM', 'WST', 'APT', 'SAR', 'CTH', 'WSP',
                    'ARE', 'BEL', 'PSU', 'LRL', 'PAR', 'PED', 'SHN', 
                    'CDR', 'VPI', 'MKG', 'KEF']
region = 'northeast'
# # # # open CTM output # # # #
# open HindcastMR2 run
mr2_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2', 
                                          years, sampling_months, 
                                          sampling_hours)
# open HindcastMR2 run with temperature-dependent emissions 
egu_castnet, egu_o3, egu_no, egu_no2, egu_co, egu_gmi_sites_fr = \
commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'EGU_T', years, 
                                          sampling_months, sampling_hours)    
# open HindcastFFIgac2 (emissions variable) run 
ffigac2_castnet, ffigac2_o3, ffigac2_no, ffigac2_no2, ffigac2_co, ffigac2_gmi_sites_fr = \
commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'HindcastFFIgac2', 
                                          years_strode, sampling_months, 
                                          sampling_hours)
# open Hindcast3Igac2 (emissions fixed at 2000 levels) run 
emfix_castnet, emfix_o3, emfix_no, emfix_no2, emfix_co, emfix_gmi_sites_fr = \
commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'Hindcast3Igac2', 
                                          years_strode, sampling_months, 
                                          sampling_hours)
# # # # open emissions inventories (NO from fossil fuel sector) # # # #
# for HindcastMR2 run
mr2_no_emiss = \
commensurability.commensurate_emiss_unperturbed(mr2_castnet, castnet_sites_fr, 
                                                years, sampling_months, 
                                                '1x1.25_IAVanthGFED4')
# for HindcastMR2 run with temperature-dependent emissions run
egu_no_emiss, emiss_lats_fr, emiss_lons_fr = \
commensurability.commensurate_emiss_perturbed(mr2_castnet, castnet_sites_fr, 
                                              years, sampling_months)
# for FFIgac2 run 
ffigac2_no_emiss = \
commensurability.commensurate_emiss_unperturbed(ffigac2_castnet, castnet_sites_fr, 
                                                years_strode, sampling_months, 
                                                '2x2.5_IAVanthGFED3gcEF')
# # # # open 2 meter temperature from MERRA # # # #
t2m, merra_lats_fr, merra_lons_fr = \
commensurability.commensurate_t2m(mr2_castnet, castnet_sites_fr, years, 
                                  sampling_months)
# # # # calculate regional averages # # # #
# for HindcastMR2 run
mr2_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co = \
calculate_regional_average(mr2_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, 1)
# for HindcastMR2 run with temperature-dependent emissions run
egu_castnet, egu_o3, egu_no, egu_no2, egu_co = \
calculate_regional_average(egu_castnet, egu_o3, egu_no, egu_no2, egu_co, 1)
# for HindcastFFIgac2 (emissions variable) run 
ffigac2_castnet, ffigac2_o3, ffigac2_no, ffigac2_no2, ffigac2_co = \
calculate_regional_average(ffigac2_castnet, ffigac2_o3, ffigac2_no, 
                           ffigac2_no2, ffigac2_co, 1)
# for Hindcast3Igac2 (emissions fixed at 2000 levels) run 
emfix_castnet, emfix_o3, emfix_no, emfix_no2, emfix_co = \
calculate_regional_average(emfix_castnet, emfix_o3, emfix_no, emfix_no2, 
                           emfix_co, 1)
# for emission inventories
egu_no_emiss = np.nanmean(egu_no_emiss, axis = 1)
mr2_no_emiss = np.nanmean(mr2_no_emiss, axis = 1)
ffigac2_no_emiss = np.nanmean(ffigac2_no_emiss, axis = 1)
# for regionally-averaged 2 meter temperatures 
t2m = np.nanmean(t2m, axis = 1)
# # # # load CEMS NOx for NEUS # # # #
states_ab = ['CT', 'DC', 'DE', 'MA', 'MD', 'ME', 'NH', 'NJ', 'NY', 'PA', 
             'RI', 'VA', 'VT', 'WV']
cems, nox_lat, nox_lon = AQSCEMSobs.cems_specifystates_dailymean(
        '/Volumes/GAIGEKERR/emissions/CEMS/', states_ab, sampling_months)
 # # # visualizations for daily-average O3 # # # # 
scatter_controlnoeguno_bysite(mr2_no, egu_no, castnet_sites_fr, years, region)
scatter_controlno2egunno2_bysite(mr2_no2, egu_no2, castnet_sites_fr, years, region)
scatter_controlo3eguno3_bysite(mr2_o3, egu_o3, castnet_sites_fr, years, region)
scatter_inventory_ctm(egu_no_emiss, mr2_no_emiss, egu_no, 
                      mr2_no, egu_no2, mr2_no2, egu_o3, mr2_o3, years, region)
scatter_deltanoxdeltao3(egu_no, mr2_no, ffigac2_no, emfix_no, egu_no2, 
                        mr2_no2, ffigac2_no2, emfix_no2, egu_o3, mr2_o3,
                        ffigac2_o3, emfix_o3, years, region)
scatter_deltanoxdeltao3_percent(egu_no, mr2_no, ffigac2_no, emfix_no, 
                                egu_no2, mr2_no2, ffigac2_no2, emfix_no2, 
                                egu_o3, mr2_o3, ffigac2_o3, emfix_o3, 
                                years, region)
timeseries_noemiss_t2m(egu_no_emiss, mr2_no_emiss, t2m, years, region)    
scatter_dEmissdNOxdO3(ffigac2_no_emiss, mr2_no_emiss, egu_no_emiss, 
                      ffigac2_no, emfix_no, ffigac2_no2, emfix_no2, mr2_no, 
                      egu_no, mr2_no2, egu_no2, ffigac2_o3, emfix_o3, mr2_o3, 
                      egu_o3, years, years_strode, cems, ffigac2_castnet, t2m)
pdf_controlo3eguo3(egu_o3, mr2_o3, mr2_castnet, years, region)
timeseries_controlo3eguo3(egu_o3, mr2_o3, 2008, years, region)
# # # # open time- and regionally-averaged diurnal trace gases # # # # 
castnet_d, mr2_o3_d, mr2_no_d, mr2_no2_d, mr2_co_d = \
commensurability.commensurate_castnet_gmi_diurnal(castnet_sites_fr, 'HindcastMR2', 
                                                  years, sampling_months)
temp, egu_o3_d, egu_no_d, egu_no2_d, egu_co_d = \
commensurability.commensurate_castnet_gmi_diurnal(castnet_sites_fr, 'EGU_T', 
                                                  years, sampling_months)
# # # # calculate regional averages of diurnal curves # # # #
castnet_d, mr2_o3_d, mr2_no_d, mr2_no2_d, mr2_co_d = \
calculate_regional_average(castnet_d, mr2_o3_d, mr2_no_d, mr2_no2_d, mr2_co_d, 
                           1)
temp, egu_o3_d, egu_no_d, egu_no2_d, egu_co_d = \
calculate_regional_average(temp, egu_o3_d, egu_no_d, egu_no2_d, egu_co_d, 
                           1)

# # # # visualizations for diurnal curves
diurnal_controlegunono2nox(mr2_no_d, mr2_no2_d, egu_no_d, egu_no2_d, 
                           castnet_sites_fr, years, region)
diurnal_controlegunono2o3(mr2_no_d, egu_no_d, mr2_no2_d, egu_no2_d, 
                          mr2_o3_d, egu_o3_d, t2m, years, region)
# # # # open HindcastFFIgac2 (emissions variable) run 
temp, ffigac2_o3, ffigac2_no, ffigac2_no2, ffigac2_co, ffigac2_gmi_sites_fr = \
commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'HindcastFFIgac2', 
                                          list(np.arange(2000, 2011, 1)), 
                                          sampling_months, sampling_hours)
# # # # open Hindcast3Igac2 (emissions fixed at 2000 levels) run 
temp, emfix_o3, emfix_no, emfix_no2, emfix_co, emfix_gmi_sites_fr = \
commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'Hindcast3Igac2', 
                                          list(np.arange(2000, 2011, 1)), 
                                          sampling_months, sampling_hours)
del temp
# find regionally-averaged trace gas concentrations from simulations
# for NO
emfix_no = np.nanmean(emfix_no, axis = 1)
ffigac2_no = np.nanmean(ffigac2_no, axis = 1)
# for NO2
emfix_no2 = np.nanmean(emfix_no2, axis = 1)
ffigac2_no2 = np.nanmean(ffigac2_no2, axis = 1)
# for O3
emfix_o3 = np.nanmean(emfix_o3, axis = 1)
ffigac2_o3 = np.nanmean(ffigac2_o3, axis = 1)
# # # # visualizations 
scatter_deltanoxdeltao3_byyear(ffigac2_no, emfix_no, ffigac2_no2, emfix_no2, 
                               ffigac2_o3, emfix_o3, list(np.arange(2000, 2011, 1)), region)
timeseries_deltanoxdeltao3_percent(ffigac2_no, emfix_no, mr2_no, egu_no, 
                                   ffigac2_no2, emfix_no2, mr2_no2, egu_no2, 
                                   ffigac2_o3, emfix_o3, mr2_o3, egu_o3, t2m, 
                                   list(np.arange(2000, 2011, 1)), region)