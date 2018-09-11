# -*- coding: utf-8 -*-
"""
NAME
    gmi_visualizations.py
PURPOSE
    Compare hourly trace gas observations and modeled trace gases from GMI CTM
PROGRAMMER
    Gaige Hunter Kerr
REVISION HISTORY
    01092017 -- initial version created
    04092017 -- boxplot function added
    05092017 -- CDF function added, fixed sampling_hours bug (formerly it 
                was for 11-16 LST, should be in UTC)
    07092017 -- Some FFIgac2 sites have questionable O3 concentrations, 
                so, the sampling of sites for this case was hard coded in 
                functions
    08092017 -- JJA timeseries function added
    09142017 -- given error in FFIgac2 longitidue, FFIgac2_error_megan function
                added
    09222017 -- change cdf_castnet_gmi to include easier CDF formation and 
                additional model
    03102017 -- Regional and daily averaged O3 concentrations were incorrectly
                calculated for boxplot, CDF, and PDF functions (i.e. before
                the average was taken over the year dimension rather than the 
                hour dimension); this problem was corrected
    10102017 -- ctmreplay_map function added
    02012018 -- changed code so that I only include GMI output from sites in 
                regional average if they have a co-located CASTNET site with 
                data for a particular year
    03012018 -- function gmi_castnet_region added; function gmi_castnet_noregavg 
                edited so that for a particular GMI site all co-located 
                CASTNet sites' observations are fetched and averaged
    19012018 -- modified function open_castnet_gmi to also extract modeled 
                CO, NO, and NO2 and adapted function to handle perturbed
                emissions runs; function gmi_castnet_dailyavg_singleyear added
    16042018 -- functions 'map_gmiaqscastnet', 'diurnal_castneto3gmio3' and
                'diurnaal_aqscogmico' added
    30052018 -- function 'boxplot_castneto3gmio3' renamed 'boxplot_tracegas' 
                and modified so that it could plot other trace gases besides O3
    03062018 -- code to open files moved to separate file 
                'run_gmi_visualizations.py', this file renamed 
                'gmi_visualizations'
    04062018 -- removed np.roll from function 'diurnal_castneto3gmio3' as 
                code modifications converted CASTNet observations to UTC 
                prior to input in this function eliminating the need for this
    26072018 -- function 'pdf_allgmio3' added
    29072018 -- function 'scatter_allgmio3' added
    22082018 -- function 'timeseries_mr2o3dato3eguto3' added
    25082018 -- function 'scatter_dt2m_dmr2o3dato3eguto3' added
    28082018 -- edited function 'map_gmiaqscastnet' to prepare for publication
    29082018 -- function 'timeseries_t2m_castneto3_cemsnox' added
    04092018 -- functions 'map_simulationschematic' and 
                'contourf_merra2h500_mr2o3event' added
    05082018 -- functions 'map_do3dt2mratio_conus', 
                'timeseries_t2m_datnox_dato3', and 'scatter_datnox_dato3' added
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
def scatter_t2m_castneto3(comm_t2m, comm_castnet, years, region): 
    """function plots regionally-averaged CASTNet O3 vs. MERRA 2 meter 
    temperature and finds the linear regression

    Parameters
    ----------    
    comm_t2m : numpy.ndarray    
        Regionally-averaged MERRA 2 meter temperatures co-located (or nearly 
        colocated) with corresponding CASTNet stations, units of K, [years in 
        measuring period, days in months in 'sampling_months']        
    comm_castnet : numpy.ndarray
        Existing CASTNet O3 observations at each station contained in variable 
        'castnet_sites_fr' for years in measuring period, units of ppbv, 
        [years in measuring period, days in months in 'sampling_months']
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
    from scipy import stats
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
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.scatter(comm_t2m, comm_castnet)
    ax.set_xlabel('T$_\mathregular{2m}$ (K)')
    ax.set_ylabel('CASTNet (ppbv)')
    # find slope of linear fit
    slope, intercept, r_value, p_value, std_err = stats.linregress(
            np.hstack(comm_t2m), np.hstack(comm_castnet))
    # calculate the 95% confidence interval for the slope in a linear 
    # regression model using t-statistic found by looking up t-tables for 
    # (n - 2) degrees of freedom
    t = stats.t.ppf(1 - 0.025, np.hstack(comm_t2m).shape[0] - 2)
    confidence_interval = t * std_err
    # add slope to plot
    fig.text(0.15, 0.8, '%.2f $\pm$ %.3f ppbv K$^{\mathregular{-1}}$' 
             %(slope, confidence_interval), fontsize = 14)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_t2m_castneto3_%d-%d_%s.eps'%(years[0], years[-1],
                                                      region), dpi = 300)   
    return    
# # # # # # # # # # # # #
def boxplot_tracegas(castnet, mr2_gmi, ccmi_gmi, ffigac2_gmi, ffigac2hr_gmi, 
                     years, region, ylabel, obsname, filename): 
    """function plots boxplots of O3 distribution for CASTNet 
    observations and the four GMI CTM model cases/configurations. 
    Concentrations are regionally-averaged over the region used when 
    finding commensurate observed/modeled O3 concentrations. In boxplots, 
    line in box represents the median, whiskers are at 1.5 IQR (i.e. 
    [(Q1-1.5 IQR), (Q3+1.5 IQR)]) and flier points represent data that extend 
    beyond the whiskers (fliers).
    
    Parameters
    ----------    
    castnet : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period * days in months in 'sampling_months',]        
    mr2_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ccmi_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-CCMI, units 
        of ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2hr_gmi : numpy.ndarray 
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2-HighRes, 
        units of ppbv, [years in measuring period * days in months in 
        'sampling_months',]       
    years : list
        Years of interest    
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
    ylabel : str
        Axis label for y-axis (i.e. name of trace gas and its units)
    obsname : str
        Observational network used (i.e. AQS or CASTNet )
    file : str
        Prefix of the filename of the output figure    
    
    Returns
    ----------
    None
    """ 
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.font_manager as font_manager
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    # set custom s
    prop = font_manager.FontProperties(fname = pollutants_constants.FONTPATH_LIGHT)   
    # since sites have missing data (NaN), filter data using np.isnan
    # or ValueError will be raised when making boxplots
    # from https://stackoverflow.com/questions/44305873/how-to-deal-with-nan-value-when-plot-boxplot-using-python
    # for CASTNet
    mask = ~np.isnan(np.hstack(castnet))
    castnet = [d[m] for d, m in zip(np.hstack(castnet).T, mask.T)]
    # for HindcastMR2
    mask = ~np.isnan(np.hstack(mr2_gmi))
    mr2_gmi = [d[m] for d, m in zip(np.hstack(mr2_gmi).T, mask.T)]
    # for HindcastMR2-CCMI
    mask = ~np.isnan(np.hstack(ccmi_gmi))
    ccmi_gmi = [d[m] for d, m in zip(np.hstack(ccmi_gmi).T, mask.T)]
    # for HindcastFFIgac2
    mask = ~np.isnan(np.hstack(ffigac2_gmi))
    ffigac2_gmi = [d[m] for d, m in zip(np.hstack(ffigac2_gmi).T, mask.T)]
    # for HindcastFFIgac2-HighRes
    mask = ~np.isnan(np.hstack(ffigac2hr_gmi))
    ffigac2hr_gmi = [d[m] for d, m in zip(np.hstack(ffigac2hr_gmi).T, mask.T)]
    # dictionary of data for boxplots
    data = {}
    data['a'] = np.hstack(castnet)
    data['b'] = np.hstack(mr2_gmi)
    data['c'] = np.hstack(ffigac2hr_gmi)
    data['d'] = np.hstack(ccmi_gmi)
    data['e'] = np.hstack(ffigac2_gmi)
    color_dict = {'CASTNet' : 'w',
                  'AQS' : 'w',
                  'MR2' : '#33a02c', 
                  'FFIgac2-\nHighRes' : '#1f78b4',
                  'MR2-CCMI' : '#33a02c',              
                  'FFIgac2' : '#1f78b4'}
    controls = ['%s' %obsname, 'MR2', 'FFIgac2-\nHighRes', 'MR2-CCMI', 
                'FFIgac2']
    # initialize figure, axis
    ax = plt.subplot2grid((1, 1), (0, 0))
    boxplot_dict = ax.boxplot(
            [data[x] for x in ['a', 'b', 'c', 'd', 'e']],
            positions = [1, 1.5, 2, 2.5, 3],
            labels = controls, 
            patch_artist = True,
            widths = 0.25, showfliers = True)
    # aesthetics
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
        #b.set_markevery(every = 2)
        b.set_markersize(1)
        b.set_markerfacecolor('k')
    for b in boxplot_dict['whiskers']:
        b.set_linewidth(2.)
    for b in boxplot_dict['caps']:
        b.set_linewidth(2.)
    # axis colors corresponding to high vs. coarse resolution
    # https://stackoverflow.com/questions/23248435/
    # fill-between-two-vertical-lines-in-matplotlib    
    ax.axvspan(1.25, 2.25, alpha = 0.7, color = 'gainsboro', 
               label = '1.0$^{\circ}$ x 1.25$^{\circ}$')    
    ax.axvspan(2.25, 3.25, alpha = 0.3, color = 'lightcoral', 
               label = '2.0$^{\circ}$ x 2.5$^{\circ}$')       
    ax.set_xlim([0.75, 3.25])
    # add custom legend, patches corresponding to colors on plot
    patch_lr = mpatches.Patch(color = 'gainsboro', alpha = 0.7, 
                             label = '1.25$^{\circ}$ x 1.0$^{\circ}$')    
    patch_hr = mpatches.Patch(color = 'lightcoral', alpha = 0.3, 
                             label = '2.5$^{\circ}$ x 2.0$^{\circ}$')       
    patch_m1 = mpatches.Patch(color = '#1f78b4', label = 'MERRA-1')
    patch_m2 = mpatches.Patch(color = '#33a02c', label = 'MERRA-2')
    leg = ax.legend(handles = [patch_lr, patch_hr, patch_m1, patch_m2], 
                    bbox_to_anchor = (1.02, 0.65))
    plt.setp(leg.get_texts(), fontproperties = prop)
    leg.get_frame().set_linewidth(0.0)
    #ax.set_ylim([0, 160])
    ax.set_ylabel('%s' %ylabel, fontproperties = prop)
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop)
    plt.subplots_adjust(right = 0.75)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'boxplot_%s_%d-%d_%s.eps' %(filename, years[0], years[-1], 
                                            region), dpi = 300)
    return 
# # # # # # # # # # # # #    
def cdf_castneto3gmio3(castnet, mr2_gmi, ccmi_gmi, ffigac2_gmi, ffigac2hr_gmi, 
    years, region):
    """function plots cumulative distribution functions ("S-curves") of O3 
    distribution for CASTNet observations and the four GMI CTM model cases/
    configurations. Concentrations are regionally-averaged over the region used 
    when finding commensurate observed/modeled O3 concentrations. CDFs are 
    places two subplots: one for configurations using MERRA-1 and one for 
    configurations using MERRA-2.
    
    Parameters
    ----------    
    castnet : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period * days in months in 'sampling_months',]        
    mr2_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ccmi_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-CCMI, units 
        of ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2hr_gmi : numpy.ndarray 
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2-HighRes, 
        units of ppbv, [years in measuring period * days in months in 
        'sampling_months',]       
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
    # from https://stackoverflow.com/questions/10640759/how-to-get-the-
    # cumulative-distribution-function-with-numpy
    N = len(np.hstack(castnet))
    F = np.array(range(N))/float(N)
    castnet = np.sort(np.hstack(castnet))
    mr2_gmi = np.sort(np.hstack(mr2_gmi))
    ccmi_gmi = np.sort(np.hstack(ccmi_gmi))
    ffigac2_gmi = np.sort(np.hstack(ffigac2_gmi))
    ffigac2hr_gmi = np.sort(np.hstack(ffigac2hr_gmi))
    # initialize figure, axes (one axis for MERRA-1, one for MERRA-2)
    fig = plt.figure()
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    # plot obs/MERRA-1
    ax1.plot(castnet, F, lw = 2.0, color = 'k', label = 'CASTNet')
    ax1.plot(ffigac2_gmi, F, lw = 1.5, color = pollutants_constants.COLOR_FFIGAC2, 
             label = 'FFIgac2')
    ax1.plot(ffigac2hr_gmi, F, lw = 1.5, 
             color = pollutants_constants.COLOR_FFIGAC2HIGHRES, 
             label = 'FFIgac2-\nHighRes')
    ax1.legend(loc = 'lower right', frameon = False, fontsize = 8)
    # plot obs/MERRA-2
    ax2.plot(castnet, F, lw = 2.0, color = 'k', label = 'CASTNet')
    ax2.plot(mr2_gmi, F, lw = 1.5, 
             color = pollutants_constants.COLOR_MR2, label = 'MR2')
    ax2.plot(ccmi_gmi, F, lw = 1.5, color = pollutants_constants.COLOR_MR2CCMI, 
             label = 'MR2-CCMI')
    ax2.legend(loc = 'lower right', frameon = False, fontsize = 8)
    ax1.set_xlim([20, 80])
    ax2.set_xlim([20, 80])
    ax2.set_yticklabels([''])
    # add axis labels
    ax1.set_ylabel('Cumulative Probability', fontsize = 12)
    ax1.set_xlabel('Ozone [ppbv]', x = 1.1, fontsize = 12)
    ax1.set_title('MERRA-1', fontsize = 12)
    ax2.set_title('MERRA-2', fontsize = 12)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'cdf_castneto3gmio3_%d-%d_%s.eps' %(years[0], years[-1],
                                                    region), dpi = 300)
    return 
# # # # # # # # # # # # #
def pdf_castneto3gmio3(castnet, mr2_gmi, ccmi_gmi, ffigac2_gmi, 
                       ffigac2hr_gmi, years, region):
    """function plots probability density functions of O3 distribution for 
    CASTNet observations and the four GMI CTM model cases/configurations.
    Concentrations are regionally-averaged over the region used when 
    finding commensurate observed/modeled O3 concentrations. 
    
    Parameters
    ----------    
    castnet : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period * days in months in 'sampling_months',]        
    mr2_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ccmi_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-CCMI, units 
        of ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2hr_gmi : numpy.ndarray 
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2-HighRes, 
        units of ppbv, [years in measuring period * days in months in 
        'sampling_months',]       
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
    import matplotlib.mlab as mlab
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
    fig = plt.figure()
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    ax1.set_title('MERRA-1', fontsize = 12)
    ax2.set_title('MERRA-2', fontsize = 12)
    # obs
    n, bins, patches = ax1.hist(np.hstack(castnet), 25, normed = 1, alpha = 0.)
    pdf_castnet = mlab.normpdf(bins, np.mean(castnet), np.std(castnet))
    ax1.plot(bins, pdf_castnet, lw = 2., color = 'k', label = 'CASTNet')
    ax2.plot(bins, pdf_castnet, lw = 2., color = 'k', label = 'CASTNet')
    # HindcastFFIgac2 
    n, bins, patches = ax1.hist(np.hstack(ffigac2_gmi), 25, normed = 1, alpha = 0.)
    pdf_ffigac2 = mlab.normpdf(bins, np.mean(ffigac2_gmi), np.std(ffigac2_gmi))
    ax1.plot(bins, pdf_ffigac2, lw = 1.5, 
             color = pollutants_constants.COLOR_FFIGAC2, label = 'FFIgac2')
    # HindcastFFIgac2-HighRes
    n, bins, patches = ax1.hist(np.hstack(ffigac2hr_gmi), 25, normed = 1, alpha = 0.)
    pdf_ffigac2hr = mlab.normpdf(bins, np.mean(ffigac2hr_gmi), np.std(ffigac2hr_gmi))
    ax1.plot(bins, pdf_ffigac2hr, lw = 1.5, 
             color = pollutants_constants.COLOR_FFIGAC2HIGHRES, 
             label = 'FFIgac2-HighRes')
    # HindcastMR2
    n, bins, patches = ax2.hist(np.hstack(mr2_gmi), 25, normed = 1, alpha = 0.)
    pdf_mr2 = mlab.normpdf(bins, np.mean(mr2_gmi), np.std(mr2_gmi))
    ax2.plot(bins, pdf_mr2, lw = 1.5, 
             color = pollutants_constants.COLOR_MR2, label = 'MR2')   
    ax1.legend(loc = 'upper left', frameon = False, fontsize = 10)
    # HindcastMR2-CCMI
    n, bins, patches = ax2.hist(np.hstack(ccmi_gmi), 25, normed = 1, alpha = 0.)
    pdf_ccmi = mlab.normpdf(bins, np.mean(ccmi_gmi), np.std(ccmi_gmi))
    ax2.plot(bins, pdf_ccmi, lw = 1.5, 
             color = pollutants_constants.COLOR_MR2CCMI, label = 'MR2-CCMI')    
    ax2.legend(loc = 'upper left', frameon = False, fontsize = 10)
    ax1.set_xlim([20, 80])
    ax2.set_xlim([20, 80])
    ax1.set_ylim([0.0, 0.09])
    ax2.set_ylim([0.0, 0.09])    
    ax2.set_yticklabels([''])
    # add axis labels
    ax1.set_ylabel('Probability', fontsize = 12)
    ax1.set_xlabel('Ozone [ppbv]', x = 1.1, fontsize = 12)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'pdf_castneto3gmio3_%d-%d_%s.eps' %(years[0], years[-1],
                                                    region), dpi = 300)
    return
# # # # # # # # # # # # #  
def scatter_castneto3_gmio3(castnet, mr2_gmi, ccmi_gmi, ffigac2_gmi, 
                            ffigac2hr_gmi, years, region): 
    """function plots scatterplots of regionally-averaged GMI O3 vs. CASTNet
    O3 for each of the four GMI CTM model cases/configurations for each year
    in the measuring period. Concentrations are regionally-averaged over the 
    region used when finding commensurate observed/modeled O3 concentrations. 
    A 1:1 line is also plotted, and the slope of the linear regression for each 
    year is found and included on each plot.
    
    Parameters
    ----------    
    castnet : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period * days in months in 'sampling_months',]        
    mr2_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ccmi_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-CCMI, units 
        of ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2hr_gmi : numpy.ndarray 
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2-HighRes, 
        units of ppbv, [years in measuring period * days in months in 
        'sampling_months',]       
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
    from scipy import stats
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
    # make colors same dimension as years (with current configuration, 
    # function wouldn't break until > 6 years are passed in)
    colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']
    colors = colors[0 : len(years)]
    # initialize figure, axes
    fig = plt.figure()
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax3 = plt.subplot2grid((2, 2), (1, 0))
    ax4 = plt.subplot2grid((2, 2), (1, 1))    
    axs = [ax1, ax2, ax3, ax4]    
    labels = ['MR2', 'FFIgac2-HighRes', 'MR2-CCMI', 'FFIgac2']
    i = 0
    # look through model configurations 
    for model in [mr2_gmi, ffigac2hr_gmi, ccmi_gmi, ffigac2_gmi]: 
        j = 0 
        k = 0.88
        for year in years: 
            axs[i].scatter(castnet[j], model[j], color = colors[j], 
                           marker = 'o', s = 15, label = year)
            # find slope of linear fit 
            b = stats.linregress(castnet[j], model[j]).slope    
            axs[i].text(0.05, k, '%.2f' %b, color = colors[j], 
                        transform = axs[i].transAxes)
            k  = k - 0.1
            j = j + 1
        # plot 1:1 line
        axs[i].plot(np.arange(20, 80, 0.25), np.arange(20, 80, 0.25), 
                    '--k', linewidth = 0.5)
        axs[i].set_ylabel('%s' %(labels[i]), fontsize = 12)
        axs[i].set_xlim([25, 75])
        axs[i].set_ylim([25, 75])
        i = i + 1
    ax1.set_xticklabels([''])
    ax2.set_xticklabels([''])
    ax3.set_xlabel('CASTNet', fontsize = 12)
    ax4.set_xlabel('CASTNet', fontsize = 12)
    ax2.set_yticklabels([''])
    ax4.set_yticklabels([''])
    plt.tight_layout()
    plt.subplots_adjust(bottom = 0.2)
    # add legend
    ax3.legend(bbox_to_anchor = (1.9, -0.25), ncol = 3, fontsize = 12, 
               scatterpoints = 3, frameon = False)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_castneto3_gmio3_%d-%d_%s.eps'%(years[0], years[-1],
                                                        region), dpi = 300)   
    return
# # # # # # # # # # # # #  
def timeseries_castneto3gmio3(castnet, mr2_gmi, ccmi_gmi, ffigac2_gmi, 
                              ffigac2hr_gmi, years, region): 
    """function plots timeseries of regionally-averaged GMI O3 vs. CASTNet
    O3 for each of the four GMI CTM model cases/configurations for all years
    in the measuring period. Concentrations are regionally-averaged over the 
    region used when finding commensurate observed/modeled O3 concentrations. 
    The Pearson correlation coefficient is also calculated.
    
    Parameters
    ----------    
    castnet : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period * days in months in 'sampling_months',]        
    mr2_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ccmi_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-CCMI, units 
        of ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2_gmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2hr_gmi : numpy.ndarray 
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2-HighRes, 
        units of ppbv, [years in measuring period * days in months in 
        'sampling_months',]       
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
    from scipy import stats
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
    fig = plt.figure(figsize = (9, 6.5))
    ax1 = plt.subplot2grid((4, 4), (0, 0), colspan = 4)
    ax2 = plt.subplot2grid((4, 4), (1, 0), colspan = 4)
    ax3 = plt.subplot2grid((4, 4), (2, 0), colspan = 4)
    ax4 = plt.subplot2grid((4, 4), (3, 0), colspan = 4)
    axs = [ax1, ax2, ax3, ax4]
    # CASTNet/HindcastMR2
    ax1.plot(np.hstack(castnet), lw = 2., color = 'k')
    ax1.plot(np.hstack(mr2_gmi), lw = 1.5, color = pollutants_constants.COLOR_MR2)
    ax1.set_xlim([0, len(np.hstack(castnet))])
    ax1.set_xticks([x * castnet.shape[1] for x in np.arange(0, castnet.shape[0])])
    ax1.set_xticklabels([''])
    r = stats.linregress(np.hstack(castnet), np.hstack(mr2_gmi)).rvalue
    ax1.set_title('MR2, r = %.2f' %r, fontsize = 12, 
                  color = pollutants_constants.COLOR_MR2)
    # CASTNet/HindcastFFIgac2-HighRes
    ax2.plot(np.hstack(castnet), lw = 2., color = 'k')
    ax2.plot(np.hstack(ffigac2hr_gmi), lw = 1.5, 
             color = pollutants_constants.COLOR_FFIGAC2HIGHRES)
    ax2.set_xlim([0, len(np.hstack(castnet))])
    ax2.set_xticks([x * castnet.shape[1] for x in np.arange(0, castnet.shape[0])])
    ax2.set_xticklabels([''])
    r = stats.linregress(np.hstack(castnet), np.hstack(ffigac2hr_gmi)).rvalue
    ax2.set_title('FFIgac2-HighRes, r = %.2f' %r, fontsize = 12,
                  color = pollutants_constants.COLOR_FFIGAC2HIGHRES)
    # CASTNet/HindcastMR2-CCMI
    ax3.plot(np.hstack(castnet), lw = 2., color = 'k')
    ax3.plot(np.hstack(ccmi_gmi), lw = 1.5, 
             color = pollutants_constants.COLOR_MR2CCMI)
    ax3.set_xlim([0, len(np.hstack(castnet))])
    ax3.set_xticks([x * castnet.shape[1] for x in np.arange(0, castnet.shape[0])])
    ax3.set_xticklabels([''])
    r = stats.linregress(np.hstack(castnet), np.hstack(ccmi_gmi)).rvalue
    ax3.set_title('MR2-CCMI, r = %.2f' %r, fontsize = 12, 
                  color = pollutants_constants.COLOR_MR2CCMI)
    # CASTNet/HindcastFFIgac2
    ax4.plot(np.hstack(castnet), lw = 2., color = 'k')
    ax4.plot(np.hstack(ffigac2_gmi), lw = 1.5, 
             color = pollutants_constants.COLOR_FFIGAC2)
    ax4.set_xlim([0, len(np.hstack(castnet))])
    ax4.set_xticks([x * castnet.shape[1] for x in np.arange(0, castnet.shape[0])])
    ax4.set_xticklabels(years)      
    r = stats.linregress(np.hstack(castnet), np.hstack(ffigac2_gmi)).rvalue
    ax4.set_title('FFIgac2, r = %.2f' %r, fontsize = 12, 
                  color = pollutants_constants.COLOR_FFIGAC2)
    plt.tight_layout()
    for ax in axs:
        ax.set_ylim([25, 75])
    ax3.set_ylabel('Ozone [ppbv]', y = 1.2, fontsize = 12)
    plt.subplots_adjust(left = 0.1)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'timeseries_castneto3gmio3_%d-%d_%s.eps' %(years[0], years[-1],
                                                           region), dpi = 300)
    return
# # # # # # # # # # # # #
def map_gmiaqscastnet(castnet, castnet_sites_fr, aqs_co_coords, aqs_no2_coords, 
                      years, sampling_months, sampling_hours, region, timezone):        
    """for the CASTNet sites specified in variable 'castnet_sites_fr,' function
    opens the locations (latitude/longitude coordinates) of those CASTNet 
    sites and their co-located (or nearly co-located) MERRA grid cells, GMI 
    sites, and AQS CO and NO2 sites and plots their locations. 
    
    Parameters
    ----------    
    castnet : numpy.ndarray
        Existing CASTNet O3 observations at each station contained in variable 
        'castnet_sites_fr' for years in measuring period, units of ppbv, 
        [years in measuring period, stations in 'castnet_sites_fr', 
        days in months in 'sampling_months']
    castnet_sites_fr : list
        CASTNET site names in focus region
    aqs_co_coords : list
        Coordinates of unique AQS stations measuring CO within bounding boxes 
        defined by CASTNet stations
    aqs_no2_coords : list
        Coordinates of unique AQS stations measuring NO2 within bounding boxes 
        defined by CASTNet stations
    years : list
        Years of interest
    sampling_months : list 
        Months of interest
    sampling_hours : list 
        Hours during which trace gas concentrations are returned
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
    
    Returns
    ----------
    None    
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from matplotlib.patches import Polygon 
    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
    from matplotlib.path import Path
    import matplotlib.patches as patches
    import shapely.geometry as sg
    import shapely.ops as so
    from descartes import PolygonPatch    
    import matplotlib as mpl    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/GMI/')
    import commensurability
    # change hatch width
    mpl.rcParams['hatch.linewidth'] = 2.0
    # open CASTNet and GMI sites
    castnet_lats, castnet_lons, gmi_lats, gmi_lons = \
    commensurability.commensurate_castnet_gmi_locations(castnet_sites_fr, 
        sampling_months, sampling_hours, timezone)
    # open MERRA-2 (requires a year's worth of GMI CTM output)
    mr2_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
    commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2', 
        years, sampling_months, sampling_hours, timezone)        
    (comm_t2m, comm_t10m, comm_u2m, comm_u10m, comm_v2m, comm_v10m, 
     merra_lats_fr, merra_lons_fr) = commensurability.commensurate_MERRA2(
        mr2_castnet, castnet_sites_fr, years, sampling_months, sampling_hours)        
    # initialize figure, axes
    fig = plt.figure(figsize = (15,7))
    ax = plt.subplot2grid((1, 1), (0, 0))
    # focus region map 
    llcrnrlon = -84.
    llcrnrlat = 36.
    urcrnrlon = -66.3
    urcrnrlat = 48.
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
                llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, 
                urcrnrlat = urcrnrlat, resolution = 'h', area_thresh = 1000)
    m.drawmapboundary(color = '#888888', fill_color = '#dcf0fa')
    m.fillcontinents(color = '#f9f6d8', lake_color = '#dcf0fa')
    m.drawlsmask(ocean_color = '#dcf0fa')
    m.drawcoastlines(color = '#888888')
    # plot GMI sites 
    x_gmi, y_gmi = m(gmi_lons, gmi_lats)
    for xi, yi in zip(x_gmi, y_gmi):
        # find longitude, latitude of every unique grid cell
        xinverse, yinverse = m(xi, yi, inverse = True)
        # change bounds of grid cell to map projection (MERRA-2 resolution is 
        # 0.625˚ x 0.5˚)
        x1, y1 = m(xinverse + 0.625, yinverse - 0.5)
        x2, y2 = m(xinverse + 0.625, yinverse + 0.5)
        x3, y3 = m(xinverse - 0.625, yinverse + 0.5)
        x4, y4 = m(xinverse - 0.625, yinverse - 0.5)
        # add patch to map
        p = Polygon([(x1, y1),(x2, y2),(x3,y3),(x4, y4)], facecolor = 
                     pollutants_constants.COLOR_CHEMISTRY, alpha = 1., 
                     edgecolor = 'none', linewidth = 0.0, zorder = 5)
        plt.gca().add_patch(p) 
#    # plot MERRA-2 reanalysis
#    x_merra, y_merra = m(merra_lons_fr, merra_lats_fr)
#    for xi, yi in zip(x_merra, y_merra):
#        # find longitude, latitude of every unique grid cell
#        xinverse, yinverse = m(xi, yi, inverse = True)
#        # change bounds of grid cell to map projection (MERRA-2 resolution is 
#        # 0.625˚ x 0.5˚)
#        x1, y1 = m(xinverse + 0.3125, yinverse - 0.25)
#        x2, y2 = m(xinverse + 0.3125, yinverse + 0.25)
#        x3, y3 = m(xinverse - 0.3125, yinverse + 0.25)
#        x4, y4 = m(xinverse - 0.3125, yinverse - 0.25)
#        # add patch to map
#        p = Polygon([(x1, y1),(x2, y2),(x3,y3),(x4, y4)], fill = False, 
#                     color = '#ff7f00', hatch = '///', linewidth = 0.0, 
#                     zorder = 10)
#        plt.gca().add_patch(p) 
    # plot CASTNet stations
    x_castnet, y_castnet = m(castnet_lons, castnet_lats)
    castnet_loc = m.scatter(x_castnet, y_castnet, 50, color = 'k', 
                            marker = 's', edgecolor = 'k', 
                            linewidth = 0.5, zorder = 15, label = 'CASTNet')
    ## plot AQS stations measuring NO2
    #aqs_no2 = np.array([np.vstack(aqs_no2_coords)[:, 1], 
    #                    np.vstack(aqs_no2_coords)[:, 0]], dtype = float).T
    #aqs_no2 = pd.DataFrame(aqs_no2)
    #aqs_no2 = aqs_no2.drop_duplicates()
    #x_aqs_no2, y_aqs_no2 = m(aqs_no2[0].values, 
    #                         aqs_no2[1].values)
    #aqs_no2_loc = m.scatter(x_aqs_no2, y_aqs_no2, 8, color = '#33a02c', 
    #                        marker = 'o', edgecolor = 'none', linewidth = 0.5, 
    #                        zorder = 21, label = 'AQS NO$_{2}$')
    ## plot AQS stations measuring CO
    #aqs_co = np.array([np.vstack(aqs_co_coords)[:, 1], 
    #                   np.vstack(aqs_co_coords)[:, 0]], dtype = float).T
    #aqs_co = pd.DataFrame(aqs_co)
    #aqs_co = aqs_co.drop_duplicates()
    #x_aqs_co, y_aqs_co = m(aqs_co[0].values, aqs_co[1].values)
    #aqs_co_loc = m.scatter(x_aqs_co, y_aqs_co, 8, color = '#b2df8a', marker = 'd', 
    #                       edgecolor = 'none', linewidth = 0.5, zorder = 22, 
    #                       label = 'AQS CO')
    # outline focus region
    m.readshapefile(pollutants_constants.PATH_SHAPEFILES + 
                    'cb_2015_us_state_20m', name = 'states', 
                    drawbounds = True, color = '#888888')
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
        'New Mexico' : 35, 'Nevada': 32, 'New York' : 36, 'Ohio' : 39, 
        'Oklahoma' : 40, 'Oregon' : 41, 'Pennsylvania' : 42, 
        'Rhode Island' : 44, 'South Carolina' : 45, 'South Dakota' : 46, 
        'Tennessee' : 47, 'Texas' : 48, 'Utah' : 49, 'Virginia' : 51, 
        'Vermont' : 50, 'Washington' : 53, 'Wisconsin' : 55, 
        'West Virginia' : 54, 'Wyoming' : 56}
    # iterate through states, if states are in the Northeastern United States, 
    # append shapely.geometry.polygon.Polygon obejcts to list
    patches_union = [] 
    for key, value in state_fips_code_listing.items():
        if key in pollutants_constants.NORTHEAST_STATES: 
                for info, shape in zip(m.states_info, m.states):
                    if info['NAME'] == key:
                        patches_union.append(sg.Polygon(shape))
    # cascaded union can work on a list of shapes, adapted from 
    # https://stackoverflow.com/questions/34475431/plot-unions-of-polygons-in-matplotlib
    neus = so.cascaded_union(patches_union) 
    ax.add_patch(PolygonPatch(neus, fill = False, ec = '#888888', 
                              zorder = 2, linewidth = 4.0)) 
    # generate inset axes
    ax_small = zoomed_inset_axes(ax, 0.055, loc = 'lower right')
    m_small = Basemap(projection = 'ortho', lat_0 = 30, lon_0 = -80,
                      resolution = 'c', area_thresh = 10000) 
    m_small.drawcoastlines(color = '#888888')
    m_small.drawcountries(color = '#888888')
    m_small.fillcontinents(color = '#f9f6d8')
    m_small.drawmapboundary(color = '#888888')
    # draw the zoom rectangle around focus region
    lbx, lby = m_small(*m(m.xmin, m.ymin, inverse= True))
    ltx, lty = m_small(*m(m.xmin, m.ymax, inverse= True))
    rtx, rty = m_small(*m(m.xmax, m.ymax, inverse= True))
    rbx, rby = m_small(*m(m.xmax, m.ymin, inverse= True))
    verts = [(lbx, lby), # left, bottom
        (ltx, lty), # left, top
        (rtx, rty), # right, top
        (rbx, rby), # right, bottom
        (lbx, lby), # ignored
        ]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor = 'r', edgecolor = 'r', 
                              fill = False, lw = 2)
    ax_small.text(0.72, 0.53, 'Map\nArea', horizontalalignment = 'center',
                  fontsize = 14, verticalalignment = 'center', 
                  transform = ax_small.transAxes)
    ax_small.add_patch(patch)
    # create legend
    merra_patch = patches.Patch(fill = False, color = '#ff7f00', 
                                 hatch = '///', linewidth = 0.0, 
                                 label = 'MERRA-2')
    gmi_patch = patches.Patch(facecolor = pollutants_constants.COLOR_CHEMISTRY, 
                               alpha = 1., edgecolor = 'none', linewidth = 0.0,
                               label = 'GMI')
    leg = ax.legend(handles = [castnet_loc, gmi_patch], 
                    loc = 'upper left', ncol = 2, 
                    fontsize = 16, scatterpoints = 1, framealpha = 0.0)
    leg.get_frame().set_linewidth(0.0)    
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_gmiaqscastnet_%s_UPDATED.png' %(region), dpi = 500)
    return 
# # # # # # # # # # # # #    
def diurnal_castneto3gmio3(castnet_o3_d, mr2_o3_d, ccmi_o3_d, ffigac2_o3_d,
                           ffigac2hr_o3_d, years, region): 
    """function plots the diurnal cycle of O3 averaged over all sites in 
    region and over days and years in measuring period. Diurnal cycles 
    represent O3 observations from CASTNet (overnight hours 0000-0300 missing?)
    and modelled O3 concentrations from 4 configurations of GMI CTM. 

    Parameters
    ----------    
    castnet_o3_d : numpy.ndarray
        Regionally- and time-averaged CASTNet O3 diurnal cycle observations, 
        units of ppbv, [24,]    
    mr2_o3_d : numpy.ndarray
        Regionally- and time-averaged O3 diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2, units of ppbv, [24,]    
    ccmi_o3_d : numpy.ndarray
        Regionally- and time-averaged O3 diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2-CCMI, units of ppbv, [24,]    
    ffigac2_o3_d : numpy.ndarray
        Regionally- and time-averaged O3 diurnal cycle concentrations from GMI 
        CTM model configuration HindcastFFIgac2, units of ppbv, [24,]    
    ffigac2hr_o3_d : numpy.ndarray
        Regionally- and time-averaged O3 diurnal cycle concentrations from GMI 
        CTM model configuration HindcastFFIgac2-HighRes, units of ppbv, [24,]    
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
    fig = plt.figure(figsize = (9, 4))
    ax = plt.subplot2grid((1, 2), (0, 0), colspan = 2)
    ax.plot(castnet_o3_d, lw = 2.0, color = 'k', label = 'CASTNet')
    ax.plot(ffigac2_o3_d, lw = 1.5, color = pollutants_constants.COLOR_FFIGAC2, 
            label = 'FFIgac2')
    ax.plot(ffigac2hr_o3_d, lw = 1.5, 
            color = pollutants_constants.COLOR_FFIGAC2HIGHRES, 
            label = 'FFIgac2-\nHighRes')
    ax.plot(mr2_o3_d, lw = 1.5, color = pollutants_constants.COLOR_MR2, 
            label = 'MR2')
    ax.plot(ccmi_o3_d, lw = 1.5, color = pollutants_constants.COLOR_MR2CCMI, 
            label = 'MR2-CCMI')
    ax.legend(loc = 'lower right', frameon = False, fontsize = 10)
    ax.set_xlim([0, 23])
    ax.set_xlabel('Hour [GMT]')
    ax.set_ylabel('O$_{3}$ [ppbv]')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'diurnal_castneto3gmio3_%d-%d_%s.eps' %(years[0], years[-1],
                                                        region), dpi = 300)
    return
# # # # # # # # # # # # #
def diurnal_aqscogmico(aqs_co_d, mr2_co_d, ccmi_co_d, ffigac2_co_d, 
                       ffigac2hr_co_d, years, region): 
    """function plots the diurnal cycle of CO averaged over all sites in 
    region and over days and years in measuring period. Diurnal cycles 
    represent CO observations from AQS which are within 1 deg of CASTNet sites    
    and modelled CO concentrations from 4 configurations of GMI CTM. 

    Parameters
    ----------    
    aqs_co_d : numpy.ndarray
        Regionally- and time-averaged AQS CO diurnal cycle observations, 
        units of ppbv, [24,]    
    mr2_co_d : numpy.ndarray
        Regionally- and time-averaged CO diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2, units of ppmv, [24,]    
    ccmi_co_d : numpy.ndarray
        Regionally- and time-averaged CO diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2-CCMI, units of ppmv, [24,]    
    ffigac2_co_d : numpy.ndarray
        Regionally- and time-averaged CO diurnal cycle concentrations from GMI 
        CTM model configuration HindcastFFIgac2, units of ppmv, [24,]    
    ffigac2hr_co_d : numpy.ndarray
        Regionally- and time-averaged CO diurnal cycle concentrations from GMI 
        CTM model configuration HindcastFFIgac2-HighRes, units of ppmv, [24,]    
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
    fig = plt.figure(figsize = (9, 4))
    ax = plt.subplot2grid((1, 2), (0, 0), colspan = 2)
    ax.plot(aqs_co_d, lw = 2.0, color = 'k', label = 'AQS')
    ax.plot(ffigac2_co_d, lw = 1.5, color = pollutants_constants.COLOR_FFIGAC2, 
            label = 'FFIgac2')
    ax.plot(ffigac2hr_co_d, lw = 1.5, 
            color = pollutants_constants.COLOR_FFIGAC2HIGHRES, 
            label = 'FFIgac2-\nHighRes')
    ax.plot(mr2_co_d, lw = 1.5, color = pollutants_constants.COLOR_MR2, 
            label = 'MR2')
    ax.plot(ccmi_co_d, lw = 1.5, color = pollutants_constants.COLOR_MR2CCMI, 
            label = 'MR2-CCMI')
    ax.legend(loc = 'upper right', frameon = False, fontsize = 10, 
              ncol = 2)
    ax.set_xlim([0, 23])
    ax.set_xlabel('Hour [GMT]')
    ax.set_ylabel('CO [ppmv]')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'diurnal_aqscogmico_%d-%d_%s.eps' %(years[0], years[-1],
                                                    region), dpi = 300)    
    return
# # # # # # # # # # # # # 
def diurnal_aqsno2gmino2(aqs_no2_d, mr2_no2_d, ccmi_no2_d, ffigac2_no2_d, 
                         ffigac2hr_no2_d, years, region): 
    """function plots the diurnal cycle of NO2 averaged over all sites in 
    region and over days and years in measuring period. Diurnal cycles 
    represent NO2 observations from AQS which are within 1 deg of CASTNet sites    
    and modelled NO2 concentrations from 4 configurations of GMI CTM. 

    Parameters
    ----------    
    aqs_no2_d : numpy.ndarray
        Regionally- and time-averaged AQS CO diurnal cycle observations, 
        units of ppbv, [24,]    
    mr2_no2_d : numpy.ndarray
        Regionally- and time-averaged CO diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2, units of ppmv, [24,]    
    ccmi_no2_d : numpy.ndarray
        Regionally- and time-averaged CO diurnal cycle concentrations from GMI 
        CTM model configuration HindcastMR2-CCMI, units of ppmv, [24,]    
    ffigac2_no2_d : numpy.ndarray
        Regionally- and time-averaged CO diurnal cycle concentrations from GMI 
        CTM model configuration HindcastFFIgac2, units of ppmv, [24,]    
    ffigac2hr_no2_d : numpy.ndarray
        Regionally- and time-averaged CO diurnal cycle concentrations from GMI 
        CTM model configuration HindcastFFIgac2-HighRes, units of ppmv, [24,]    
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
    fig = plt.figure(figsize = (9, 4))
    ax = plt.subplot2grid((1, 2), (0, 0), colspan = 2)
    ax.plot(aqs_no2_d, lw = 2.0, color = 'k', label = 'AQS')
    ax.plot(ffigac2_no2_d, lw = 1.5, color = pollutants_constants.COLOR_FFIGAC2, 
            label = 'FFIgac2')
    ax.plot(ffigac2hr_no2_d, lw = 1.5, 
            color = pollutants_constants.COLOR_FFIGAC2HIGHRES, 
            label = 'FFIgac2-\nHighRes')
    ax.plot(mr2_no2_d, lw = 1.5, color = pollutants_constants.COLOR_MR2, 
            label = 'MR2')
    ax.plot(ccmi_no2_d, lw = 1.5, color = pollutants_constants.COLOR_MR2CCMI, 
            label = 'MR2-CCMI')
    ax.legend(loc = 'upper right', frameon = False, fontsize = 10, 
              ncol = 2)
    ax.set_xlim([0, 23])
    ax.set_xlabel('Hour [GMT]')
    ax.set_ylabel('NO$_{2}$ [ppbv]')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'diurnal_aqsno2gmino2_%d-%d_%s.eps' %(years[0], years[-1],
                                                      region), dpi = 300) 
    return
# # # # # # # # # # # # #    
def pdf_aqscogmico(aqs_co, mr2_gmi, ccmi_gmi, ffigac2_gmi, ffigac2hr_gmi, 
                   years, region):
    """function plots probability density functions of CO distribution for 
    AQS observations and the four GMI CTM model cases/configurations.
    Concentrations are regionally-averaged over the region used when 
    finding commensurate observed/modeled CO concentrations. 
    
    Parameters
    ----------    
    aqs_co : numpy.ndarray
        AQS CO observations in region, units of ppbm, [years in measuring 
        period * days in months in 'sampling_months',]        
    mr2_gmi : numpy.ndarray
        GMI CTM CO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ccmi_gmi : numpy.ndarray
        GMI CTM CO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-CCMI, units 
        of ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2_gmi : numpy.ndarray
        GMI CTM CO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2hr_gmi : numpy.ndarray 
        GMI CTM CO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2-HighRes, 
        units of ppbv, [years in measuring period * days in months in 
        'sampling_months',]       
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
    import matplotlib.mlab as mlab
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
    fig = plt.figure()
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    ax1.set_title('MERRA-1', fontsize = 12)
    ax2.set_title('MERRA-2', fontsize = 12)
    # obs
    n, bins, patches = ax1.hist(np.hstack(aqs_co), 25, normed = 1, alpha = 0.)
    pdf_aqs = mlab.normpdf(bins, np.mean(aqs_co), np.std(aqs_co))
    ax1.plot(bins, pdf_aqs, lw = 2., color = 'k', label = 'AQS')
    ax2.plot(bins, pdf_aqs, lw = 2., color = 'k', label = 'AQS')
    ## HindcastFFIgac2
    n, bins, patches = ax1.hist(np.hstack(ffigac2_gmi), 25, normed = 1, alpha = 0.)
    pdf_ffigac2 = mlab.normpdf(bins, np.mean(ffigac2_gmi), np.std(ffigac2_gmi))
    ax1.plot(bins, pdf_ffigac2, lw = 1.5, 
             color = pollutants_constants.COLOR_FFIGAC2, label = 'FFIgac2')
    # HindcastFFIgac2-HighRes
    n, bins, patches = ax1.hist(np.hstack(ffigac2hr_gmi), 25, normed = 1, alpha = 0.)
    pdf_ffigac2hr = mlab.normpdf(bins, np.mean(ffigac2hr_gmi), np.std(ffigac2hr_gmi))
    ax1.plot(bins, pdf_ffigac2hr, lw = 1.5, 
             color = pollutants_constants.COLOR_FFIGAC2HIGHRES, 
             label = 'FFIgac2-HighRes')
    ax1.legend(loc = 'upper left', frameon = False, fontsize = 10)
    # HindcastMR2
    n, bins, patches = ax2.hist(np.hstack(mr2_gmi), 25, normed = 1, alpha = 0.)
    pdf_mr2 = mlab.normpdf(bins, np.mean(mr2_gmi), np.std(mr2_gmi))
    ax2.plot(bins, pdf_mr2, lw = 1.5, 
             color = pollutants_constants.COLOR_MR2, label = 'MR2')   
    # HindcastMR2-CCMI
    n, bins, patches = ax2.hist(np.hstack(ccmi_gmi), 25, normed = 1, alpha = 0.)
    pdf_ccmi = mlab.normpdf(bins, np.mean(ccmi_gmi), np.std(ccmi_gmi))
    ax2.plot(bins, pdf_ccmi, lw = 1.5, 
             color = pollutants_constants.COLOR_MR2CCMI, label = 'MR2-CCMI')    
    ax2.legend(loc = 'upper left', frameon = False, fontsize = 10)
    ax1.set_xlim([0.05, 0.40])
    ax2.set_xlim([0.05, 0.40])
    ax1.set_ylim([0.0, 25.])
    ax2.set_ylim([0.0, 25.])
    ax2.set_yticklabels([''])
    # add axis labels
    ax1.set_ylabel('Probability', fontsize = 12)
    ax1.set_xlabel('CO [ppmv]', x = 1.1, fontsize = 12)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'pdf_aqscogmico_%d-%d_%s.eps' %(years[0], years[-1], region), 
                dpi = 300)
    return
# # # # # # # # # # # # #  
def pdf_aqsno2gmino2(aqs_no2, mr2_gmi, ccmi_gmi, ffigac2_gmi, ffigac2hr_gmi, 
                     years, region):
    """function plots probability density functions of NO2 distribution for 
    AQS observations and the four GMI CTM model cases/configurations.
    Concentrations are regionally-averaged over the region used when 
    finding commensurate observed/modeled NO2 concentrations. 
    
    Parameters
    ----------    
    aqs_no2 : numpy.ndarray
        AQS NO2 observations in region, units of ppbm, [years in measuring 
        period * days in months in 'sampling_months',]        
    mr2_gmi : numpy.ndarray
        GMI CTM NO2 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ccmi_gmi : numpy.ndarray
        GMI CTM NO2 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-CCMI, units 
        of ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2_gmi : numpy.ndarray
        GMI CTM NO2 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2, units of 
        ppbv, [years in measuring period * days in months in 
        'sampling_months',]      
    ffigac2hr_gmi : numpy.ndarray 
        GMI CTM NO2 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastFFIgac2-HighRes, 
        units of ppbv, [years in measuring period * days in months in 
        'sampling_months',]       
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
    import matplotlib.mlab as mlab
    import matplotlib.font_manager as font_manager    
    import sys
    sys.path.append('/Users/ghkenrr/phd/')
    import pollutants_constants
    # set custom font
    path = pollutants_constants.FONTPATH_LIGHT
    prop = font_manager.FontProperties(fname=path)
    mpl.rcParams['font.family'] = prop.get_name()
    path = pollutants_constants.FONTPATH_BOLD
    prop = font_manager.FontProperties(fname = path)
    mpl.rcParams['mathtext.bf'] = prop.get_name()
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    ax1.set_title('MERRA-1', fontsize = 12)
    ax2.set_title('MERRA-2', fontsize = 12)
    # obs
    n, bins, patches = ax1.hist(np.hstack(aqs_no2), 25, normed = 1, alpha = 0.)
    pdf_aqs = mlab.normpdf(bins, np.mean(aqs_no2), np.std(aqs_no2))
    ax1.plot(bins, pdf_aqs, lw = 2., color = 'k', label = 'AQS')
    ax2.plot(bins, pdf_aqs, lw = 2., color = 'k', label = 'AQS')
    ## HindcastFFIgac2
    n, bins, patches = ax1.hist(np.hstack(ffigac2_gmi), 25, normed = 1, alpha = 0.)
    pdf_ffigac2 = mlab.normpdf(bins, np.mean(ffigac2_gmi), np.std(ffigac2_gmi))
    ax1.plot(bins, pdf_ffigac2, lw = 1.5, 
             color = pollutants_constants.COLOR_FFIGAC2, label = 'FFIgac2')
    # HindcastFFIgac2-HighRes
    n, bins, patches = ax1.hist(np.hstack(ffigac2hr_gmi), 25, normed = 1, alpha = 0.)
    pdf_ffigac2hr = mlab.normpdf(bins, np.mean(ffigac2hr_gmi), np.std(ffigac2hr_gmi))
    ax1.plot(bins, pdf_ffigac2hr, lw = 1.5, 
             color = pollutants_constants.COLOR_FFIGAC2HIGHRES, 
             label = 'FFIgac2-HighRes')
    ax1.legend(loc = 'upper left', frameon = False, fontsize = 10)
    # HindcastMR2
    n, bins, patches = ax2.hist(np.hstack(mr2_gmi), 25, normed = 1, alpha = 0.)
    pdf_mr2 = mlab.normpdf(bins, np.mean(mr2_gmi), np.std(mr2_gmi))
    ax2.plot(bins, pdf_mr2, lw = 1.5, 
             color = pollutants_constants.COLOR_MR2, label = 'MR2')   
    # HindcastMR2-CCMI
    n, bins, patches = ax2.hist(np.hstack(ccmi_gmi), 25, normed = 1, alpha = 0.)
    pdf_ccmi = mlab.normpdf(bins, np.mean(ccmi_gmi), np.std(ccmi_gmi))
    ax2.plot(bins, pdf_ccmi, lw = 1.5, 
             color = pollutants_constants.COLOR_MR2CCMI, label = 'MR2-CCMI')    
    ax2.legend(loc = 'upper left', frameon = False, fontsize = 10)
    ax1.set_xlim([0, 8])
    ax2.set_xlim([0, 8])
    ax1.set_ylim([0.0, 6.])
    ax2.set_ylim([0.0, 6.])
    ax2.set_yticklabels([''])
    # add axis labels
    ax1.set_ylabel('Probability', fontsize = 12)
    ax1.set_xlabel('NO2 [ppbv]', x = 1.1, fontsize = 12)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'pdf_aqsno2gmino2_%d-%d_%s.eps' %(years[0], years[-1], region), 
                dpi = 300)
    return
# # # # # # # # # # # # #
def pdf_allgmio3(castnet, mr2, ccmi, egu, dat, mechanism, region):
    """function plots PDFs of regionally-averaged O3 for different GMI 
    simulations.
    
    Parameters
    ----------
    castnet : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period, days in months in 'sampling_months']   
    mr2 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    ccmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-CCMI, units 
        of ppbv, [years in measuring period, days in months in 
        'sampling_months']         
    egu : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2 with daily-
        varying NO emissions, units of ppbv, [years in measuring period, days 
        in months in 'sampling_months']                 
    dat : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case haa, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']          
    mechanism : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMERRA, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months']                  
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
        
    Returns
    ----------      
    None     
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.stats
    castnet = np.hstack(castnet)
    mr2 = np.hstack(mr2)
    ccmi = np.hstack(ccmi)
    egu = np.hstack(egu)
    dat = np.hstack(dat)
    mechanism = np.hstack(mechanism)
    # initialize figure, axes
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    # obs
    n, bins, patches = ax.hist(np.hstack(castnet), 25, density = 1, alpha = 0.)
    pdf_castnet = scipy.stats.norm.pdf(bins, np.mean(castnet), np.std(castnet))
    ax.plot(bins, pdf_castnet, lw = 2.5, color = 'k', label = 'Observations')
    # HindcastMR2
    n, bins, patches = ax.hist(np.hstack(castnet), 25, density = 1, alpha = 0.)
    pdf_mr2 = scipy.stats.norm.pdf(bins, np.mean(mr2), np.std(mr2))
    ax.plot(bins, pdf_mr2, lw = 1.5, color = '#e41a1c', label = 'Control')    
    # HindcastMR2-CCMI
    n, bins, patches = ax.hist(np.hstack(ccmi), 25, density = 1, alpha = 0.)
    pdf_ccmi = scipy.stats.norm.pdf(bins, np.mean(ccmi), np.std(ccmi))
    ax.plot(bins, pdf_ccmi, lw = 1.5, color = '#377eb8', label = 'Resolution') 
    # EGU_T
    n, bins, patches = ax.hist(np.hstack(egu), 25, density = 1, alpha = 0.)
    pdf_egu = scipy.stats.norm.pdf(bins, np.mean(egu), np.std(egu))
    ax.plot(bins, pdf_egu, lw = 1.5, color = '#4daf4a', label = 'Emissions')
    # HindcastMERRA        
    n, bins, patches = ax.hist(np.hstack(mechanism), 25, density = 1, alpha = 0.)
    pdf_egu = scipy.stats.norm.pdf(bins, np.mean(mechanism), np.std(mechanism))
    ax.plot(bins, pdf_egu, lw = 1.5, color = '#984ea3', label = 'Mechanism')        
    # HindcastMR2-DiurnalAvgT
    n, bins, patches = ax.hist(np.hstack(dat), 25, density = 1, alpha = 0.)
    pdf_mad = scipy.stats.norm.pdf(bins, np.mean(dat), np.std(dat))
    ax.plot(bins, pdf_mad, lw = 1.5, color = '#ff7f00', label = 'Chemistry')
    ax.set_xlabel('Ozone [ppbv]', fontsize = 12)
    ax.set_ylabel('Probability [ppbv$^{\mathregular{-1}}$]', fontsize = 12)
    plt.subplots_adjust(top = 0.85)
    ax.legend(ncol = 3, frameon = False, fontsize = 10, 
              bbox_to_anchor = [0.93, 1.2])
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'pdf_allgmio3_%s.eps' %(region), dpi = 300)
    return 
# # # # # # # # # # # # #
def scatter_allgmio3(castnet, mr2, ccmi, egu, dat, mechanism, t2m_ra, region):
    """function plots the relationship of O3 from regionally-averaged CASTNet
    observations and co-located GMI output versus MERRA-2 2-meter temperatures. 
    O3 from GMI is plotted for all simulations. The lines of best fit are found
    and plotted along with the slope of these lines. 
    
    Parameters
    ----------
    castnet : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period, days in months in 'sampling_months']   
    mr2 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    ccmi : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-CCMI, units 
        of ppbv, [years in measuring period, days in months in 
        'sampling_months']         
    egu : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2 with daily-
        varying NO emissions, units of ppbv, [years in measuring period, days 
        in months in 'sampling_months']                 
    dat : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']          
    mechanism : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMERRA, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    t2m_ra : numpy.ndarray
        Regionally-averaged MERRA-2 2-meter temperatures co-located (or nearly 
        colocated) with corresponding CASTNet stations, units of K, [years in 
        measuring period, days in months in 'sampling_months']   
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
        
    Returns
    ----------      
    None         
    """
    import numpy as np
    import matplotlib.pyplot as plt
    castnet = np.hstack(castnet)
    mr2 = np.hstack(mr2)
    ccmi = np.hstack(ccmi)
    egu = np.hstack(egu)
    dat = np.hstack(dat)
    mechanism = np.hstack(mechanism)    
    # initialize figure, axes
    fig = plt.figure()
    ax = plt.subplot2grid((4, 4), (0, 0), rowspan = 4, colspan = 3)
    # obs
    ax.plot(t2m_ra, castnet, marker = 'o', lw = 0.0, color = 'k', 
            markersize = 3, label = 'Observations')
    # HindcastMR2-CCMI
    ax.plot(t2m_ra, ccmi, marker = 'o', lw = 0.0, color = '#377eb8', 
            markersize = 2, label = 'Resolution')
    # HindcastMERRA        
    ax.plot(t2m_ra, mechanism , marker = 'o', lw = 0.0, color = '#984ea3', 
            markersize = 2, label = 'Mechanism')
    # HindcastMR2
    ax.plot(t2m_ra, mr2, marker = 'o', lw = 0.0, color = '#e41a1c', 
            markersize = 2, label = 'Control')
    # EGU_T
    ax.plot(t2m_ra, egu, marker = 'o', lw = 0.0, color = '#4daf4a', 
            markersize = 2, label = 'Emissions')
    # HindcastMR2-DiurnalAvgT
    ax.plot(t2m_ra, dat, marker = 'o', lw = 0.0, color = '#ff7f00', 
            markersize = 2, label = 'Chemistry')
    ax.set_xlabel('T$_{\mathregular{2 m}}$ [K]', fontsize = 12, x = 0.75)
    ax.set_ylabel('Ozone [ppbv]', fontsize = 12)
    # axis for lines of best fit
    axr = plt.subplot2grid((4, 4), (0, 3), rowspan = 4, colspan = 1)
    # obs
    castnet_m = np.poly1d(np.polyfit(t2m_ra, castnet, 1))[1]
    axr.plot(np.unique(t2m_ra), np.poly1d(np.polyfit(t2m_ra, castnet, 
             1))(np.unique(t2m_ra)), color = 'k', lw = 2.5)
    # HindcastMR2
    mr2_m = np.poly1d(np.polyfit(t2m_ra, mr2, 1))[1]
    axr.plot(np.unique(t2m_ra), np.poly1d(np.polyfit(t2m_ra, mr2, 
             1))(np.unique(t2m_ra)), color = '#e41a1c')
    # HindcastMR2-CCMI
    ccmi_m = np.poly1d(np.polyfit(t2m_ra, ccmi, 1))[1]
    axr.plot(np.unique(t2m_ra), np.poly1d(np.polyfit(t2m_ra, ccmi, 
             1))(np.unique(t2m_ra)), color = '#ff7f00')
    # EGU_T
    egu_m = np.poly1d(np.polyfit(t2m_ra, ccmi, 1))[1]
    axr.plot(np.unique(t2m_ra), np.poly1d(np.polyfit(t2m_ra, egu, 
             1))(np.unique(t2m_ra)), color = '#4daf4a')
    # HindcastMR2-DiurnalAvgT
    dat_m = np.poly1d(np.polyfit(t2m_ra, dat, 1))[1]
    axr.plot(np.unique(t2m_ra), np.poly1d(np.polyfit(t2m_ra, dat, 
             1))(np.unique(t2m_ra)), color = '#377eb8')
    # HindcastMERRA     
    mechanism_m = np.poly1d(np.polyfit(t2m_ra, mechanism, 1))[1]
    axr.plot(np.unique(t2m_ra), np.poly1d(np.polyfit(t2m_ra, mechanism, 
             1))(np.unique(t2m_ra)), color = '#984ea3')
    # set axis limits to be the same as the larger axis
    axr.set_xlim(ax.get_xlim())
    axr.set_ylim(ax.get_ylim())
    axr.set_yticklabels([''])
    pos = 0.95
    for m, col in zip([castnet_m, mr2_m, ccmi_m, egu_m, dat_m, mechanism_m], 
        ['k', '#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3']):
        axr.text(0.2, pos, '%.2f' %m,
                va = 'center', ha = 'center', transform = axr.transAxes,
                color = col, fontsize = 10)
        pos = pos - 0.06
    # add legend
    plt.subplots_adjust(top = 0.85)
    ax.legend(ncol = 3, frameon = False, fontsize = 10, 
              bbox_to_anchor = [1.28, 1.2])
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatterplot_allgmio3_%s.eps' %(region), dpi = 300)    
    return 
# # # # # # # # # # # # #    
def scatter_t2m_dato3eguo3_highlowdays(mr2, dat, egu, t2m_ra, region):
    """function identifies hot and cold days (> 80th percentile and < 20th 
    percentile) from regionally-averaged MERRA-2 2-meter temperatures and then
    calculates the difference in O3 between the (1) MR2 and EGU-T simulations and
    (2) MR2 and DAT simulations on these hot and cold days and plots the 
    differences as a function of 2-meter temperature. 
    
    Parameters
    ----------
    mr2 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    egu : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2 with daily-
        varying NO emissions, units of ppbv, [years in measuring period, days 
        in months in 'sampling_months']                 
    dat : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']          
    t2m_ra : numpy.ndarray
        Regionally-averaged MERRA-2 2-meter temperatures co-located (or nearly 
        colocated) with corresponding CASTNet stations, units of K, [years in 
        measuring period, days in months in 'sampling_months']   
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
        
    Returns
    ----------      
    None      
    """
    import numpy as np
    import matplotlib.pyplot as plt
    mr2 = np.hstack(mr2)
    dat = np.hstack(dat)
    egu = np.hstack(egu)
    t2m_ra = np.hstack(t2m_ra)
    # find warm and cool days 
    cool = np.where(t2m_ra < np.percentile(t2m_ra, 20))[0]
    warm = np.where(t2m_ra > np.percentile(t2m_ra, 80))[0]
    # initialize figure, axis
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    # for DAT
    ax.plot(t2m_ra[warm], mr2[warm] - dat[warm], 'ko', markersize = 3, 
            label = 'Chemistry') 
    ax.plot(t2m_ra[cool], mr2[cool] - dat[cool], 'ko', markersize = 3)
    # for EGU-T
    ax.plot(t2m_ra[warm], mr2[warm] - egu[warm], 'k+', label = 'Emissions') 
    ax.plot(t2m_ra[cool], mr2[cool] - egu[cool], 'k+')
    # color "warm" and "cool" sectors of plot
    p20 = np.percentile(t2m_ra, 20)
    p80 = np.percentile(t2m_ra, 80)
    ax.set_xlim([t2m_ra.min() - 0.5, t2m_ra.max() + 0.5])
    ax.axvspan(ax.get_xlim()[0], p20, alpha = 0.5, color = '#2166ac')
    ax.axvspan(p80, ax.get_xlim()[1], alpha = 0.5, color = '#b2182b')
    ax.set_xlabel('T$_{\mathregular{2 m}}$ [K]', fontsize = 12)
    ax.set_ylabel('$\mathregular{\Delta}$ Ozone [ppbv]', fontsize = 12)
    ax.text(0.08, 0.9, 'Cool days', va = 'center', ha = 'left', 
            transform = ax.transAxes, fontsize = 16)
    ax.text(0.7, 0.1, 'Warm days', va = 'center', ha = 'left', 
            transform = ax.transAxes, fontsize = 16)
    ax.legend(ncol = 1, frameon = False, fontsize = 10,
              bbox_to_anchor = [0.35, 0.55])
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_t2m_dato3eguo3_highlowdays_%s.eps' %(region), 
                dpi = 300) 
    return 
# # # # # # # # # # # # #
def scatter_matrix_mr2o3meteorology(times, mr2, H500, U500, V500, comm_t2m,
    comm_u10m, comm_v10m, sampling_hours, region):
    """using 6-hourly MERRA-2 output, function finds output in 'sampling_hours'
    and finds regionally-averaged geopotential height and U/V wind components
    at 500 hPa along with regionally-averaged 2-meter temperature and 10-meter
    U/V wind components. A scatterplot matrix with histograms on the diagonal
    is plotted for these variables on all JJA days during 2008-2010 and for 
    only warm days. 
    
    Parameters
    ----------
    times : numpy.ndarray
        Datetime objects corresponding to each 6-hourly timestep of MERRA-2,
        [time,]
    mr2 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    H500 : numpy.ndarray
        Geopotential height at 500 hPa, units of m, [time, lat, lon]
    U500 : numpy.ndarray 
        U wind component at 500 hPa, units of m/s, [time, lat, lon]
    V500 : numpy.ndarray
        V wind component at 500 hPa, units of m/s, [time, lat, lon]        
    mr2 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    comm_t2m : numpy.ndarray
        MERRA-2 2-meter temperatures co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of K, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']
    comm_u2m : numpy.ndarray
        MERRA-2 2-meter eastward wind co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of m s-1, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']    
    comm_u10m : numpy.ndarray
        MERRA-2 2-meter eastward wind co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of m s-1, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']    
    comm_v2m : numpy.ndarray
        MERRA-2 2-meter northward wind co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of m s-1, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']    
    comm_v10m : numpy.ndarray
        MERRA-2 1m-meter northward wind co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of m s-1, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']    
    sampling_hours : list 
        Hours during which trace gas concentrations are returned
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
        
    Returns
    ----------      
    None        
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    # for 6-hourly fields, find fields within 'sampling_hours'
    times = pd.to_datetime(times)
    times_ish = np.where(np.in1d(
            times.hour, np.array(sampling_hours)) == True)[0]
    H500_ish = H500[times_ish]
    U500_ish = U500[times_ish] 
    V500_ish = V500[times_ish] 
    # find regionally-averaged quantities, n.b. if 'sampling_hours' is 
    # expanded to include more hours, H500, U/V500 fields will need to be
    # reshaped to calculate a daily average
    H500_ish = np.nanmean(H500_ish, axis = tuple((1, 2)))
    U500_ish = np.nanmean(U500_ish, axis = tuple((1, 2)))
    V500_ish = np.nanmean(V500_ish, axis = tuple((1, 2)))
    t2m_ra = np.hstack(np.nanmean(comm_t2m, axis = 1))
    u10m_ra = np.hstack(np.nanmean(comm_u10m, axis = 1))
    v10m_ra = np.hstack(np.nanmean(comm_v10m, axis = 1))
    mr2 = np.hstack(mr2)
    # make pandas DataFrame with meteorological/chemistry variables
    df = pd.DataFrame(data = [H500_ish, U500_ish, V500_ish, u10m_ra, v10m_ra, 
                              t2m_ra, mr2])
    df = df.T
    df.columns = ['H$_{\mathregular{500}}$', 'U$_{\mathregular{500}}$', 
                  'V$_{\mathregular{500}}$', 'U$_{\mathregular{10}}$', 
                  'V$_{\mathregular{10}}$',  'T$_{\mathregular{2 m}}$', 
                  'O$_{\mathregular{3}}$']
    fig = plt.figure()
    sm = pd.plotting.scatter_matrix(df, diagonal = 'hist', color = 'k', 
                                    marker = '.')
    [s.get_xaxis().set_label_coords(0.5, -1.) for s in sm.reshape(-1)]
    [s.get_yaxis().set_label_coords(-0.5,0.5) for s in sm.reshape(-1)]
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_matrix_mr2o3meteorology_%s.eps' %(region), 
                dpi = 300)     
    plt.show()    
    # same as above but for hot days 
    warm = np.where(t2m_ra > np.percentile(t2m_ra, 80))[0]
    H500_ish_warm = H500_ish[warm]
    U500_ish_warm = U500_ish[warm]
    V500_ish_warm = V500_ish[warm]
    u10m_ra_warm = u10m_ra[warm]
    v10m_ra_warm = v10m_ra[warm]
    t2m_ra_warm = t2m_ra[warm]
    mr2_warm = mr2[warm] 
    df_warm = pd.DataFrame(data = [H500_ish_warm, U500_ish_warm, V500_ish_warm, 
                                   u10m_ra_warm, v10m_ra_warm, t2m_ra_warm, 
                                   mr2_warm])
    df_warm = df_warm.T
    df_warm.columns = df.columns
    smwarm = pd.plotting.scatter_matrix(df_warm, diagonal = 'hist', color = 'k', 
                                        marker = '.')
    [s.get_xaxis().set_label_coords(0.5, -1.) for s in smwarm.reshape(-1)]
    [s.get_yaxis().set_label_coords(-0.5,0.5) for s in smwarm.reshape(-1)]
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_matrix_mr2o3meteorology_warm_%s.eps' %(region), 
                dpi = 300)         
    return
# # # # # # # # # # # # #    
def timeseries_mr2o3dato3(mr2, dat, region):
    """function plots time series of O3 from MR2 and DAT. The threshold for 
    high (> 80th %-ile) and low (< 20 %-ile) events from the MR2 simulation 
    are plotted and the number of events in both simulations are written in 
    text on the plot, along with the correlation coefficient. 
    
    Parameters
    ----------       
    mr2 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    dat : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
        
    Returns
    ----------      
    None         
    """
    import numpy as np
    import matplotlib.pyplot as plt
    # stack yearly time series into one continuous time series
    mr2 = np.hstack(mr2)
    dat = np.hstack(dat)
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.plot(mr2, '-', color = '#ff7f00', label = 'MR2')
    ax.plot(dat, '-', color = '#377eb8', label = 'DAT')
    # calculate and graph threshold for extreme events
    mr2_80 = np.percentile(mr2, 80)
    mr2_20 = np.percentile(mr2, 20)
    ax.hlines(mr2_80, xmin = ax.get_xlim()[0], xmax = ax.get_xlim()[1])
    ax.hlines(mr2_20, xmin = ax.get_xlim()[0], xmax = ax.get_xlim()[1])
    # find percentage of extreme events 
    where_high_mr2 = np.where(mr2 > mr2_80)[0]
    where_high_dat = np.where(dat > mr2_80)[0]
    where_low_mr2 = np.where(mr2 < mr2_20)[0]
    where_low_dat = np.where(dat < mr2_20)[0]
    # percentage of co-occurring MR2 and DAT high extremes
    both_high = np.in1d(where_high_mr2, where_high_dat)
    both_high = len(np.where(both_high == True)[0])/len(where_high_mr2)
    both_low = np.in1d(where_low_mr2, where_low_dat)
    both_low = len(np.where(both_low == True)[0])/len(where_low_dat)
    # add text corresponding to the number of total high/low events for the
    # MR2 simulation (n = 55, by definition) and the number of events in the 
    # DAT simulation using the thresholds from the MR2 simulation
    ax.text(0.05, 0.95, 'No. high events MR2 = %d' %len(where_high_mr2), 
            transform = ax.transAxes)
    ax.text(0.05, 0.90, 'No. high events DAT = %d' %len(where_high_dat), 
            transform = ax.transAxes)
    ax.text(0.05, 0.85, 'No. low events MR2 = %d' %len(where_low_mr2), 
            transform = ax.transAxes)
    ax.text(0.05, 0.80, 'No. low events DAT = %d' %len(where_low_dat), 
            transform = ax.transAxes)    
    # calculate correlation coefficient and plot
    r = np.corrcoef(mr2, dat)[0, 1]
    ax.set_title('r = %.2f' %r, x = 0.1)
    # aesthetics
    ax.set_xlim([0, len(mr2)])
    # legend
    ax.legend(loc = 4, frameon = False)
    ax.set_xticks([0, 30, 61, 92, 122, 153, 184, 214, 245])
    ax.set_xticklabels(['Jun 2008', 'Jul', 'Aug', 'Jun 2009', 'Jul', 'Aug', 
                         'Jun 2010', 'Jul', 'Aug'])
    ax.set_ylabel('Ozone [ppbv]')
    ax.text(len(mr2) + 6, mr2_80, '80$^{\mathregular{th}}$ %-ile', 
            va = 'center')
    ax.text(len(mr2) + 6, mr2_20, '20$^{\mathregular{th}}$ %-ile', 
            va = 'center')
    plt.subplots_adjust(right = 0.85)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'timeseries_mr2o3dato3_%s.eps' %(region), dpi = 300)     
    return 
# # # # # # # # # # # # #    
def histogram_mr2o3dato3(mr2, dat, region):
    """function plots histograms of O3 from MR2 and DAT simulations and 
    denotes the 20th and 80th percentiles. This was done because PDFs smooth 
    out some of the detail on the tails of the distribution.
    
    Parameters
    ----------       
    mr2 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    dat : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
        
    Returns
    ----------      
    None         
    """
    import numpy as np
    import matplotlib.pyplot as plt
    # stack yearly time series into one continuous time series
    mr2 = np.hstack(mr2)
    dat = np.hstack(dat)
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.hist(mr2, bins = 25, histtype = 'stepfilled', density = True, 
            color = '#ff7f00', label = 'MR2', zorder = 18)
    ax.hist(dat, bins = 25, histtype = 'step', density = True, lw = 3., 
            color = '#377eb8', label = 'DAT', zorder = 20)
    # calculate and graph threshold for extreme events
    mr2_80 = np.percentile(mr2, 80)
    mr2_20 = np.percentile(mr2, 20)
    ymin, ymax = ax.get_ylim()[0], ax.get_ylim()[1]
    ax.vlines(mr2_80, ymin = ymin, ymax = ymax, zorder = 16)
    ax.vlines(mr2_20, ymin = ymin, ymax = ymax, zorder = 16)
    ax.set_ylim([0, 0.095])
    ax.text(mr2_20, ax.get_ylim()[1] + 0.002, '20$^{\mathregular{th}}$ %-ile', 
            ha = 'center')
    ax.text(mr2_80, ax.get_ylim()[1] + 0.002, '80$^{\mathregular{th}}$ %-ile', 
            ha = 'center')
    plt.legend(frameon = False)
    ax.set_xlabel('Ozone [ppbv]')
    ax.set_ylabel('Probability')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'histogram_mr2o3dato3_%s.eps' %(region), dpi = 300)         
    return 
# # # # # # # # # # # # #    
def scatter_t2m_mr2o3dato3_deltat2m_delta_o3(mr2, dat, t2m, region):
    """function plots the relationship between O3 and MERRA-2 2-meter 
    temperature. The left subplot is O3 versus temperature from the MR2 
    and DAT simulations (where the DAT simulation temperature is represented
    by montly mean values). The right subplot is delta O3 versus delta T with 
    the slope is as the subplot's title and a line of best fit plotted. 

    Parameters
    ----------       
    mr2 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    dat : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    t2m : numpy.ndarray    
        Regionally-averaged MERRA-2 2-meter temperatures co-located (or nearly 
        colocated) with corresponding CASTNet stations, units of K, [years in 
        measuring period, days in months in 'sampling_months'] 
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
        
    Returns
    ----------      
    None    
    """
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt
    # calculate regionally-averaged temperature 
    t2m = np.nanmean(t2m, axis = 1)
    # stack yearly time series into one continuous time series
    mr2 = np.hstack(mr2)
    dat = np.hstack(dat)
    t2m = np.hstack(t2m)
    dat_m, t2m_m = [], []
    # for the DAT (fixed temperature) run, find average temperature in each 
    # month 
    md = [0, 30, 61, 92, 122, 153, 184, 214, 245, 276]
    for i, mds in enumerate(md[:-1]):
        dat_m.append(dat[mds:md[i + 1]])
        # mean temperature in month
        t2m_t = t2m[mds:md[i + 1]]
        t2m_m.append(np.repeat(np.mean(t2m_t), len(dat[mds:md[i + 1]])))
    t2m_m = np.hstack(t2m_m)
    dat_m = np.hstack(dat_m)
    fig = plt.figure()
    ax1 = plt.subplot2grid((1, 2), (0, 0))    
    ax1.plot(t2m_m, dat_m, 'o', markersize = 3, color = '#377eb8', 
            label = 'DAT')
    ax1.plot(t2m, mr2, 'o', markersize = 3, color = '#ff7f00', label = 'MR2')
    ## calculate and graph threshold for extreme events
    #mr2_80 = np.percentile(mr2, 80)
    #mr2_20 = np.percentile(mr2, 20)
    ## find percentage of extreme events 
    #where_high_mr2 = np.where(mr2 > mr2_80)[0]
    #where_high_dat = np.where(dat > mr2_80)[0]
    #where_low_mr2 = np.where(mr2 < mr2_20)[0]
    #where_low_dat = np.where(dat < mr2_20)[0]
    ## find and plot high extremes in DAT that are not in MR2
    #high_only_dat = np.where(np.in1d(where_high_dat, where_high_mr2) == 
    #                         False)[0]
    #high_only_dat = where_high_dat[high_only_dat]
    #ax1.plot(t2m_m[high_only_dat], dat_m[high_only_dat], 'o', 
    #        markerfacecolor = None, markeredgecolor = 'k', label = 'DAT event')
    plt.legend(frameon = False)
    ax1.set_xlabel('T$_\mathregular{2m}$ [K]')
    ax1.set_ylabel('Ozone [ppbv]')
    # calculate delta O3 and delta T (where delta = MR2 - DAT)
    dO3 = mr2 - dat
    dT = t2m - t2m_m
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    ax2.plot(dT, dO3, 'ko', markersize = 3)
    ax2.set_xlabel('$\mathregular{\Delta}$ T$_\mathregular{2m}$ [K]')
    ax2.set_ylabel('$\mathregular{\Delta}$ Ozone [ppbv]')
    # find slope and plot line of best fit
    m = np.poly1d(np.polyfit(dT, dO3, 1))[1]
    ax2.plot(np.unique(dT), np.poly1d(np.polyfit(dT, dO3,
             1))(np.unique(dT)), '--', color = 'r', lw = 2.5, zorder = 1)
    ax2.set_title('m = %.2f' %m, x = 0.2)
    plt.tight_layout()
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_t2m_mr2o3dato3_deltat2m_delta_o3_%s.eps' %(region), 
                dpi = 300)          
    return 
# # # # # # # # # # # # #    
def timeseries_mr2o3dato3eguto3(transport, chemistry, emissions, obs, region,
                                year, years):
    """function plots time series of O3 from GMI transport simulations with 
    variability from (1) only transport; (2) transport and chemistry; and (3) 
    transport, chemistry, and emissions. Observations from CASTNet are plotted
    alongside GMI O3.
    
    Parameters
    ----------       
    transport : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    chemistry : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    emissions : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2 with daily-
        varying NO emissions, units of ppbv, [years in measuring period, days 
        in months in 'sampling_months']
    obs : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period, days in months in 'sampling_months']        
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
    # stack yearly time series into one continuous time series
    transport = np.hstack(transport)
    chemistry = np.hstack(chemistry)
    emissions = np.hstack(emissions)
    obs = np.hstack(obs)
    # for single year
    yearpos = np.where(np.array(years) == year)[0][0]    
    fig = plt.figure(figsize = (10, 4))
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.plot(obs[yearpos*92:(yearpos+1)*92], '-', color = 'lightgrey', 
             zorder = 1, lw = 2, label = 'CASTNet')
    ax.plot(emissions[yearpos*92:(yearpos+1)*92], '-', 
            color = pollutants_constants.COLOR_EMISSIONS, 
            label = '+$\:$Emissions')
    ax.plot(chemistry[yearpos*92:(yearpos+1)*92],
            '-', color = pollutants_constants.COLOR_CHEMISTRY, 
            label = '+$\:$Chemistry')
    ax.plot(transport[yearpos*92:(yearpos+1)*92], '-', 
            color = pollutants_constants.COLOR_TRANSPORT, label = 'Transport')
    ax.set_xlim([0, 91])
    # axis labels
    ax.set_xticks([0, 30, 61])
    ax.set_xticklabels(['1 June %d' %year, '1 July', '1 Aug'], ha = 'left', 
                       fontsize = 12)
    ax.set_ylabel('Ozone [ppbv]', fontsize = 16)
    for tl in ax.get_yticklabels():
        tl.set_fontsize(12)
    # reorder and place legend
    handles, labels = ax.get_legend_handles_labels()
    handles = [handles[0], handles[3], handles[2], handles[1]]
    labels = [labels[0], labels[3], labels[2], labels[1]]
    ax.legend(handles,labels, bbox_to_anchor = (0.5, 1.08), loc = 'center', 
              ncol = 4, frameon = False, fontsize = 16)
    print('CASTNet-DAT correlation = %.3f' %(
            np.corrcoef(obs[yearpos*92:(yearpos+1)*92], 
            transport[yearpos*92:(yearpos+1)*92])[0, 1]))
    print('CASTNet-MR2 correlation = %.3f' %(
            np.corrcoef(obs[yearpos*92:(yearpos+1)*92], 
            chemistry[yearpos*92:(yearpos+1)*92])[0, 1]))
    print('CASTNet-EGU correlation = %.3f' %(
            np.corrcoef(obs[yearpos*92:(yearpos+1)*92], 
            emissions[yearpos*92:(yearpos+1)*92])[0, 1]))     
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'timeseries_mr2o3dato3eguto3_%d_%s.eps' %(year, region), 
                dpi = 300)
#    # for all three years
#    f, (ax, ax2, ax3) = plt.subplots(1, 3, sharey = True, figsize=(10, 6))
#    for axt in [ax, ax2, ax3]:    
#        axt.plot(obs, '-', color = 'lightgrey', zorder = 1, lw = 2,
#                 label = 'CASTNet')
#        axt.plot(emissions, '-', color = pollutants_constants.COLOR_EMISSIONS, 
#                 label = '+$\:$Emissions')
#        axt.plot(chemistry, '-', color = pollutants_constants.COLOR_CHEMISTRY, 
#                 label = '+$\:$Chemistry')
#        axt.plot(transport, '-', color = pollutants_constants.COLOR_TRANSPORT, 
#                 label = 'Transport')
#    ax.set_xlim([0, 92 + 1])
#    ax2.set_xlim([92, 92*2 + 1])
#    ax3.set_xlim([92*2, 92*3])
#    # hide the spines between axes
#    ax.spines['right'].set_visible(False)
#    ax2.spines['left'].set_visible(False)
#    ax2.spines['right'].set_visible(False)
#    ax3.spines['left'].set_visible(False)
#    ax.yaxis.tick_left()
#    ax2.tick_params(axis = 'y', which = 'both', left = False, right = False)
#    ax3.yaxis.tick_right()
#    # add discontinuous lines to axis (diagonal lines)
#    d = .015
#    kwargs = dict(transform = ax.transAxes, color='k', lw = 0.75, clip_on=False)
#    ax.plot((1-d,1+d), (-d,+d), **kwargs)
#    ax.plot((1-d,1+d),(1-d,1+d), **kwargs)
#    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
#    ax2.plot((-d,+d), (1-d,1+d), **kwargs)
#    ax2.plot((-d,+d), (-d,+d), **kwargs)
#    ax2.plot((1-d,1+d), (-d,+d), **kwargs)
#    ax2.plot((1-d,1+d),(1-d,1+d), **kwargs)
#    kwargs.update(transform=ax3.transAxes)
#    ax3.plot((-d,+d), (1-d,1+d), **kwargs)
#    ax3.plot((-d,+d), (-d,+d), **kwargs)
#    # move axes closer to each other
#    f.subplots_adjust(wspace = 0.05)
#    # axis labels
#    ax.set_xticks([0, 30, 61])
#    ax.set_xticklabels(['Jun', 'Jul', 'Aug'], ha = 'left', fontsize = 12)
#    ax.set_xlabel('2008', x = 0., ha = 'left', fontsize = 14)
#    ax2.set_xticks([92, 122, 153])
#    ax2.set_xticklabels(['Jun', 'Jul', 'Aug'], ha = 'left', fontsize = 12)
#    ax2.set_xlabel('2009', x = 0., ha = 'left', fontsize = 14)
#    ax3.set_xticks([184, 214, 245])
#    ax3.set_xticklabels(['Jun', 'Jul', 'Aug'], ha = 'left', fontsize = 12)
#    ax3.set_xlabel('2010', x = 0., ha = 'left', fontsize = 16)
#    ax.set_ylabel('Ozone [ppbv]', fontsize = 16)
#    for tl in ax.get_yticklabels():
#        tl.set_fontsize(12)
#    # remove axis tick labels overlap discontinuous lines
#    for axt in [ax2, ax3]:
#        i = 0
#        for t in axt.xaxis.get_ticklines(): 
#            if i == 0:
#                t.set_visible(False) 
#            i = i + 1
#    # reorder and place legend
#    handles, labels = ax.get_legend_handles_labels()
#    handles = [handles[0], handles[3], handles[2], handles[1]]
#    labels = [labels[0], labels[3], labels[2], labels[1]]
#    ax.legend(handles,labels, bbox_to_anchor = (1.55, 1.05), loc = 'center', 
#              ncol = 4, frameon = False, fontsize = 16)
#    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
#                'timeseries_mr2o3dato3eguto3_%s.eps' %(region), 
#                dpi = 300)   
    return
# # # # # # # # # # # # #    
def scatterhist_mr2o3dato3eguto3(transport, chemistry, emissions, obs, t2m, 
                                 region):
    """function plots a scatterplot of regionally-averaged O3 versus MERRA-2
    2-meter temperatures for the transport, chemistry, and emissions 
    simulations and observations. KDEs of distributions are plotted along the 
    spines of the scatterplot. 
    
    Parameters
    ----------       
    transport : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    chemistry : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    emissions : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2 with daily-
        varying NO emissions, units of ppbv, [years in measuring period, days 
        in months in 'sampling_months']
    obs : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period, days in months in 'sampling_months']   
    t2m : numpy.ndarray
        MERRA-2 2-meter temperatures co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of K, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']       
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
        
    Returns
    ----------      
    None     
    """
    import numpy as np
    import scipy.stats as stats
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.ticker import NullFormatter
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['xtick.minor.width'] = 1.5
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams['ytick.minor.width'] = 1.5    
    t2m = np.nanmean(t2m, axis = 1)
    t2m = np.hstack(t2m)
    transport = np.hstack(transport)
    chemistry = np.hstack(chemistry)
    emissions = np.hstack(emissions)
    obs = np.hstack(obs)
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
    axScatter.scatter(t2m, obs, label = 'CASTNet', s = 20, 
                      color = 'lightgrey', zorder = 2)    
    axScatter.scatter(t2m, transport, label = 'Transport', s = 20, 
                      color = pollutants_constants.COLOR_TRANSPORT, zorder = 9)
    axScatter.scatter(t2m, chemistry, label = '+$\:$Chemistry', s = 20, 
                      color = pollutants_constants.COLOR_CHEMISTRY, zorder = 6)
    axScatter.scatter(t2m, emissions, label = '+$\:$Emissions', s = 20, 
                      color = pollutants_constants.COLOR_EMISSIONS, zorder = 3)
    # add in lines of best fit 
    obs_m = np.poly1d(np.polyfit(t2m, obs, 1))[1]
    transport_m = np.poly1d(np.polyfit(t2m, transport, 1))[1]
    chemistry_m = np.poly1d(np.polyfit(t2m, chemistry, 1))[1]
    emissions_m = np.poly1d(np.polyfit(t2m, emissions, 1))[1]
    print('m_obs = ', '%.3f' %obs_m)
    print('m_transport = ', '%.3f' %transport_m)
    print('m_chemistry = ', '%.3f' %chemistry_m)
    print('m_emissions = ', '%.3f' %emissions_m)
    axScatter.plot(np.unique(t2m), np.poly1d(np.polyfit(t2m, obs, 
        1))(np.unique(t2m)), zorder = 30, color = 'lightgrey', lw = 2.5)
    axScatter.plot(np.unique(t2m), np.poly1d(np.polyfit(t2m, transport, 
        1))(np.unique(t2m)), zorder = 30,
        color = pollutants_constants.COLOR_TRANSPORT, lw = 2.5)
    axScatter.plot(np.unique(t2m), np.poly1d(np.polyfit(t2m, chemistry, 
        1))(np.unique(t2m)), zorder = 20,
        color = pollutants_constants.COLOR_CHEMISTRY, lw = 2.5)
    axScatter.plot(np.unique(t2m), np.poly1d(np.polyfit(t2m, emissions, 
        1))(np.unique(t2m)), zorder = 10, 
        color = pollutants_constants.COLOR_EMISSIONS, lw = 2.5)    
    axScatter.set_xlabel('T$_{\mathregular{2\:m}}$ [K]', fontsize = 16)
    axScatter.set_ylabel('Ozone [ppbv]', fontsize = 16)
    for tl in axScatter.get_xticklabels():
        tl.set_fontsize(12)
    for tl in axScatter.get_yticklabels():
        tl.set_fontsize(12)    
    axScatter.get_xaxis().set_label_coords(0.5, -0.07)        
    axScatter.get_yaxis().set_label_coords(-0.09, 0.5)        
    # histograms
    nbins = 18
    n, x, _  = axHistx.hist(t2m, bins = nbins, histtype = u'step', 
                            density = True, color = 'k', lw = 0.)
    density = stats.gaussian_kde(t2m)
    axHistx.plot(x, density(x), lw = 2., color = 'k')
    axHistx.set_ylabel('Density', fontsize = 16)
    axHistx.text(0.25, 0.7, 'T$_{\mathregular{2\:m}}$', transform = axHistx.transAxes, 
                 fontsize = 20)    
    axHistx.get_yaxis().set_label_coords(-0.09, 0.5)        
    # observations
    n, x, _  = axHisty.hist(obs, bins = nbins, histtype = u'step', 
                            orientation = 'horizontal', density = True, lw = 0.)
    density = stats.gaussian_kde(obs)
    axHisty.plot(density(x), x, zorder = 2, lw = 2.,color = 'lightgray')
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
    axHisty.text(0.25, 0.8, 'Ozone', transform = axHisty.transAxes, 
                 fontsize = 20)
    axHisty.get_xaxis().set_label_coords(0.5, -0.07)        
    axHistx.set_xlim(axScatter.get_xlim())
    axHistx.spines['top'].set_visible(False)
    axHistx.spines['right'].set_visible(False)
    axHisty.set_ylim(axScatter.get_ylim())
    axHisty.spines['top'].set_visible(False)
    axHisty.spines['right'].set_visible(False)
    axScatter.legend(ncol = 1, frameon = False, fontsize = 16)
    for ax in [axScatter, axHistx, axHisty]:
        for tl in ax.get_xticklabels():
            tl.set_fontsize(12)
        for tl in ax.get_yticklabels():
            tl.set_fontsize(12)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatterhist_mr2o3dato3eguto3_%s.eps' %(region), dpi = 300)     
    return    
# # # # # # # # # # # # #    
def scatter_dt2m_dmr2o3dato3eguto3(transport, chemistry, emissions, t2m, region):
    """function calculates the change in temperature between Transport and 
    + Chemistry simulations (that is, the difference between daily 2-meter
    temperature and monthly mean 2-meter temperatures). The difference in O3 
    between + Chemistry and Transport is plotted against this temperature 
    difference with the line of best fit.
    
    Parameters
    ----------       
    transport : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    chemistry : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    emissions : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2 with daily-
        varying NO emissions, units of ppbv, [years in measuring period, days 
        in months in 'sampling_months']
    t2m : numpy.ndarray
        MERRA-2 2-meter temperatures co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of K, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']       
    region : str
        Region over which regionally-averaged concentrations are supplied to 
        function
        
    Returns
    ----------      
    None         
    """
    import numpy as np
    import matplotlib.pyplot as plt
    # calculate regionally-averaged T2m
    t2m = np.nanmean(t2m, axis = 1)
    # stack along time axis
    chemistry = np.hstack(chemistry)
    emissions = np.hstack(emissions)
    transport = np.hstack(transport)
    t2m = np.hstack(t2m)
    # find monthly mean temperatures 
    t2m_m = []
    month_idx = [0, 30, 61, 92, 122, 153, 184, 214, 245, 276]
    for i, mds in enumerate(month_idx[:-1]):
        # mean temperature in month
        t2m_t = t2m[mds:month_idx[i + 1]]
        t2m_m.append(np.repeat(np.mean(t2m_t), len(t2m_t)))
    t2m_m = np.hstack(t2m_m)
    dO3 = chemistry - transport 
    # initialize figure, axis 
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.scatter((t2m - t2m_m), dO3, s = 20, color = 'k')
    ax.set_xlabel('T$_{\mathregular{2\:m}}$ - <T$_{\mathregular{2\:m}}$> [K]', 
                  fontsize = 14)
    ax.set_ylabel('O$_{\mathregular{3, +\:Chemistry}}$ - O$_{\mathregular{3, Transport}}$ [ppbv]', 
                  fontsize = 14)
    for tl in ax.get_xticklabels():
        tl.set_fontsize(12)
    for tl in ax.get_yticklabels():
        tl.set_fontsize(12)
    # calculate best fit line
    m = np.poly1d(np.polyfit(t2m - t2m_m, dO3, 1))[1]
    ax.plot(np.unique(t2m - t2m_m), np.poly1d(np.polyfit(t2m - 
            t2m_m, dO3, 1))(np.unique(t2m - t2m_m)), color = 'k', lw = 2.5)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 'scatter_dt2m_do3_%s.eps' 
                %(region), dpi = 300)     
    return 
# # # # # # # # # # # # #    
def meano3_t2mpercentile(transport, chemistry, emissions, obs, t2m):
    """function determines on which summer days the regionally-averaged 
    temperatures fall into 0-10th, 10-20th, 45-55th, 80-90th, and 90-100th 
    percentiles, and - on these days - calculates the average O3 from GMI 
    transport, chemistry, and emissions simulations and observations on 
    these days.
    
    Parameters
    ----------       
    transport : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    chemistry : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    emissions : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2 with daily-
        varying NO emissions, units of ppbv, [years in measuring period, days 
        in months in 'sampling_months']
    obs : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period, days in months in 'sampling_months']        
    t2m : numpy.ndarray
        MERRA-2 2-meter temperatures co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of K, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']       
        
    Returns
    ----------      
    None     
    """
    import numpy as np
    # calculate regionally-averaged T2m
    t2m = np.nanmean(t2m, axis = 1)
    # stack along time axis
    obs = np.hstack(obs)
    chemistry = np.hstack(chemistry)
    emissions = np.hstack(emissions)
    transport = np.hstack(transport)
    t2m = np.hstack(t2m)
    # low, mid-range, and high percentiles for T2m 
    p010 = np.where((t2m > np.percentile(t2m, 0)) & 
                    (t2m <= np.percentile(t2m, 10)))[0]
    p1020 = np.where((t2m > np.percentile(t2m, 10)) & 
                    (t2m <= np.percentile(t2m, 20)))[0]
    p4555 = np.where((t2m > np.percentile(t2m, 45)) & 
                    (t2m <= np.percentile(t2m, 55)))[0]
    p8090 = np.where((t2m > np.percentile(t2m, 80)) & 
                    (t2m <= np.percentile(t2m, 90)))[0]
    p90100 = np.where((t2m > np.percentile(t2m, 90)) & 
                    (t2m <= np.percentile(t2m, 100)))[0]
    # O3 for 0th %-ile < T2m <= 10th %-ile
    print('Obs., 0-10 %-ile = ', '%.3f' %(obs[p010].mean()))
    print('Transport, 0-10 %-ile = ', '%.3f' %(transport[p010].mean()))
    print('Chemistry, 0-10 %-ile = ', '%.3f' %(chemistry[p010].mean()))
    print('Emissions, 0-10 %-ile = ', '%.3f' %(emissions[p010].mean()))
    # for 10th %-ile < T2m <= 20th %-ile
    print('Obs., 10-20 %-ile = ', '%.3f' %(obs[p1020].mean()))
    print('Transport, 10-20 %-ile = ', '%.3f' %(transport[p1020].mean()))
    print('Chemistry, 10-20 %-ile = ', '%.3f' %(chemistry[p1020].mean()))
    print('Emissions, 10-20 %-ile = ', '%.3f' %(emissions[p1020].mean()))
    # for 45th %-ile < T2m <= 55th %-ile
    print('Obs., 45-55 %-ile = ', '%.3f' %(obs[p4555].mean()))
    print('Transport, 45-55 %-ile = ', '%.3f' %(transport[p4555].mean()))
    print('Chemistry, 45-55 %-ile = ', '%.3f' %(chemistry[p4555].mean()))
    print('Emissions, 45-55 %-ile = ', '%.3f' %(emissions[p4555].mean()))
    # for 80th %-ile < T2m <= 90th %-ile
    print('Obs., 80-90 %-ile = ', '%.3f' %(obs[p8090].mean()))
    print('Transport, 80-90 %-ile = ', '%.3f' %(transport[p8090].mean()))
    print('Chemistry, 80-90 %-ile = ', '%.3f' %(chemistry[p8090].mean()))
    print('Emissions, 80-90 %-ile = ', '%.3f' %(emissions[p8090].mean()))
    # for 90th %-ile < T2m <= 100th %-ile
    print('Obs., 90-100 %-ile = ', '%.3f' %(obs[p90100].mean()))
    print('Transport, 90-100 %-ile = ', '%.3f' %(transport[p90100].mean()))
    print('Chemistry, 90-100 %-ile = ', '%.3f' %(chemistry[p90100].mean()))
    print('Emissions, 90-100 %-ile = ', '%.3f' %(emissions[p90100].mean()))
    return 
# # # # # # # # # # # # #    
def scatter_dt2m_o3(transport, chemistry, emissions, t2m, region):
    """function calculates the change in temperature between Transport and 
    + Chemistry simulations (that is, the difference between daily 2-meter
    temperature and monthly mean 2-meter temperatures). O3 from GMI 
    Transport, + Chemistry, and + Emissions simulations are plotted against
    this change in temperature along with the line of best fit. 
    
    Parameters
    ----------       
    transport : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    chemistry : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months'] 
    emissions : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2 with daily-
        varying NO emissions, units of ppbv, [years in measuring period, days 
        in months in 'sampling_months']
    t2m : numpy.ndarray
        MERRA-2 2-meter temperatures co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of K, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']       
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
    # calculate regionally-averaged T2m
    t2m = np.nanmean(t2m, axis = 1)
    # stack along time axis
    chemistry = np.hstack(chemistry)
    emissions = np.hstack(emissions)
    transport = np.hstack(transport)
    t2m = np.hstack(t2m)
    # find monthly mean temperatures 
    t2m_m = []
    month_idx = [0, 30, 61, 92, 122, 153, 184, 214, 245, 276]
    for i, mds in enumerate(month_idx[:-1]):
        # mean temperature in month
        t2m_t = t2m[mds:month_idx[i + 1]]
        t2m_m.append(np.repeat(np.mean(t2m_t), len(t2m_t)))
    t2m_m = np.hstack(t2m_m)
    # initialize figure, axis 
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.scatter((t2m - t2m_m), transport, label = 'Transport', s = 20, 
               color = pollutants_constants.COLOR_TRANSPORT, zorder = 9, 
               alpha = 0.4)
    ax.scatter((t2m - t2m_m),chemistry, label = '+$\:$Chemistry', s = 20, 
               color = pollutants_constants.COLOR_CHEMISTRY, zorder = 6, 
               alpha = 0.4)
    ax.scatter((t2m - t2m_m), emissions, label = '+$\:$Emissions', s = 20, 
               color = pollutants_constants.COLOR_EMISSIONS, zorder = 3, 
               alpha = 0.4)
    # calculate and plot best fit lines 
    transport_m = np.poly1d(np.polyfit(t2m - t2m_m, transport, 1))[1]
    chemistry_m = np.poly1d(np.polyfit(t2m - t2m_m, chemistry, 1))[1]
    emissions_m = np.poly1d(np.polyfit(t2m - t2m_m, emissions, 1))[1]
    print('m_transport = ', '%.3f' %transport_m)
    print('m_chemistry = ', '%.3f' %chemistry_m)
    print('m_emissions = ', '%.3f' %emissions_m)
    ax.plot(np.unique(t2m - t2m_m), np.poly1d(np.polyfit(t2m - t2m_m, transport, 
             1))(np.unique(t2m - t2m_m)), zorder = 30,
             color = pollutants_constants.COLOR_TRANSPORT, lw = 2.5)
    ax.plot(np.unique(t2m - t2m_m), np.poly1d(np.polyfit(t2m - t2m_m, chemistry, 
             1))(np.unique(t2m - t2m_m)), zorder = 20,
             color = pollutants_constants.COLOR_CHEMISTRY, lw = 2.5)
    ax.plot(np.unique(t2m - t2m_m), np.poly1d(np.polyfit(t2m - t2m_m, emissions, 
             1))(np.unique(t2m - t2m_m)), zorder = 10,
             color = pollutants_constants.COLOR_EMISSIONS, lw = 2.5)
    ax.set_xlabel('T$_{\mathregular{2\:m}}$ - <T$_{\mathregular{2\:m}}$> [K]', 
                  fontsize = 14)
    ax.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 14)
    for tl in ax.get_xticklabels():
        tl.set_fontsize(12)
    for tl in ax.get_yticklabels():
        tl.set_fontsize(12)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 'scatter_dt2m_o3_%s.eps' 
                %(region), dpi = 300)     
    return 
# # # # # # # # # # # # #    
def timeseries_t2m_castneto3_cemsnox(castnet, t2m, year, years, 
    sampling_months, region):
    """for a given summer (JJA) function plots regionally-averaged 2-meter
    temperatures from MERRA-2, regionally-averaged O3 from CASTNet, and 
    regionally-summed NOx emissions from CEMS and determines the correlation 
    coefficients between these variables. 

    Parameters
    ----------       
    cast et : numpy.ndarray
        CASTNet O3 observations in region, units of ppbv, [years in measuring 
        period, days in months in 'sampling_months']   
    t2m : numpy.ndarray
        MERRA-2 2-meter temperatures co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of K, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']  
    year : int
        Year of interest
    years : list
        Years in measuring period
    sampling_months : list 
        Months of interest
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
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/emissions/')
    import AQSCEMSobs      
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['xtick.minor.width'] = 1.5
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams['ytick.minor.width'] = 1.5
    # load CEMS NOx emissions in NEUS
    states_ab = ['CT', 'DC', 'DE', 'MA', 'MD', 'ME', 'NH', 'NJ', 'NY', 'PA', 
                 'RI', 'VA', 'VT', 'WV']
    nox_state, nox_lat, nox_lon = AQSCEMSobs.cems_specifystates_dailymean(
            '/Volumes/GAIGEKERR/emissions/CEMS/', states_ab, sampling_months)
    nox = nox_state['%d-0%d-01'%(year, sampling_months[0]):
                    '%d-0%d-31'%(year, sampling_months[-1])].values
    # select year of interest
    yearpos = np.where(np.array(years) == year)[0][0]
    t2m = t2m[yearpos]
    castnet = castnet[yearpos]
    # calculate regionally-averaged temperature 
    t2m = np.nanmean(t2m, axis = 0)
    # initialize figure, axes
    fig = plt.figure(figsize = (10, 6))
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan = 3)
    ax2 = plt.subplot2grid((3, 3), (1, 0), colspan = 3)
    ax3 = plt.subplot2grid((3, 3), (2, 0), colspan = 3)
    # MERRA-2 2-meter temperature plot
    ax1.text(0.03, 0.84, '(a)', ha = 'center', va = 'center', 
             transform = ax1.transAxes, fontsize = 20)
    ax1.plot(t2m, lw = 2., color = '#ff7f00')
    ax1.set_xlim([0, len(t2m) - 1])
    ax1.set_xticks([0, 30, 61])
    ax1.set_xticklabels([''])
    for t in ax1.get_yticklabels():
        t.set_fontsize(12)        
    ax1.get_yaxis().set_label_coords(-0.07, 0.5)
    ax1.set_ylabel('T$_{\mathregular{2\:m}}$ [K]', fontsize = 16)
    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_ticks_position('both')
    # CASTNet O3
    ax2.text(0.03, 0.84, '(b)', ha = 'center', va = 'center', 
             transform = ax2.transAxes, fontsize = 20)    
    ax2.plot(castnet, lw = 2., color = 'k')
    ax2.set_xlim([0, len(castnet) - 1])
    ax2.set_xticks([0, 30, 61])
    for t in ax2.get_yticklabels():
        t.set_fontsize(12)
    ax2.set_xticklabels([''])
    ax2.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    ax2.get_yaxis().set_label_coords(-0.07, 0.5)
    ax2.xaxis.set_ticks_position('both')
    ax2.yaxis.set_ticks_position('both')
    # CEMS NOx plot
    ax3.text(0.03, 0.84, '(c)', ha = 'center', va = 'center', 
             transform = ax3.transAxes, fontsize = 20)    
    ax3.plot(nox, lw = 2., color = '#377eb8')
    ax3.set_xlim([0, len(nox) - 1])
    ax3.set_xticks([0, 30, 61])
    ax3.set_xticklabels(['1 June %d' %year, '1 July', '1 Aug'], fontsize = 12)
    for t in ax3.get_yticklabels():
        t.set_fontsize(12)
    ax3.set_ylabel('NO$_{x}$ [tons]', fontsize = 16)
    ax3.get_yaxis().set_label_coords(-0.07, 0.5)
    #
    pc = (nox - np.mean(nox))/np.mean(nox) * 100.
    ax3t = ax3.twinx()
    ax3t.plot(pc, lw = 1.0)
    for t in ax3t.get_yticklabels():
        t.set_fontsize(12)    
    ax3t.set_ylabel('Change [%]', rotation = 270, fontsize = 16)   
    ax3t.get_yaxis().set_label_coords(1.09, 0.5)
    ax3.xaxis.set_ticks_position('both')
    ax3.yaxis.set_ticks_position('both')
    plt.subplots_adjust(hspace = 0.3)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'timeseries_t2m_castneto3_cemsnox_%d_%s.eps' %(year, region), 
                dpi = 300)    
    print('T2m-O3 correlation = %.3f' %(np.corrcoef(t2m, castnet)[0, 1]))
    print('T2m-NOx correlation = %.3f' %(np.corrcoef(t2m, nox)[0, 1]))
    print('O3-NOx correlation = %.3f' %(np.corrcoef(castnet, nox)[0, 1])) 
    return 
# # # # # # # # # # # # #
def map_simulationschematic(lat_merra, lon_merra, T0, times, castnet_sites_fr,
                            years, sampling_months, sampling_hours):     
    """function plots a map of the region over which anthropogenic NO 
    emissions were varied. In this region the locations of industrial 
    facilities which report to CEMS were plotted, and their colors/scatterpoint
    sizes were determined by the JJA 2008-2010 cumulative NOx emissions. 
    The Northeastern U.S. is outlined, and two inset (zoomed) plots 
    illustrate "Transport" and "+ Emissions" simulations. n.b. The grid cells 
    referenced in the inset plots do NOT correspond to the data shown in the 
    insets. 
    
    Parameters
    ----------  
    lat_merra : numpy.ndarray
        MERRA-2 latitude coordinates of focus region, units of degrees north, 
        [lat,]
    lon_merra : numpy.ndarray
        MERRA-2 longitude coordinates of focus region, units of degrees east, 
        [lon,]
    T0 : numpy.ndarray 
        Surface air temperature, units of K, [time, lat, lon]
    times : numpy.ndarray
        Datetime objects corresponding to each 6-hourly timestep of MERRA-2,
        [time,]
    castnet_sites_fr : list
        CASTNET site names in focus region
    years : list
        Years of interest
    sampling_months : list 
        Months of interest
    sampling_hours : list 
        Hours during which trace gas concentrations are returned           

    Returns
    ----------
    None        
    """     
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
    from mpl_toolkits.axes_grid.inset_locator import mark_inset
    import matplotlib as mpl
    from matplotlib.lines import Line2D
    import shapely.geometry as sg
    import shapely.ops as so
    from descartes import PolygonPatch    
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    from geo_idx import geo_idx
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/emissions/')
    import cems_nei
    sys.path.append('/Users/ghkerr/phd/GMI/')
    import commensurability
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['xtick.minor.width'] = 1.5
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams['ytick.minor.width'] = 1.5
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
                colors.append('#fb9a99')                           
        return sizes, colors
    # # # # 
    # open 12Z gridded HindcastMR2 to get lat, lon coordinates
    (lat, lon, pressure, times_mr2, co, no, no2, o3) = \
    commensurability.open_gridded_idailyCTM('HindcastMR2', years)
    # open emissions inventories for HindcastMR2 run and HindcastMR2 run with 
    # temperature-dependent emissions run
    mr2_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
    commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2', 
        years, sampling_months, sampling_hours)
    mr2_no_emiss = \
    commensurability.commensurate_emiss_unperturbed(mr2_castnet, castnet_sites_fr, 
        years, sampling_months, '1x1.25_IAVanthGFED4')
    egu_no_emiss, emiss_lats_fr, emiss_lons_fr = \
    commensurability.commensurate_emiss_perturbed(mr2_castnet, castnet_sites_fr, 
        years, sampling_months)
    # convert longitude
    lon = np.mod(lon - 180.0, 360.0) - 180.0
    # initialize figure, axis
    fig = plt.figure(figsize = (8, 10))
    ax = plt.subplot2grid((1, 1), (0, 0))
    # focus region map 
    llcrnrlon = -93.
    llcrnrlat = 24.
    urcrnrlon = -66.3
    urcrnrlat = 48.
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
                llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, 
                urcrnrlat = urcrnrlat, resolution = 'h', area_thresh = 1000)
    m.drawmapboundary(color = '#888888', fill_color = '#dcf0fa')
    m.fillcontinents(color = '#f9f6d8', lake_color = '#dcf0fa')
    m.drawlsmask(ocean_color = '#dcf0fa')
    m.drawcoastlines(color = '#888888')
    # plot shapefile of NEUS
    m.readshapefile(pollutants_constants.PATH_SHAPEFILES + 
                    'cb_2015_us_state_20m', name = 'states', 
                    drawbounds = True, color = '#888888')
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
    patches_union = [] 
    for key, value in state_fips_code_listing.items():
        if key in pollutants_constants.NORTHEAST_STATES: 
                for info, shape in zip(m.states_info, m.states):
                    if info['NAME'] == key:
                        patches_union.append(sg.Polygon(shape))
    # cascaded union can work on a list of shapes, adapted from 
    # https://stackoverflow.com/questions/34475431/plot-unions-of-polygons-in-matplotlib
    neus = so.cascaded_union(patches_union) 
    ax.add_patch(PolygonPatch(neus, fill = False, ec = '#888888', 
                              zorder = 2, linewidth = 4.0)) 
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
              edgecolor = colors, zorder = 10)
    lonres = 1.25
    latres = 1.              
    # lon/lat index that will be highlighted
    latidx, lonidx = 5, 8
    # FIRST INSET, Transport simulation
    axins = zoomed_inset_axes(ax, 4.5, loc = 3, bbox_to_anchor = 
                              (-0.02, -0.02, 1,1), bbox_transform = ax.transAxes)
    x, y = m(lon[lonidx-4], lat[latidx-1])
    x2, y2 = m(lon[lonidx-4]+lonres, lat[latidx-1]+latres)
    axins.set_xlim(x, x2)
    axins.set_ylim(y, y2)
    axins.set_yticks([])
    axins.set_yticklabels([])
    mark_inset(ax, axins, loc1 = 2, loc2 = 4, fc = 'none', lw = 1.5)
    # find 6 hourly temperature at highlighted lon/lat index for a given month
    lat_where = geo_idx(lat[latidx], lat_merra)
    lon_where = geo_idx(lon[lonidx] + 360., lon_merra)
    times = pd.to_datetime(times)
    times_month = np.where((times.year == 2008) & (times.month == 6))[0]
    T0_gc = T0[times_month, lat_where, lon_where]
    # reshape to daily 
    T0_gc = np.reshape(T0_gc, (-1, 4)) # 4 for six hourly
    T0_gc = np.roll(T0_gc, 2, axis = 1)
    # plot in twinned inset axis
    axinst = axins.twinx()
    axinst.plot(np.linspace(x, x2, 4), T0_gc.T, color = 'grey', lw = 0.5)
    axinst.plot(np.linspace(x, x2, 4), np.mean(T0_gc.T, axis = 1), 
                color = 'k', lw = 2.5)
    axinst.set_xlim([x, x2])
    axinst.set_xticks(np.linspace(x, x2, 8))
    axinst.yaxis.tick_left()
    axins.set_xticklabels(['6', '', '12', '', '18', '', '0', ''])
    for t in axins.get_xticklabels():
        t.set_fontsize(12)
    axinst.text(0.96, -0.17, 'UTC', fontsize = 12, transform = axinst.transAxes)
    axinst.set_ylim([int(T0_gc.min()) - 2, int(T0_gc.max()) + 3])
    axinst.set_yticks(np.linspace(int(T0_gc.min()) - 2, 
                                  int(T0_gc.max()) + 3, 5))
    for t in axinst.get_yticklabels():
        t.set_fontsize(12)
    axins.set_ylabel('T$_{\mathregular{2\:m}}$ [K]', fontsize =16)
    axins.get_yaxis().set_label_coords(-0.35, 0.5)
    axins.set_xlabel('Transport, June 2008', fontsize = 16)
    # SECOND INSET, + Emissions simulation
    axins2 = zoomed_inset_axes(ax, 6, loc = 4, bbox_to_anchor = (0.02, -0.02, 1,1),
                               bbox_transform = ax.transAxes)
    x, y = m(lon[lonidx], lat[latidx])
    x2, y2 = m(lon[lonidx]+lonres, lat[latidx]+latres)
    axins2.set_xlim(x, x2)
    axins2.set_ylim(y, y2)
    axins2.set_yticks([])
    axins2.set_yticklabels([])
    mark_inset(ax, axins2, loc1 = 1, loc2 = 3, fc = 'none', lw = 1.5)
    axins2t = axins2.twinx()
    axins2t.plot(np.linspace(x, x2, 92), np.nanmean(egu_no_emiss, axis = 1)[0], 
                 linewidth = 1.5, label = 'Perturbed', color = '#984ea3', 
                 linestyle = ':')
    # plot unperturbed emissions (n.b. loop through indices of months)
    month_idx = [0, 30, 61, 92]
    bm_idx = list(np.linspace(x, x2, 4))
    i = 0
    for a, b in zip(month_idx[:-1], month_idx[1:]):
        axins2t.plot(np.linspace(bm_idx[i], bm_idx[i + 1], b-a), 
                     np.nanmean(mr2_no_emiss, axis = 1)[0, a:b], linewidth = 2.,  
                     color = '#984ea3')   
        i = i + 1
    axins2.set_xlim([x, x2])
    axins2.set_xticks(np.linspace(x, x2, 4))
    axins2.set_xticklabels(['June', 'July', 'Aug'], fontsize = 12)
    for t in axins2t.get_yticklabels():
        t.set_fontsize(12)  
    axins2t.set_ylabel('NO [kg s$^{-1}$ cell$^{-1}$]', fontsize = 16, 
                       rotation = 270)
    axins2t.get_yaxis().set_label_coords(1.45, 0.5)
    axins2.set_xlabel('+$\:$Emissions, 2008      ', 
                      ha = 'center', fontsize = 16)
    # legend for a scatter plot using a proxy artists 
    # see http://matplotlib.sourceforge.net/users/legend_guide.html#using-proxy-artist
    circle1 = Line2D(range(1), range(1), color = 'w', marker = 'o', 
                     markersize = np.log(18)**1.5, 
                     markerfacecolor = '#a6cee3', markeredgecolor = '#a6cee3',
                     label = '$\mathregular{\Sigma}$$\,$NO$_x$$\,$<$\,$25$^{\mathregular{th}}$')
    circle2 = Line2D(range(1), range(1), color = 'w', marker = 'o', 
                     markersize = np.log(26)**1.5, markerfacecolor = '#1f78b4', 
                     markeredgecolor = '#1f78b4', 
                     label = '25$^{\mathregular{th}}$$\,$$\mathregular{\leq}$$\,$$\mathregular{\Sigma}$$\,$NO$_x$$\,$<$\,$50$^{\mathregular{th}}$')
    circle3 = Line2D(range(1), range(1), color = 'w', marker = 'o', 
                     markersize = np.log(34)**1.5, 
                     markerfacecolor = '#b2df8a', markeredgecolor = '#b2df8a', 
                     label = '50$^{\mathregular{th}}$$\,$$\mathregular{\leq}$$\,$$\mathregular{\Sigma}$$\,$NO$_x$$\,$<$\,$75$^{\mathregular{th}}$')
    circle4 = Line2D(range(1), range(1), color = 'w', marker = 'o', 
                     markersize = np.log(42)**1.5, 
                     markerfacecolor = '#33a02c', markeredgecolor = '#33a02c', 
                     label = '75$^{\mathregular{th}}$$\,$$\mathregular{\leq}$$\,$$\mathregular{\Sigma}$ NO$_x$$\,$<$\,$90$^{\mathregular{th}}$')
    circle5 = Line2D(range(1), range(1), color = 'w', marker = 'o', 
                     markersize = np.log(50)**1.5, 
                     markerfacecolor = '#fb9a99', markeredgecolor = '#fb9a99', 
                     label = '$\mathregular{\Sigma}$$\,$NO$_x$$\,$$\mathregular{\geq}$$\,$90$^{\mathregular{th}}$')
    # add legend 
    leg = ax.legend(loc = 9, bbox_to_anchor = (0.5, 1.2),
                    handles = [circle1, circle2, circle3, circle4, circle5], 
                    ncol = 2, fontsize = 16, numpoints = 1, facecolor = 'w')
    leg.get_frame().set_linewidth(0.0)
    plt.subplots_adjust(right = 0.85)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 'map_simulationschematic.eps', 
                dpi = 350)
    return
# # # # # # # # # # # # #    
def contourf_merra2h500_mr2o3event(chemistry, H500, lon, lat, times, 
    sampling_hours):
    """function identifies days in which regionally-averaged O3 from Hindcast
    MR2 (+ Chemistry) simulation exceeds the 90th percentile for O3 for JJA
    2008-2010. MERRA-2 500 hPa geopotential heights on these days are aveaged 
    and shown as the average anomaly on event days for the first output figure
    and the geopotential height anomaly is shown for all individual O3 events 
    for the second output figure. 

    Parameters
    ----------      
    chemistry : numpy.ndarray            
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2, units of 
        ppbv, [years in measuring period, days in months in 'sampling_months']     
    H500 : numpy.ndarray
         MERRA-2 geopotential height at 500 hPa, units of m, [time, lat, lon]
    lon : numpy.ndarray
        Longitude coordinates of focus region, units of degrees east, [lon,]        
    lat : numpy.ndarray
        Latitude coordinates of focus region, units of degrees north, [lat,]
    times : numpy.ndarray
        Datetime objects corresponding to each 6-hourly timestep of MERRA-2,
        [time,]
    sampling_hours : list 
        Hours during which trace gas concentrations are returned        

    Returns
    ----------
    None       
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    chemistry = np.hstack(chemistry)
    # find indices of MERRA-2 meteorological fields in 'sampling_hours;' n.b.
    # this is meant to be for when 'sampling_hours' is 15-20 local time
    times = pd.to_datetime(times)
    times_intersect = np.in1d(times.hour, sampling_hours)
    H500_intersect = H500[times_intersect]
    # find O3 events in + Chemistry simulation (days where Northeast-averaged
    # O3 exceeds 80th percentile)
    events = np.where(chemistry > np.percentile(chemistry, 90))[0]
    # geopotential height on days with O3 events
    H500_events = H500_intersect[events]
    # colorscale and map 
    cmap = plt.get_cmap('Spectral_r')
    vmin, vmax = -150, 150
    clevs = np.linspace(vmin, vmax, 21)
    m = Basemap(projection = 'mill', llcrnrlon = -85., llcrnrlat = 35, 
                urcrnrlon = -66., urcrnrlat = 49., resolution = 'i', 
                area_thresh = 10000) 
    # convert longitude
    lon = np.mod(lon - 180.0, 360.0) - 180.0
    x, y = np.meshgrid(lon, lat)
    x, y = m(x, y)    
    # initialize figure for mean event
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    event = m.contourf(x, y, np.mean(H500_events, axis = 0) - 
            np.mean(H500_intersect, axis = 0), clevs, cmap = cmap, 
            vmin = vmin, vmax = vmax, extend = 'both')   
    m.drawstates(color = '#888888', linewidth = 0.5)
    m.drawcountries(color = '#888888', linewidth = 0.5)
    m.drawcoastlines(color = '#888888', linewidth = 0.5)
    #fig.subplots_adjust(right = 0.85)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7]) 
    fig.colorbar(event, cax = cbar_ax, label = '$\mathregular{\Delta}$Z [m]')   
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'contourf_merra2h500_mr2meano3event.eps', dpi = 300)
    # initialize figure for individual events 
    fig = plt.figure()
    for day in np.arange(0, len(events), 1):
        # n.b. will only work for 90th %-ile/28 events 
        ax = plt.subplot(4, 7, day + 1)
        ax.set_title('%s' %(str(times[np.where(times_intersect == 
            True)[0]][events][day])[:10]), fontsize = 8)
        event = m.contourf(x, y, H500_events[day] - 
            np.mean(H500_intersect, axis = 0), clevs, cmap = cmap, 
            vmin = vmin, vmax = vmax, extend = 'both')   
        m.drawstates(color = '#888888', linewidth = 0.5)
        m.drawcountries(color = '#888888', linewidth = 0.5)
        m.drawcoastlines(color = '#888888', linewidth = 0.5)
    # add central colorbar
    fig.subplots_adjust(left = 0.1, right = 0.80)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7]) 
    fig.colorbar(event, cax = cbar_ax, label = '$\mathregular{\Delta}$Z [m]')   
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'contourf_merra2h500_mr2allo3event.eps', dpi = 300)
    return    
# # # # # # # # # # # # #
def timeseries_t2m_datnox_dato3(transportno, transportno2, transporto3, 
    t2m, year, years):
    """function plots a timeseries of 2-meter temperatures from MERRA-2 and 
    NOx and O3 from the Transport simulation for a year of interest. 

    Parameters
    ----------  
    transportno : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    transportno2 : numpy.ndarray
        GMI CTM NO2 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    transporto3 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']            
    t2m : numpy.ndarray
        MERRA-2 2-meter temperatures co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of K, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']     
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
    # regionally-averaged 2-meter temperature
    t2m = np.nanmean(t2m, axis = 1)
    transporto3 = np.hstack(transporto3)
    transportno = np.hstack(transportno)
    transportno2 = np.hstack(transportno2)
    t2m = np.hstack(t2m)
    # initialize figure, axes
    fig = plt.figure()
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan = 3)
    ax2 = plt.subplot2grid((3, 3), (1, 0), colspan = 3) 
    ax3 = plt.subplot2grid((3, 3), (2, 0), colspan = 3)
    # plot trace gases and temperature for a particular summer
    where_year = np.where(np.array(years) == year)[0][0]
    # 2-meter temperatures
    ax1.plot(t2m[92*where_year:92*(where_year+1)], lw = 1.5, color = 'k')
    ax1.set_xlim([0, 91])
    ax1.set_xticks([0, 30, 61])
    ax1.set_xticklabels([''])
    for t in ax1.get_yticklabels():
        t.set_fontsize(12)  
    ax1.set_ylabel('T$_{\mathregular{2\:m}}$ [K]', fontsize = 16)
    ax1.get_yaxis().set_label_coords(-0.09, 0.5)
    # NOx
    ax2.plot(transportno[92*where_year:92*(where_year+1)] + 
                         transportno2[92*where_year:92*(where_year+1)], 
                         lw = 1.5, color = 'k')
    ax2.set_xlim([0, 91])
    ax2.set_xticks([0, 30, 61])
    ax2.set_xticklabels([''])
    for t in ax2.get_yticklabels():
        t.set_fontsize(12)  
    ax2.set_ylabel('NO$_{x}$ [ppbv]', fontsize = 16)
    ax2.get_yaxis().set_label_coords(-0.09, 0.5)
    # O3
    ax3.plot(transporto3[92*where_year:92*(where_year+1)], lw = 1.5, color = 'k')
    ax3.set_xlim([0, 91])
    ax3.set_xticks([0, 30, 61])
    ax3.set_xticklabels(['1 June', '1 July', '1 Aug %d' %year])
    for t in ax3.get_xticklabels():
        t.set_fontsize(12)  
    for t in ax3.get_yticklabels():
        t.set_fontsize(12)      
    ax3.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    ax3.get_yaxis().set_label_coords(-0.09, 0.5)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'timeseries_t2m_datnox_dato3_%d.eps' %year, dpi = 300)
    return
# # # # # # # # # # # # #
def scatter_datnox_dato3(transportno, transportno2, transporto3, t2m):
    """function plots regionally-averaged O3 versus NOx with the colorscale
    representing 2-meter temperatures.
    
    Parameters
    ----------  
    transportno : numpy.ndarray
        GMI CTM NO concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    transportno2 : numpy.ndarray
        GMI CTM NO2 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']    
    transporto3 : numpy.ndarray
        GMI CTM O3 concentrations co-located (or nearly colocated) with 
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
        units of ppbv, [years in measuring period, days in months in 
        'sampling_months']            
    t2m : numpy.ndarray
        MERRA-2 2-meter temperatures co-located (or nearly colocated) with 
        corresponding CASTNet stations, units of K, [years in measuring 
        period, stations in 'castnet_sites_fr', days in months in 
        'sampling_months']     

    Returns
    ----------
    None        
    """
    import numpy as np
    import matplotlib.pyplot as plt
    # regionally-averaged 2-meter temperature
    t2m = np.nanmean(t2m, axis = 1)
    transporto3 = np.hstack(transporto3)
    transportno = np.hstack(transportno)
    transportno2 = np.hstack(transportno2)
    t2m = np.hstack(t2m)
    # initialize figure, axis
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    mpble = ax.scatter(transportno + transportno2, transporto3, c = t2m, 
                       cmap = plt.get_cmap('gnuplot2'))
    plt.colorbar(mpble, label = 'T$_{\mathregular{2\:m}}$ [K]')
    ax.set_xlabel('NO$_{x}$ [ppbv]', fontsize = 16)
    ax.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize = 16)
    for t in ax.get_xticklabels():
        t.set_fontsize(12)
    for t in ax.get_yticklabels():
        t.set_fontsize(12)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'scatter_datnox_dato3.eps', dpi = 300)
# # # # # # # # # # # # #    
def map_do3dt2mratio_conus():
    """function opens GMI output for the continential U.S. (CONUS) and finds 
    the ratio between the slope of the linear fit through Transport O3 versus
    T2m and + Chemistry O3 versus T2m. These ratios are calculated for each 
    GMI site with hourly profile data rather than as regional-averages. Four 
    outplot maps include (1) O3-T2m correlation coefficients, (2) O3-T2m 
    sensitivity for the Transport simulation, (3) O3-T2m sensitivity for the 
    + Chemistry simulation, and (4) the aforementioned ratio for all GMI
    sites. 
    
    Parameters
    ----------       
    None
        
    Returns
    ----------      
    None             
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    # find ratio between Transport and + Chemistry O3 in each region 
    def dO3_dT(o3, t2m):
        """calculate slope through ozone versus 2-meter temperature, the ozone-
        temperature sensitivity at individual GMI/CASTNet sites. 
        
        Parameters
        ----------  
        o3 : numpy.ndarray
            O3 concentrations co-located (or nearly colocated) with 
            corresponding CASTNet stations, units of ppbv, [years in measuring 
            period, days in months in 'sampling_months']            
        t2m : numpy.ndarray        
            MERRA-2 2-meter temperatures co-located (or nearly colocated) with 
            corresponding CASTNet stations, units of K, [years in measuring 
            period, stations in 'castnet_sites_fr', days in months in 
            'sampling_months'] 
            
        Returns
        ----------               
        slopes : list
            The O3-T2m sensitivity for each GMI/CASTNet site
        """
        import numpy as np
        from scipy.stats import linregress
        slopes = []
        # swap axes so that dimensionality is (no. summers, no. days in summer, 
        # no. sites)
        o3 = np.swapaxes(o3, 1, 2)
        t2m = np.swapaxes(t2m, 1, 2)
        # stack year and day axes
        o3 = np.vstack(o3)
        t2m = np.vstack(t2m)
        # loop through sites and find dO3/dT
        for site in np.arange(0, t2m.shape[1], 1):
            slope = linregress(t2m[:, site], o3[:, site]).slope
            slopes.append(slope)
        return np.array(slopes)
    # # # # 
    def r_O3_dT(o3, t2m):
        """calculate Pearson correlation coefficient between ozone and 2-meter 
        temperature at individual GMI/CASTNet sites. 
       
        Parameters
        ----------  
        o3 : numpy.ndarray
            O3 concentrations co-located (or nearly colocated) with 
            corresponding CASTNet stations, units of ppbv, [years in measuring 
            period, days in months in 'sampling_months']            
        t2m : numpy.ndarray        
            MERRA-2 2-meter temperatures co-located (or nearly colocated) with 
            corresponding CASTNet stations, units of K, [years in measuring 
            period, stations in 'castnet_sites_fr', days in months in 
            'sampling_months'] 
            
        Returns
        ----------               
        correlations : list
            The O3-T2m correlation coefficient for each GMI/CASTNet site 
        """
        import numpy as np
        from scipy.stats import linregress
        correlations = []
        # swap axes so that dimensionality is (no. summers, no. days in summer, 
        # no. sites)
        o3 = np.swapaxes(o3, 1, 2)
        t2m = np.swapaxes(t2m, 1, 2)
        # stack year and day axes
        o3 = np.vstack(o3)
        t2m = np.vstack(t2m)
        # loop through sites and find dO3/dT
        for site in np.arange(0, t2m.shape[1], 1):
            correlation = linregress(t2m[:, site], o3[:, site]).rvalue
            correlations.append(correlation)
        return correlations
    # # # #    
    # CASTNet sites in continental U.S.
    castnet_sites_neus = ['ASH', 'WFM', 'WST', 'APT', 'SAR', 'CTH', 'WSP', 
       'ARE', 'BEL', 'PSU', 'LRL', 'PAR', 'PED', 'SHN', 'CDR', 'VPI', 'MKG', 'KEF']
    castnet_sites_south = ['ESP', 'SPD', 'GRS', 'BLR', 'COW', 'PNF', 'GRS', 'SND',
        'GAS', 'CVL', 'CAD', 'SUM']
    castnet_sites_midwest = ['PRK', 'BVL', 'VIN', 'ALH', 'OXF', 'SAL', 'LYK', 
        'DCP', 'UVL', 'ANA']
    castnet_sites_west = ['LAV', 'PIN', 'SEK', 'JOT', 'GLR', 'CHA', 'GRC', 'YEL', 
        'PND', 'CNT', 'ROM', 'GTH']
    # times to sample model
    sampling_hours_neus = [15, 16, 17, 18, 19, 20]
    sampling_hours_south = [15, 16, 17, 18, 19, 20]
    sampling_hours_midwest = [16, 17, 18, 19, 20, 21]
    sampling_hours_west = [17, 18, 19, 20, 21, 22]
    # open daily mean (11-16 LST average) CASTNet O3 and GMI trace gases
    # for the Northeast
    castnet_neus, mr2_o3_neus, mr2_no_neus, mr2_no2_neus, mr2_co_neus, \
    mr2_gmi_sites_neus = commensurability.commensurate_castnet_gmi(
        castnet_sites_neus, 'HindcastMR2', years, sampling_months, 
        sampling_hours_neus, 'GMT+4')
    temp, dat_o3_neus, dat_no_neus, dat_no2_neus, dat_co_neus, \
    dat_gmi_sites_neus = commensurability.commensurate_castnet_gmi(
        castnet_sites_neus, 'HindcastMR2-DiurnalAvgT', years, sampling_months, 
        sampling_hours_neus, 'GMT+4')
    # for the South
    castnet_south, mr2_o3_south, mr2_no_south, mr2_no2_south, mr2_co_south, \
    mr2_gmi_sites_south = commensurability.commensurate_castnet_gmi(
        castnet_sites_south, 'HindcastMR2', years, sampling_months, 
        sampling_hours_south, 'GMT+4')
    temp, dat_o3_south, dat_no_south, dat_no2_south, dat_co_south, \
    dat_gmi_sites_south = commensurability.commensurate_castnet_gmi(
        castnet_sites_south, 'HindcastMR2-DiurnalAvgT', years, sampling_months, 
        sampling_hours_south, 'GMT+4')
    # for the Midwest
    castnet_midwest, mr2_o3_midwest, mr2_no_midwest, mr2_no2_midwest, mr2_co_midwest, \
    mr2_gmi_sites_midwest = commensurability.commensurate_castnet_gmi(
        castnet_sites_midwest, 'HindcastMR2', years, sampling_months, 
        sampling_hours_midwest, 'GMT+5')
    temp, dat_o3_midwest, dat_no_midwest, dat_no2_midwest, dat_co_midwest, \
    dat_gmi_sites_midwest = commensurability.commensurate_castnet_gmi(
        castnet_sites_midwest, 'HindcastMR2-DiurnalAvgT', years, sampling_months, 
        sampling_hours_midwest, 'GMT+5')
    # for the West
    castnet_west, mr2_o3_west, mr2_no_west, mr2_no2_west, mr2_co_west, \
    mr2_gmi_sites_west = commensurability.commensurate_castnet_gmi(
        castnet_sites_west, 'HindcastMR2', years, sampling_months, 
        sampling_hours_west, 'GMT+6')
    temp, dat_o3_west, dat_no_west, dat_no2_west, dat_co_west, \
    dat_gmi_sites_west = commensurability.commensurate_castnet_gmi(
        castnet_sites_west, 'HindcastMR2-DiurnalAvgT', years, sampling_months, 
        sampling_hours_west, 'GMT+6')
    # look up CASTNet/GMI site locations in each region 
    castnet_lats_neus, castnet_lons_neus, gmi_lats_neus, gmi_lons_neus = \
    commensurability.commensurate_castnet_gmi_locations(castnet_sites_neus, 
        sampling_months, sampling_hours_neus, 'GMT+4')
    castnet_lats_south, castnet_lons_south, gmi_lats_south, gmi_lons_south = \
    commensurability.commensurate_castnet_gmi_locations(castnet_sites_south, 
        sampling_months, sampling_hours_south, 'GMT+4')
    castnet_lats_midwest, castnet_lons_midwest, gmi_lats_midwest, gmi_lons_midwest = \
    commensurability.commensurate_castnet_gmi_locations(castnet_sites_midwest, 
        sampling_months, sampling_hours_midwest, 'GMT+5')
    castnet_lats_west, castnet_lons_west, gmi_lats_west, gmi_lons_west = \
    commensurability.commensurate_castnet_gmi_locations(castnet_sites_west, 
        sampling_months, sampling_hours_west, 'GMT+6')
    # find 2-meter temperatures in each region
    (t2m_neus, t10m_neus, u2m_neus, u10m_neus, v2m_neus, v10m_neus, 
     merra_lats_neus, merra_lons_neus) = commensurability.commensurate_MERRA2(
        castnet_neus, castnet_sites_neus, years, sampling_months, 
        sampling_hours_neus)        
    (t2m_south, t10m_south, u2m_south, u10m_south, v2m_south, v10m_south, 
     merra_lats_south, merra_lons_south) = commensurability.commensurate_MERRA2(
        castnet_south, castnet_sites_south, years, sampling_months, 
        sampling_hours_south)  
    (t2m_midwest, t10m_midwest, u2m_midwest, u10m_midwest, v2m_midwest, v10m_midwest, 
     merra_lats_midwest, merra_lons_midwest) = commensurability.commensurate_MERRA2(
        castnet_midwest, castnet_sites_midwest, years, sampling_months, 
        sampling_hours_midwest)  
    (t2m_west, t10m_west, u2m_west, u10m_west, v2m_west, v10m_west, 
     merra_lats_west, merra_lons_west) = commensurability.commensurate_MERRA2(
        castnet_west, castnet_sites_west, years, sampling_months, 
        sampling_hours_west)  
    # find O3-T2m sensitivity 
    dat_slopes_neus = dO3_dT(dat_o3_neus, t2m_neus)
    mr2_slopes_neus = dO3_dT(mr2_o3_neus, t2m_neus)
    dat_slopes_south = dO3_dT(dat_o3_south, t2m_south)
    mr2_slopes_south = dO3_dT(mr2_o3_south, t2m_south)
    dat_slopes_midwest = dO3_dT(dat_o3_midwest, t2m_midwest)
    mr2_slopes_midwest = dO3_dT(mr2_o3_midwest, t2m_midwest)
    dat_slopes_west = dO3_dT(dat_o3_west, t2m_west)
    mr2_slopes_west = dO3_dT(mr2_o3_west, t2m_west)
    # find O3-T2m correlations
    mr2_correl_neus = r_O3_dT(mr2_o3_neus, t2m_neus)
    mr2_correl_south = r_O3_dT(mr2_o3_south, t2m_south)
    mr2_correl_midwest = r_O3_dT(mr2_o3_midwest, t2m_midwest)
    mr2_correl_west = r_O3_dT(mr2_o3_west, t2m_west)    
    # concatenate for CONUS
    all_lons = np.concatenate([gmi_lons_neus, gmi_lons_south, 
                               gmi_lons_midwest, gmi_lons_west])
    all_lats = np.concatenate([gmi_lats_neus, gmi_lats_south, 
                               gmi_lats_midwest, gmi_lats_west])
    slope_transport = np.concatenate([dat_slopes_neus, dat_slopes_south, 
                                      dat_slopes_midwest, dat_slopes_west])   
    slope_chemistry = np.concatenate([mr2_slopes_neus, mr2_slopes_south, 
                                      mr2_slopes_midwest, mr2_slopes_west])
    correl_chemistry = np.concatenate([mr2_correl_neus, mr2_correl_south, 
                                       mr2_correl_midwest, mr2_correl_west]) 
    # ratio of temperature sensitivities
    ratio_neus = (dat_slopes_neus/mr2_slopes_neus)
    ratio_south = (dat_slopes_south/mr2_slopes_south)
    ratio_midwest = (dat_slopes_midwest/mr2_slopes_midwest)
    ratio_west = (dat_slopes_west/mr2_slopes_west)
    all_ratio = np.concatenate([ratio_neus, ratio_south, ratio_midwest, ratio_west])
    # Eastern U.S. focus region map
    llcrnrlon = -94.
    llcrnrlat = 34.
    urcrnrlon = -66.3
    urcrnrlat = 50.
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
                llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, 
                urcrnrlat = urcrnrlat, resolution = 'h', area_thresh = 1000)
    x, y = m(all_lons, all_lats)
    fig = plt.figure(figsize = (8, 10))
    ax = plt.subplot2grid((1, 1), (0, 0))    
    m.drawmapboundary(color = '#888888', fill_color = '#dcf0fa')
    m.fillcontinents(color = '#f9f6d8', lake_color = '#dcf0fa')
    m.drawlsmask(ocean_color = '#dcf0fa')
    m.drawcoastlines(color = '#888888')
    m.drawstates(color = '#888888')
    m.drawcountries(color = '#888888')    
    cmap = plt.get_cmap('brg')        
    # find where O3-T2m correlation is less than threshold and set to nan.
    correl_chemistry_color = []
    for i in correl_chemistry:
        if i > 0.3:
            correl_chemistry_color.append(150)
        if i < 0.3:
            correl_chemistry_color.append(0)
    m.scatter(x, y, c = all_ratio[~np.isnan(all_ratio)], 
              s = (correl_chemistry_color), vmin = 0.3, vmax = 1.0, 
              cmap = cmap, zorder = 20)
    cb = m.colorbar(extend = 'both')
    cb.ax.tick_params(labelsize = 12)
    cb.set_label(label = 'Ratio O$_{3}$-T$_{\mathregular{2\:m}}$', size = 16)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_do3dt2mratio_conus_highcorrel.eps', dpi = 300)
    # Continental U.S. focus region map 
    llcrnrlon = -130.
    llcrnrlat = 24.8
    urcrnrlon = -66.3
    urcrnrlat = 50.
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
                llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, 
                urcrnrlat = urcrnrlat, resolution = 'h', area_thresh = 1000)
    x, y = m(all_lons, all_lats)
    # # # #
    # plot O3-T2m correlation
    fig = plt.figure(figsize = (8, 10))
    ax = plt.subplot2grid((1, 1), (0, 0))    
    m.drawmapboundary(color = '#888888', fill_color = '#dcf0fa')
    m.fillcontinents(color = '#f9f6d8', lake_color = '#dcf0fa')
    m.drawlsmask(ocean_color = '#dcf0fa')
    m.drawcoastlines(color = '#888888')
    m.drawstates(color = '#888888')
    m.drawcountries(color = '#888888')    
    cmap = plt.get_cmap('brg')    
    m.scatter(x, y, c = correl_chemistry[~np.isnan(correl_chemistry)], 
              vmin = -1., vmax = 1.0, cmap = cmap, zorder = 20)
    m.colorbar(extend = 'both', label = 'r')    
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_mr2o3t2mcorrelation_conus.eps', dpi = 300)
    # # # #
    # plot O3-T2m sensitivity in Transport simulation
    fig = plt.figure(figsize = (8, 10))
    ax = plt.subplot2grid((1, 1), (0, 0))    
    m.drawmapboundary(color = '#888888', fill_color = '#dcf0fa')
    m.fillcontinents(color = '#f9f6d8', lake_color = '#dcf0fa')
    m.drawlsmask(ocean_color = '#dcf0fa')
    m.drawcoastlines(color = '#888888')
    m.drawstates(color = '#888888')
    m.drawcountries(color = '#888888')    
    cmap = plt.get_cmap('brg')    
    m.scatter(x, y, c = slope_transport[~np.isnan(slope_transport)], 
              vmin = -0.25, vmax = 1.25, cmap = cmap, zorder = 20)
    m.colorbar(extend = 'both', 
               label = r'$\frac{\delta \mathregular{O}_' +
               '{\mathregular{3, Transport}}}' +
               '{\delta \mathregular{T}_{\mathregular{2\:m}}}$')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_ddato3dt2m_conus.eps', dpi = 300)    
    # # # # 
    # plot O3-T2m sensitivity in + Chemistry simulation 
    fig = plt.figure(figsize = (8, 10))
    ax = plt.subplot2grid((1, 1), (0, 0))    
    m.drawmapboundary(color = '#888888', fill_color = '#dcf0fa')
    m.fillcontinents(color = '#f9f6d8', lake_color = '#dcf0fa')
    m.drawlsmask(ocean_color = '#dcf0fa')
    m.drawcoastlines(color = '#888888')
    m.drawstates(color = '#888888')
    m.drawcountries(color = '#888888')    
    cmap = plt.get_cmap('brg')    
    m.scatter(x, y, c = slope_chemistry[~np.isnan(slope_chemistry)], 
              vmin = -0.25, vmax = 1.25, cmap = cmap, zorder = 20)
    m.colorbar(extend = 'both', 
               label = r'$\frac{\delta \mathregular{O}_' +
               '{\mathregular{3, +\:Chemistry}}}' +
               '{\delta \mathregular{T}_{\mathregular{2\:m}}}$')    
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_dmr2o3dt2m_conus.eps', dpi = 300)        
    # # # #
    # plot ratio of temperature sensitivities
    fig = plt.figure(figsize = (8, 10))
    ax = plt.subplot2grid((1, 1), (0, 0))
    m.drawmapboundary(color = '#888888', fill_color = '#dcf0fa')
    m.fillcontinents(color = '#f9f6d8', lake_color = '#dcf0fa')
    m.drawlsmask(ocean_color = '#dcf0fa')
    m.drawcoastlines(color = '#888888')
    m.drawstates(color = '#888888')
    m.drawcountries(color = '#888888')
    cmap = plt.get_cmap('brg')
    m.scatter(x, y, c = all_ratio[~np.isnan(all_ratio)], vmin = 0.5, 
              vmax = 1.5, cmap = cmap, zorder = 20)
    m.colorbar(extend = 'both', 
               label = 'm$_{\mathregular{Transport}}$/m$_{\mathregular{+\:Chemistry}}$')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_do3dt2mratio_conus.eps', dpi = 300)
    return
# # # # # # # # # # # # #    