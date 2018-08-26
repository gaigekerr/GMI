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
                'diurnal_aqscogmico' added
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
                      llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, years, 
                      sampling_months, sampling_hours, region):        
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
    llcrnrlon : float
        Lower left longitude value of map, units of degrees east
    llcrnrlat : float
        Lower left latitude value of map, units of degrees north    
    urcrnrlon : float
        Upper right longitude value of map, units of degrees east    
    urcrnrlat : float
        Upper right latitude value of map, units of degrees north
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
    import matplotlib.font_manager as font_manager 
    from matplotlib.patches import Polygon 
    import matplotlib as mpl    
    import matplotlib.patches as mpatches
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
    sys.path.append('/Users/ghkerr/phd/stagnation/')
    from plot_focus_region import plot_focus_regionOLD
    sys.path.append('/Users/ghkerr/phd/GMI/')
    import commensurability
    # set custom font
    path = pollutants_constants.FONTPATH_LIGHT
    prop = font_manager.FontProperties(fname=path)
    mpl.rcParams['font.family'] = prop.get_name()
    path = pollutants_constants.FONTPATH_BOLD
    prop = font_manager.FontProperties(fname = path)
    mpl.rcParams['mathtext.bf'] = prop.get_name()
    # open CASTNet and GMI sites
    castnet_lats, castnet_lons, gmi_lats, gmi_lons = \
    commensurability.commensurate_castnet_gmi_locations(castnet_sites_fr, 
        sampling_months, sampling_hours)
    # open MERRA (requires a year's worth of GMI CTM output)
    mr2_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr = \
    commensurability.commensurate_castnet_gmi(castnet_sites_fr, 'HindcastMR2', 
                                              years, sampling_months, 
                                              sampling_hours)
    comm_t2m, merra_lats, merra_lons = commensurability.commensurate_t2m(
            mr2_castnet, castnet_sites_fr, years, sampling_months)
    del mr2_castnet, mr2_o3, mr2_no, mr2_no2, mr2_co, mr2_gmi_sites_fr    
    # initialize figure, axis, and a new instance of basemap
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    m = Basemap(projection = 'merc', llcrnrlon = llcrnrlon, 
                llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, 
                urcrnrlat = urcrnrlat, resolution = 'c', area_thresh = 10000)
    # plot MERRA grid cells
    x_merra, y_merra = m(merra_lons, merra_lats)
    for xi, yi in zip(x_merra, y_merra):
        # find longitude, latitude of every unique grid cell
        xinverse, yinverse = m(xi, yi, inverse = True)
        # change bounds of grid cell to map projection (MERRA longitude as 
        # 2/3 degree resolution, latitude has 1/2 degree resolution)
        x1, y1 = m(xinverse + 1/3., yinverse - 1/4.)
        x2, y2 = m(xinverse + 1/3., yinverse + 1/4.)
        x3, y3 = m(xinverse - 1/3., yinverse + 1/4.)
        x4, y4 = m(xinverse - 1/3., yinverse - 1/4.)
        # add patch to map
        p = Polygon([(x1, y1),(x2, y2),(x3,y3),(x4, y4)], facecolor = 
                     '#fb9a99', alpha = 1., edgecolor = 'none', 
                     linewidth = 0.0, zorder = 20)
        plt.gca().add_patch(p) 
    # plot GMI sites 
    x_gmi, y_gmi = m(gmi_lons, gmi_lats)
    gmi_loc = m.scatter(x_gmi, y_gmi, 25, color = '#1f78b4', marker = 'x', 
                        edgecolor = 'none', linewidth = 0.85, zorder = 24, 
                        label = 'GMI CTM')
    # plot CASTNet stations
    x_castnet, y_castnet = m(castnet_lons, castnet_lats)
    castnet_loc = m.scatter(x_castnet, y_castnet, 20, color = 'k', marker = '^', 
                           edgecolor = 'k', linewidth = 0.5, zorder = 23, 
                           label = 'CASTNet')
    # plot AQS stations measuring NO2
    aqs_no2 = np.array([np.vstack(aqs_no2_coords)[:, 1], 
                        np.vstack(aqs_no2_coords)[:, 0]], dtype = float).T
    aqs_no2 = pd.DataFrame(aqs_no2)
    aqs_no2 = aqs_no2.drop_duplicates()
    x_aqs_no2, y_aqs_no2 = m(aqs_no2[0].values, 
                             aqs_no2[1].values)
    aqs_no2_loc = m.scatter(x_aqs_no2, y_aqs_no2, 8, color = '#33a02c', 
                            marker = 'o', edgecolor = 'none', linewidth = 0.5, 
                            zorder = 21, label = 'AQS NO$_{2}$')
    # plot AQS stations measuring CO
    aqs_co = np.array([np.vstack(aqs_co_coords)[:, 1], 
                       np.vstack(aqs_co_coords)[:, 0]], dtype = float).T
    aqs_co = pd.DataFrame(aqs_co)
    aqs_co = aqs_co.drop_duplicates()
    x_aqs_co, y_aqs_co = m(aqs_co[0].values, aqs_co[1].values)
    aqs_co_loc = m.scatter(x_aqs_co, y_aqs_co, 8, color = '#b2df8a', marker = 'd', 
                           edgecolor = 'none', linewidth = 0.5, zorder = 22, 
                           label = 'AQS CO')
    # overlay shapefile of Northeast
    plot_focus_regionOLD(ax, m)
    plt.tight_layout()
    # create patch which represents MERRA
    merra_patch = mpatches.Patch(facecolor = '#fb9a99', alpha = 1., edgecolor = 'none', 
                     linewidth = 0.0, label = 'MERRA')
    leg = ax.legend(handles = [gmi_loc, castnet_loc, aqs_no2_loc, 
                    aqs_co_loc, merra_patch], bbox_to_anchor = (0.8, 0.5), 
                    ncol = 1, fontsize = 12, scatterpoints = 1)
    leg.get_frame().set_linewidth(0.0)
    # remove axis frame from plot             
    ax.axis('off')
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'map_gmiaqscastnet_%s.eps' %(region), dpi = 300)
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
        corresponding CASTNet stations for model case HindcastMR2-DiurnalAvgT, 
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
def timeseries_mr2o3dato3eguto3(transport, chemistry, emissions, obs, region):
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
    f, (ax, ax2, ax3) = plt.subplots(1, 3, sharey = True, figsize=(10, 6))
    for axt in [ax, ax2, ax3]:    
        axt.plot(obs, '-', color = 'lightgrey', zorder = 1, lw = 2,
                 label = 'CASTNet')
        axt.plot(emissions, '-', color = pollutants_constants.COLOR_EMISSIONS, 
                 label = '+$\:$Emissions')
        axt.plot(chemistry, '-', color = pollutants_constants.COLOR_CHEMISTRY, 
                 label = '+$\:$Chemistry')
        axt.plot(transport, '-', color = pollutants_constants.COLOR_TRANSPORT, 
                 label = 'Transport')
    ax.set_xlim([0, 92 + 1])
    ax2.set_xlim([92, 92*2 + 1])
    ax3.set_xlim([92*2, 92*3])
    # hide the spines between axes
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax2.tick_params(axis = 'y', which = 'both', left = False, right = False)
    ax3.yaxis.tick_right()
    # add discontinuous lines to axis (diagonal lines)
    d = .015
    kwargs = dict(transform = ax.transAxes, color='k', lw = 0.75, clip_on=False)
    ax.plot((1-d,1+d), (-d,+d), **kwargs)
    ax.plot((1-d,1+d),(1-d,1+d), **kwargs)
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d,+d), (1-d,1+d), **kwargs)
    ax2.plot((-d,+d), (-d,+d), **kwargs)
    ax2.plot((1-d,1+d), (-d,+d), **kwargs)
    ax2.plot((1-d,1+d),(1-d,1+d), **kwargs)
    kwargs.update(transform=ax3.transAxes)
    ax3.plot((-d,+d), (1-d,1+d), **kwargs)
    ax3.plot((-d,+d), (-d,+d), **kwargs)
    # move axes closer to each other
    f.subplots_adjust(wspace = 0.05)
    # axis labels
    ax.set_xticks([0, 30, 61])
    ax.set_xticklabels(['Jun', 'Jul', 'Aug'], ha = 'left', fontsize = 12)
    ax.set_xlabel('2008', x = 0., ha = 'left', fontsize = 14)
    ax2.set_xticks([92, 122, 153])
    ax2.set_xticklabels(['Jun', 'Jul', 'Aug'], ha = 'left', fontsize = 12)
    ax2.set_xlabel('2009', x = 0., ha = 'left', fontsize = 14)
    ax3.set_xticks([184, 214, 245])
    ax3.set_xticklabels(['Jun', 'Jul', 'Aug'], ha = 'left', fontsize = 12)
    ax3.set_xlabel('2010', x = 0., ha = 'left', fontsize = 14)
    ax.set_ylabel('Ozone [ppbv]', fontsize = 14)
    for tl in ax.get_yticklabels():
        tl.set_fontsize(12)
    # remove axis tick labels overlap discontinuous lines
    for axt in [ax2, ax3]:
        i = 0
        for t in axt.xaxis.get_ticklines(): 
            if i == 0:
                t.set_visible(False) 
            i = i + 1
    # reorder and place legend
    handles, labels = ax.get_legend_handles_labels()
    handles = [handles[0], handles[3], handles[2], handles[1]]
    labels = [labels[0], labels[3], labels[2], labels[1]]
    ax.legend(handles,labels, bbox_to_anchor = (1.55, 1.05), loc = 'center', 
              ncol = 4, frameon = False, fontsize = 12)
    plt.savefig('/Users/ghkerr/phd/GMI/figs/' + 
                'timeseries_mr2o3dato3eguto3_%s.eps' %(region), 
                dpi = 300)   
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
    from matplotlib.ticker import NullFormatter
    import sys
    sys.path.append('/Users/ghkerr/phd/')
    import pollutants_constants
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
    # scatter plot
    axScatter.scatter(t2m, obs, label = 'CASTNet', s = 20, 
                      color = 'lightgrey', zorder = 2)    
    axScatter.scatter(t2m, transport, label = 'Transport', s = 20, 
                      color = pollutants_constants.COLOR_TRANSPORT, zorder = 9)
    axScatter.scatter(t2m, chemistry, label = '+$\:$Chemistry', s = 20, 
                      color = pollutants_constants.COLOR_CHEMISTRY, zorder = 6)
    axScatter.scatter(t2m, emissions, label = '+$\:$Emissions', s = 20, 
                      color = pollutants_constants.COLOR_EMISSIONS, zorder = 3)
    axScatter.set_xlabel('T$_{\mathregular{2\:m}}$ [K]', fontsize = 14)
    axScatter.set_ylabel('Ozone [ppbv]', fontsize = 14)
    # histograms
    nbins = 18
    n, x, _  = axHistx.hist(t2m, bins = nbins, histtype = u'step', 
                            density = True, color = 'k', lw = 0.)
    density = stats.gaussian_kde(t2m)
    axHistx.plot(x, density(x), lw = 2., color = 'k')
    axHistx.set_ylabel('Density', fontsize = 14)
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
    axHisty.set_xlabel('Density', fontsize = 14)
    axHistx.set_xlim(axScatter.get_xlim())
    axHistx.spines['top'].set_visible(False)
    axHistx.spines['right'].set_visible(False)
    axHisty.set_ylim(axScatter.get_ylim())
    axHisty.spines['top'].set_visible(False)
    axHisty.spines['right'].set_visible(False)
    axScatter.legend(ncol = 1, frameon = False, fontsize = 12)
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