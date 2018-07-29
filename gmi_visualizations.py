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
    09042018 -- SOMETHING HERE                
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
    # set custom font
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
def cdf_castneto3gmio3(castnet, mr2_gmi, ccmi_gmi, ffigac2_gmi, 
                       ffigac2hr_gmi, years, region):
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
             1))(np.unique(t2m_ra)), color = '#377eb8')
    # EGU_T
    egu_m = np.poly1d(np.polyfit(t2m_ra, ccmi, 1))[1]
    axr.plot(np.unique(t2m_ra), np.poly1d(np.polyfit(t2m_ra, egu, 
             1))(np.unique(t2m_ra)), color = '#4daf4a')
    # HindcastMR2-DiurnalAvgT
    dat_m = np.poly1d(np.polyfit(t2m_ra, dat, 1))[1]
    axr.plot(np.unique(t2m_ra), np.poly1d(np.polyfit(t2m_ra, dat, 
             1))(np.unique(t2m_ra)), color = '#ff7f00')
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