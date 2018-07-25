#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAME
    generate_times.py
PURPOSE
    Generate a list of datetime objects corresponding to CTM output
PROGRAMMER
    Gaige Hunter Kerr
REVISION HISTORY
    20072018 -- initial version created
"""
# # # # # # # # # # # # #
def generate_times(year, smonth, emonth, step):
    """function generates a list of datetime objects representing each hourly 
    model timestep in measuring period
  
    Parameters
    ----------   
    year : int
        Year of interest
    smonth : int
        Starting month, times will be retrieved starting at the first day 
        of the month 
    emonth : int
        Ending month, times will be retrieved until the end of the month
    step : int
        Timestep, units of hours
    
    Returns
    ----------
    times : list 
        Datetime objects for each model timestep, [step * no. days in months,]
    """
    from datetime import datetime, timedelta
    # find indicies of relevant dates/hours
    def perdelta(start, end, delta):
        curr = start
        while curr < end:
            yield curr
            curr += delta
    times = []
    # loop through days in year of interest with an hourly timestep
    for result in perdelta(datetime(year, smonth, 1, 0), 
                           datetime(year, emonth + 1, 1, 0), 
                           timedelta(hours = step)):
        times.append(result)     
    return times
# # # # # # # # # # # # #