def geo_idx(dd, dd_array):
    """Function searches for nearest decimal degree in an array of decimal 
    degrees and returns the index. np.argmin returns the indices of minimum 
    value along an axis. So subtract dd from all values in dd_array, take 
    absolute value and find index of minimum.
    
    Parameters
    ----------
    dd : int
        Latitude or longitude whose index in dd_array is being sought
    dd_array : numpy.ndarray 
        1D array of latitude or longitude 
    
    Returns
    -------
    geo_idx : int
        Index of latitude or longitude in dd_array that is closest in value to 
        dd
    """
    import numpy as np    
    geo_idx = (np.abs(dd_array - dd)).argmin()
    # if distance from closest cell to intended value is 2x the value of the
    # spatial resolution, raise error 
    res = np.diff(dd_array).mean()
    if np.abs(dd_array[geo_idx] - dd) > (2 * res):
        print('Closet index far from intended value!')
        return 
    return geo_idx