#!/usr/bin/env python

import numpy as np


######################################################################
######################################################################

def binwise_average(d, 
                    variable_name, 
                    selector_name, 
                    selector_bins, 
                    operator = np.ma.mean,
                    output_tseries = False):
    
    '''
    Calculates Average of one variable conditioned on the range of a second 
    variable.
    

    Parameters
    ----------
    d : dict
        set of cell properties

    variable_name : str
        variable name for which average is calculated

    selector_name : str
        variable name which is selected for conditioning

    selector_bins : list or tuple
        bin values for conditioning

    operator : func, optional, default = np.ma.mean
        function applied for statistics (e.g. averaging)

    output_tseries : bool, optional, default = False
        switch if ave and sem, OR, full time series is returned


    Returns
    -------
    ave : numpy array, 1d
        temporally averaged statistics

    sem : numpy array, 1d
        standard error of ave

    stats : numpy array, 2d, optional if output_tseries == True
        time series of binwise-statistics


    Notes
    -----
    First averages per configuration (i.e. per time slot) are calculated. Then, these
    averages are averaged by time and a standard error is derived from the variation in time.
    '''

    # extract variable and selector
    v = d[variable_name]
    s = d[selector_name]


    
    # get number of bins
    nbins = len( selector_bins )

    
    # get time array
    tlist = sorted( set(d['time_id']) )
    ntime = len(tlist)

    stats  = np.zeros((ntime, nbins - 1))

    for itime, t in enumerate( tlist ):
        
        # time mask
        m = d['time_id'] == t
        
        vm = v[m]
        sm = s[m]

        for i in range(nbins - 1):
                
            # set bin edges 
            s1 = selector_bins[i]
            s2 = selector_bins[i + 1]
    
            # prepare mask
            m = (sm >= s1) & (sm < s2)

            # do the averaging
            stats[itime, i] = operator( vm[m] )
            

    ave = stats.mean( axis = 0 )
    std = stats.std( axis = 0 )
    sem = 2 * std / np.sqrt( ntime )

    if output_tseries: 
        return stats
    else:
        return ave, sem

    
######################################################################
######################################################################

