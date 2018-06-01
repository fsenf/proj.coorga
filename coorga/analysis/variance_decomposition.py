#!/usr/bin/env python


######################################################################
######################################################################


def var_decomp(hin, ndays = 31, nhours = 24):

    '''
    Variance Decomposition of a time series.

    
    Parameters
    ----------
    hin : numpy array with shape = (ntimes, nvars)
       input time series to decompose

    ndays : int, optional, default = 31
       number of days to be analyzed

    nhours : int, optional, default = 24
       number of hours that contributed to diurnal cycle, 
       e.g. 24 for hourly data, 12 for two-hourly data, etc


    Returns
    --------
    day_var : numpy array with shape = (nvars)
       daily-mean variance

    cyc_var : numpy array with shape = (nvars)
       average diurnal cycle

    res_var : numpy array with shape = (nvars)
       residual variance

    '''

    h = hin.reshape(ndays, nhours, -1)

    # daily average
    hd = h.mean(axis = 1)

    # residual time series
    hdrep = hd.reshape(ndays,1,-1).repeat(nhours, axis = 1)

    hd_residual = h - hdrep

    # mean dirunal cycle
    hc = hd_residual.mean(axis = 0)


    hcrep = hc.reshape(1, nhours, -1).repeat(ndays, axis = 0)

    hres = (h - hcrep - hdrep).reshape(ndays * nhours, -1)

    day_var, cyc_var, res_var = hd.var(axis = 0), hc.var(axis = 0), hres.var(axis = 0)

    return day_var, cyc_var, res_var

######################################################################
######################################################################

