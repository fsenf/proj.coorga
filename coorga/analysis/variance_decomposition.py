#!/usr/bin/env python


######################################################################
######################################################################


def var_decomp(hin, ndays = 31, nhours = 24):

    h = hin.reshape(ndays, nhours, -1)

    # daily average
    hd = h.mean(axis = 1)

    # mean dirunal cycle
    hc = h.mean(axis = 0)


    # residual time series
    hdrep = hd.reshape(ndays,1,-1).repeat(nhours, axis = 1)
    hcrep = hc.reshape(1, nhours, -1).repeat(ndays, axis = 0)

    hres = (h - hcrep - hdrep).reshape(ndays * nhours, -1)

    return hd.var(axis = 0), hc.var(axis = 0), hres.var(axis = 0)

######################################################################
######################################################################

