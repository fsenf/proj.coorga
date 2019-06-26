#!/usr/bin/env python

import numpy as np
import xarray as xr

import tropy.analysis_tools.grid_and_interpolation as gi

######################################################################
######################################################################

def timevar_hist(d, vname, 
                 x = None, 
                 atot = 1.,
                 full_time_info = False,
                 nbins = 20,
                 logbinning = False, 
                 edge_values = False,
                 density = True,
                 do_cnsd = False, 
                 do_ccsd = False):

    

    '''
    Makes time-average 1d histograms for a certain variable.

    
    Parameters
    ----------
    d : dict
        set of cell properties

    vname : str
        variable name used for analysis

    x : numpy array, optional, default = None
        bins for histogram
        if None: either predefined bins or bins from percentiles are used
    
    atot : float, topional, default = 1.
        total area of analysis domain (for cnsd)

    full_time_info : bool, optional, default = False
        switch if full time information is returned (rather than ave and standard error)

    nbins : int, optional, default = 20
        number of bins

    logbinning : bool, optional, default = False
         switch if autogenerated bins use logarithmic bin distances

    edge_values : bool, optional, default = False
        if edge rather than mid-point values are returned
    
    density : bool, optional, default = True
        switch if PDF rather than absolute counts are returned

    do_cnsd : bool, optional, default = False
        if cloud number size distribution is calculated

    do_ccsd : bool, optional, default = False
        if cloud cover size distribution is calculated


    Returns
    --------
    xm : numpy array, 1d
        mid-points of bins
    
    pdf : numpy array, 2d 
        time-array of histograms per time slots

    hm : numpy array, 1d
        average histogram


    dhm : numpy array, 1d
        standard error of average histogram
    
    '''


    # set bins -------------------------------------------------------
    if np.any(x) == None:
        if vname in ['diameter', 'nN_diameter']:
            lx = np.linspace(np.log10(20), np.log10(1000), nbins)
            lx = np.linspace(np.log10(10), np.log10(1000), nbins)
            x = 10**lx

        elif vname == 'nN_distance':
            lx = np.linspace(np.log10(20), np.log10(1000), nbins)
            x = 10**lx

        elif vname == 'ml':
            lx = np.linspace(np.log10(1), np.log10(10000), nbins)
            x = 10**lx

        elif vname in ['shape_bt108', 'aspect']:
            x = np.linspace(0, 1, nbins)

        elif vname in ['bt108_min']:
            x = np.linspace(180, 230, nbins)
        else:
            # use percentile bins
            v1, v99 = np.percentile( d[vname], [1, 99])
            
            if logbinning: 
                lx = np.linspace(np.log10(v1), np.log10(v99), nbins)
                x = 10**lx
            else:
                x = np.linspace(v1, v99, nbins)
        

    if vname in ['diameter', 'nN_diameter', 'nN_distance'] and np.any(x) == None:
        lx = np.log10( x )
        lxm = gi.lmean(lx)
        xm = 10**lxm
    else:
        xm = gi.lmean( x )
    # ================================================================


    # collect histograms ---------------------------------------------
    tlist = sorted( set(d['time_id']) )

    ntime = len(tlist)
    nbins = len(xm)

    pdf  = np.zeros((ntime, nbins))
    cnsd = np.zeros((ntime, nbins))
    ccsd = np.zeros((ntime, nbins))

    range_mask = ( d[vname] >= x.min() ) & ( d[vname] <= x.max() )

    for itime, t in enumerate( tlist ):
        

        m = ( d['time_id'] == t )  & range_mask

        v = d[vname][m]

        Ntot = len( v )

            
        pdf[itime], x = np.histogram(v, x, density = density)
        # pdf[itime] = scipy.ndimage.gaussian_filter(pdf[itime], 0) #.5)


        cnsd[itime] = Ntot * pdf[itime] / atot
        ccsd[itime] = 100. * np.pi / (4.) * ( cnsd[itime] * xm**3)
    # ================================================================


    pdf = np.ma.masked_invalid(pdf)
    cnsd = np.ma.masked_invalid(cnsd)
    ccsd = np.ma.masked_invalid(ccsd)


    if do_cnsd:
        pdf = cnsd
    elif do_ccsd:
        pdf = ccsd

    hm = pdf.mean(axis = 0)
    hsig = pdf.std(axis = 0)

    # try an standard error approach
    ntimes, nbins = pdf.shape
    dhm = 2 * hsig / np.sqrt( ntimes )


    if edge_values:
        xout = x
    else:
        xout = xm

    if full_time_info:
        return xout, pdf
    else:
        return xout, hm, dhm

######################################################################
######################################################################

def histxr(dset, varname, bins, loop_dim = 'time', do_log = False):
    
    '''
    Calculates histogram using xarray data structure. Allows to loop one dimension (e.g. time).
    
    
    Parameters
    ----------
    dset : xarray Dataset
        input data container
        
    varname : str
        variable name for which binning is done
        
    bins : numpy array or list
        set of bin edges for histogram computations
        
    loop_dim : str, optional
        dimension name for which looping is done
        default = 'time'
        
    do_log : bool, optional
        switch for bin center calculations
        if True: bins are assume to be equi-distant in log-space
        
        
    Returns
    --------
    hset : xarray Dataset
        histogram data, incl. pdf, cdf and cumulative moments
    '''
    
    # ----------------------------
    # (1) Preparation Part
    # ----------------------------
    
    # nbin center points
    if do_log:
        logbins = np.log( bins )
        clogbins = gi.lmean( logbins )
        cbins = np.exp( clogbins )
    else:
        cbins = gi.lmean( bins )

        
    # ----------------------------
    # (2) Loop Part
    # ----------------------------
    nloop = len( dset.coords[ loop_dim ] )
    hist_list = []
    
    for i in range( nloop ):
        d = dset.isel( {loop_dim : i} )[varname].data
        
        h, xe = np.histogram( d.flatten(), bins)
        
        hist_list += [ h.copy(), ]
    
    harray = np.row_stack( hist_list )
    
    
    # ----------------------------
    # (3) Xarray Part
    # ----------------------------
    # get variable attributes
    attrs = dict( dset[varname].attrs )
    attrs['varname'] = varname
    attrs['var_longname'] = attrs.pop('longname', 'None')
    
    
    # store absolute counts
    hset = xr.Dataset({'histogram': ([loop_dim, 'cbin'],  harray, {'longname':'absolute histogram counts'}) }, 
                    coords={loop_dim: (loop_dim, dset[loop_dim]), 
                            'ebin': ('ebin', bins, dict(longname =  'bin egdes', **attrs) ),
                            'cbin': ('cbin', cbins, dict(longname = 'bin mid-points', **attrs) ),})
    
    hset['delta_bin'] = xr.DataArray( hset.ebin.diff('ebin').data, dims = 'cbin', 
                                     attrs=dict(longname = 'bin widths', **attrs) )
    
    
    # pdf
    h_relative = hset.histogram / (1. * hset.histogram.sum( dim = 'cbin'))
    pdf  =  h_relative / hset.delta_bin
    hset['pdf'] = xr.DataArray( pdf, attrs = {'longname':'probability density function'})
    
    # cdf
    cdf = ( pdf * hset.delta_bin ).cumsum( dim = 'cbin' )
    hset['cdf'] = xr.DataArray( cdf, attrs = {'longname':'cumulative distribution function'})

    # cumulative moments
    for i in range(1,3):
        mom = hset.cbin**i
        cum_mom = ( mom * pdf * hset.delta_bin ).cumsum( dim = 'cbin' )
        cum_mom_name = 'cum_mom%d' % i
        hset[cum_mom_name] = xr.DataArray( cum_mom, 
                                           attrs = {'longname':'cumulative %d. moment' % i})
    
    
    return hset

