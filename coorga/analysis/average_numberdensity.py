

##LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
## automatically generated by 18-Average_NARVAL-SMF_NumberDensites.ipynb
##TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


import numpy as np
import scipy.ndimage
import tropy.analysis_tools.grid_and_interpolation as gi

######################################################################
######################################################################

def calculate_average_numberdensity(d,
                                    smooth_sig = 2.,
                                    filter_in_logspace = False,
                                    var_names = None,
                                    bins = None,
                                    itime_range = None,
                                    dx = 20. ,
                                    dy = 20.,):

    
    '''
    Calculates average number density.
    

    Parameters
    ----------
    d : dict of arrays
       set of cell properties
    
    smooth_sig : float, optional, default = 2.
       sigma of Gaussian filter for smoothing the results

    filter_in_logspace : bool, optional, default = True,
       switch if filter is applied in log-space

    var_names : list of two strings, optional, default = None
       list of the two variables used for binning
       if None: 'x_mean' and 'y_mean' are used
       
    bins : list or tuple of two numpy arrays, optional, default = None
       sets the bins of histogram analysis
       if None: min and max is estimated from data and (dx, dy) is used
       
    itime_range : list or tuple of two int, optional, default = None
       selects a range of relative time for the analysis
       if None: all times are used
    
    dx : float, optional, default = 20. 
       interval of x-binning if bins == None

    dy : float, optional, default = 20.
       interval of y-binning if bins == None


    Returns
    --------
    egrid : tuple of two numpy arrays
       edge-based output grid

    nd : numpy array, 2dim
       number density field.


    Notes
    -----
    There might be some problems with the (x,y) versus the transformed (x,y)
    coordinates.

    ToDo:
    * implemnet super sampling
    
    '''
    
    # get the right variable names -----------------------------------
    if var_names is None:
        var_name1 = 'x_mean'
        var_name2 = 'y_mean'
    else: 
        var_name1, var_name2 = var_names

    x, y, tid, trel = d[var_name1], d[var_name2], d['time_id'], d['rel_time']
    # ================================================================

    
    # prepare time masking -------------------------------------------
    if itime_range is None:
        mask = np.ones_like( trel ).astype( np.bool )
    else:
        itime_min, itime_max = itime_range
        mask = (trel >= itime_min) & (trel <= itime_max)
        
    # ================================================================

                

    # transformed coordinate -> equator at zero
    # clon, clat = gi.xy2ll(1.*x, 1.*y)
    # xt, yt  =  gi.ll2xy(clon, clat, lon0 = 0, lat0 = 0)
    
    # mask out the spinup time
    # mask = trel > -12

    
    # get binning ----------------------------------------------------
    if bins is None:
        # take the values from the cell props
        xbins = np.arange(x.min(), x.max(), dx)
        ybins = np.arange(y.min(), y.max(), dy)

        bins = (xbins, ybins)
    else:
        xbins, ybins = bins
    # ================================================================


    # check for region -----------------------------------------------
    reg_xmask = ( x > xbins.min() ) & (x < xbins.max() )
    reg_ymask = ( y > ybins.min() ) & (y < ybins.max() )
    reg_mask = reg_xmask & reg_ymask
    mask = mask & reg_mask     
    # ================================================================


    # histogram analysis ---------------------------------------------
    h, xe, ye = np.histogram2d(x[mask], y[mask], bins)

    if filter_in_logspace:
        # smoothing ???
        logh = np.ma.log(h)
        logh.data[ logh.mask ] = -5.

        loghs = scipy.ndimage.gaussian_filter(logh, smooth_sig)
        h = np.exp(loghs)
    else:
        h = scipy.ndimage.gaussian_filter(h, smooth_sig)
    # ================================================================
    


    # grid center points ---------------------------------------------
    xm = gi.lmean(xe)
    ym = gi.lmean(ye)
    
    egrid = np.meshgrid(xe, ye, indexing = 'ij')
    # ================================================================


    
    # get the number of cells per time slot --------------------------
    target_times = sorted( set(tid[mask]) )
    ntimes = 1.* len(target_times)

    ncells =  1. * len (tid[mask] )
    ncells_per_slot = ncells /  ntimes
    # ================================================================


    # normalization --------------------------------------------------
    N0 = (h * dx * dy).sum()
    nd0 = ncells_per_slot * h / N0
    # ================================================================
    
    return egrid, nd0

######################################################################
######################################################################
