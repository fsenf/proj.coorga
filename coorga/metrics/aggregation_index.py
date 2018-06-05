#!/usr/bin/env python

# load libraries -----------------------------------------------------
import sys, os, glob

import numpy as np
import scipy.stats, scipy.ndimage

import datetime


import tropy.analysis_tools.grid_and_interpolation as gi
import tropy.analysis_tools.segmentation as seg
import tropy.io_tools.hdf as hio

from standard_config import *


######################################################################
######################################################################

def cluster_analysis( lon, lat, c, dmin = 20., 
                      do_fast = False, xya = None):

    '''
    Simple version of cluster analysis used for aggregation metrics.
    
    Parameters
    -----------
    lon : numpy array, 2dim with shape = (nrows, ncols) 
        longitude

    lat :  numpy array, 2dim with shape = (nrows, ncols) 
        latitude

    c : numpy array, 2dim with shape = (nrows, ncols), type = int
        categorical cluster field
    
    dmin : int or float, optional, default = 20
        minimum diameter kept in the cell set data

    do_fast : bool, optional, default = False
        switch to use faster calculation of cell properties
    
    xya : tuble of 3 numpy arrays, optional, default = None
        recalculated fields of x- and y-coordinate and gridbox area
        (x, y, a) = xya


    Returns
    --------
    ca : dict
        cluster analysis dict
    '''
    
    
    # get geo and area fields ----------------------------------------
    if xya == None:
        x, y = gi.ll2xy(lon, lat)
        a = np.abs( gi.simple_pixel_area(lon , lat) )
    else:
        x, y, a = xya
    # ================================================================
        

    # init output dict -----------------------------------------------
    ca = {}


    if do_fast:
        index = range(1, c.max() + 1)
        ca['xc'] = scipy.ndimage.measurements.mean(x, c, index = index)
        ca['yc'] = scipy.ndimage.measurements.mean(y, c, index = index)
        ca['area'] = scipy.ndimage.measurements.sum(a, c, index = index)
    else:
        vnames = ['area',  'xc', 'yc']
        for vname in vnames:
            ca[vname] = []
     
     
        # loop over clusters---------------------------------------------
     
        for n in range(1, c.max() + 1):
                
            # make cluster mask
            m = (c == n)
     
            ca['area'].append( a[m].sum() )
            ca['xc'].append(  x[m].mean() )
            ca['yc'].append(  y[m].mean() )
     
        for vname in ca.keys():
            ca[vname] = np.array( ca[vname] )
        # ================================================================


    # calculate diameter ---------------------------------------------
    ca['dia'] = 2 * np.sqrt( ca['area'] / np.pi)
    # ================================================================


    # do masking -----------------------------------------------------
    m = ca['dia'] > dmin
    for vname in ca.keys():
        ca[vname] = ca[vname][m]
    # ================================================================


    # get number of clusters -----------------------------------------
    ca['number'] = len( ca['dia'] )
    # ================================================================


    # calculate distances field --------------------------------------
    
    # make a mesh for fast matrix-based distance calculation
    xx, yy = np.meshgrid( ca['xc'], ca['yc'])

    # get squared direction deviations
    dxq = (xx - xx.transpose())**2
    dyq = (yy - yy.transpose())**2


    dist_matrix = np.sqrt( dxq + dyq ) 
    
    # take upper triangular matrix
    dist_matrix = np.triu( dist_matrix )
    mask_matrix = (dist_matrix != 0)
    
    
    ca['dist'] = dist_matrix[ mask_matrix ]
    # ================================================================

    
    # calculate diameters --------------------------------------------
    ca['D0'] = scipy.stats.gmean( ca['dist'] )
    ca['D1'] = np.mean( ca['dist'] )
    # ================================================================

    return ca

######################################################################
######################################################################


def tobin_scai(ca, Nmax = 10, L = 10):

    '''
    Calculates the Tobin et al. (2012) SCAI index.


    Parameters
    -----------
    ca : dict
        cluster analysis field

    Nmax : int, optional, default = 10
        number of all grid points

    L : int, optional, default = 10
        geometrical size of the domain
    

    Returns
    --------
    scai : float
        SCAI index

    

    See Also
    --------
    Tobin, I., S. Bony, and R. Roca (2012), Observational Evidence for Relationships between the Degree of Aggregation of Deep Convection, Water Vapor, Surface Fluxes, and Radiation, J. Climate, 25(20), 6885-6904, doi:10.1175/jcli-d-11-00258.1.
    
    '''

    
    scai_max = 1. * Nmax * L

    scai = ( ca['number'] * ca['D0'] ) / scai_max * 1e3

    

    return scai

######################################################################
######################################################################


def get_domain_properties(lon, lat):

    '''
    Calculates properties of the considered domain, typical domain size and number of pixels.


    
    Parameters
    ----------
    lon : numpy array
       longitude

    lat : numpy array
       latitude
    
    
    Returns
    --------
    Nmax : int
       number of pixels / grid points in the domain
    
    L : float
       typical size / diameter of the domain (square-equivalent)
    '''
    


    # get projected coordinates ......................................
    x, y = gi.ll2xy(lon, lat)


    # construct a rectangle using maximum dimensions .................
    dx = x.max() - x.min()
    dy = y.max() - y.min()

    A = dx * dy


    # and then calculate size of an equal-area suqare .................
    L = np.sqrt(A)


    # number of pixels ...............................................
    Nmax =  x.size

    return Nmax, L


######################################################################
######################################################################

def coverage(f, thresh1, thresh2, dthresh, inverse = False):

    '''
    Calculates average fraction of pixels / grid points above or below
    a certain threshold. 

    A threshold range with step interval can be given and calculations 
    are subsequently performed for all different thresholds.


    
    Parameters
    ----------
    f : numpy array
        2d or 3d field

    thresh1 : float
        lower threshold

    thresh2 : float 
        upper threshold

    dthresh : float
        interval to increase threshold 

    inverse : bool, optional, default = False
        switch to use values below threshold (if True)

    
    Returns
    --------
    cov: dict
        dictionary of coverage values depending on threshold
    '''
 

    Nmax = np.size(f)

    tlist = np.arange(thresh1, thresh2, dthresh)

    cov = {}
    for t in tlist:

        # minimum or maximum ?
        if inverse:
            m = f < t
        else:
            m = f > t

        # get number of pixel above (below) threshold
        npix = len( f[m] )

        cov[t] = (100. * npix) / Nmax


    return cov

######################################################################
######################################################################

    

if __name__ == "__main__":


    # ----------------------------------------------------------------
    date_stamp = sys.argv[1]
    # ================================================================


    # input data -----------------------------------------------------
    ad = '%s/SEVIRI/tseries' % local_data_path
    region_flag = 'de'    
    fname = '%s/tseries-%s-%s-rss-%s.h5' % (ad, region_flag, '*', date_stamp)
    fname = glob.glob(fname)[0]
    
    # bts ............................................................
    bt = hio.read_var_from_hdf(fname, 'bt_108') / 100.
    time_list = hio.read_var_from_hdf(fname, 'time') 

    ntime = bt.shape[0]


    # and geo ref ....................................................
    geo =  hio.read_dict_from_hdf('%s/georef-de.h5' % ad)
    lon, lat = geo['lon'], geo['lat']
    Nmax, L = get_domain_properties(lon, lat)
    # ================================================================


    
    # loop over cluster analysis -------------------------------------
    dmins = [10., 20.]
    threshs = [210., 240.]
    out = {}
    for n in range(ntime):

        tstamp = time_list[n]
        print n, tstamp


        # cloud cloud coverage .......................................
        cov = coverage(bt[n], 200, 240, 5, inverse = True)

        cluster = {}
        for it, thresh in enumerate( threshs ):
        
            tname = 'thresh_%d' % int(thresh)

            # do segmentation ............................................
            c = seg.connectivity_clustering(-bt[n], -thresh)


            # cluster analysis and scai calculation ......................
            ca = cluster_analysis(lon, lat, c, dmin = dmins[it])
            scai = tobin_scai(ca, Nmax = Nmax, L = L)

            
            cluster[tname] = ca
            cluster[tname]['scai'] = scai
        
        out[tstamp] = {}
        out[tstamp]['cluster'] = cluster
        out[tstamp]['coverage'] = cov
    # ================================================================


    # save output ----------------------------------------------------
    outname = '%s/hdcp2/aggregation/msevi-rss-de_scai_tobin_%s.h5' % (
        local_data_path, date_stamp)
    print outname
    hio.save_dict2hdf(outname, out)
    # ================================================================

