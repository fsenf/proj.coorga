#!/usr/bin/env python

import os, sys, glob

import numpy as np
import scipy.ndimage

import analysis_tools.grid_and_interpolation as gi
from coorga.metrics.pair_correlation_function import radius_ring_edges, pairCorrelationFunction_2D


######################################################################
######################################################################

def calculate_pcf4clust( c, egrid, nd_ref, weight = 1., rMax = 500., dr = 20.):


    '''
    Calculates Pair-correlation function for a cluster set.


    Parameters
    -----------
    c : dict
        set o cluster properties

    egrid : tuple or list of two numpy arrays
        edge-based grid of reference number density

    nd_ref : numpy array, 2dim
        reference number density

    weight : float, optional, default = 1.
        weight between equi-distant and equal-area range binning
        0 = equal-area
        1 = equi-distant

    rMax : float, optional, default = 500.
        largest range bin

    dr : float, optional, default = 20.
        range ring interval


    Returns
    --------
    c : dict
        set of cluster props with pcf field updated
    '''


    # read input data ------------------------------------------------
    ncells = len( c['abs_time'] )
    # ================================================================



    # get cell coordinates and cutout region -------------------------
    lon, lat = c['lon_mean'], c['lat_mean']
    xt, yt  =  gi.ll2xy(lon, lat, lon0 = 0, lat0 = 0)
    # ================================================================


    # loop over all time ids -----------------------------------------
    if 'time_id' in c.keys():
        tname = 'time_id'
    else:
        tname = 'rel_time'

    time_ids = sorted( set(c[tname]) )
    
    rbins = radius_ring_edges(rMax, dr, weight = weight)
    nbins = len(rbins)

    
    
    c['pcf'] = np.nan * np.ma.ones( (ncells, nbins - 1) )
    

    for tid in time_ids:

        # masking 
        m = (c[tname] == tid) 
        ccount = m.sum()
        # ccount_inside_region = (m & rmask).sum()

        pos = np.array([xt[m], yt[m]]).T
        # print tid, m.sum()


        # pair-correlation analysis with variable BG density
        gfull, g, rbins, indices = pairCorrelationFunction_2D(pos, egrid, nd_ref, rMax, dr, weight = weight)

        gfull_with_edge = np.nan *  np.ma.ones( (ccount, nbins - 1) )
        gfull_with_edge[indices] = gfull
        c['pcf'][m] = gfull_with_edge

    c['rbins_pcf'] = rbins
    # ================================================================

    return c
    # ================================================================
