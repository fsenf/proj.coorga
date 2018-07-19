#!/usr/bin/env python

import os, sys, glob

import numpy as np
import scipy.ndimage

import tropy.analysis_tools.grid_and_interpolation as gi
from coorga.metrics.pair_correlation_function import radius_ring_edges, pairCorrelationFunction_2D


######################################################################
######################################################################

def calculate_pcf4clust( c, egrid, nd_ref, 
                         weight = 1., 
                         rMax = 500., 
                         pcf_sensitivity = True,
                         dr = 20.,
                         verbose = False):


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
    
    verbose : bool, optional, default = False
        switch if more standard output is provided


    Returns
    --------
    c : dict
        set of cluster props with pcf field updated
    '''


    # read input data ------------------------------------------------
    ncells = len( c['abs_time'] )
    # ================================================================

    
    # get norm of cell density function ------------------------------
    dx = gi.lmean(np.diff(egrid[0], axis = 0), axis = 1)
    dy = gi.lmean(np.diff(egrid[1], axis = 1), axis = 0)

    ncells_per_slot = (dx * dy * nd).sum()
    nd_normed = nd_ref / ncells_per_slot
    # ================================================================

    
    # and also define homogeneous cell density function --------------
    nd_const = 1. / (dx * dy * np.size(nd_ref))
    # ================================================================
    
    
    # get cell coordinates and cutout region -------------------------
    xt, yt = c['x_mean'], c['y_mean']
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
        

        pos = np.array([xt[m], yt[m]]).T

        if verbose:
            print tid, m.sum()

        # correct for temporal variation of cell number
        nd_corr = ccount * nd_normed
        
        # pair-correlation analysis with variable BG density
        gfull, g, rbins, indices = pairCorrelationFunction_2D(pos, egrid, nd_corr, rMax, dr, weight = weight)

        gfull_with_edge = np.nan *  np.ma.ones( (ccount, nbins - 1) )
        gfull_with_edge[indices] = gfull
        c['pcf'][m] = gfull_with_edge
        c['pcf'] = np.ma.masked_invalid( c['pcf'] )
        
        if pcf_sensitivity:
            
            nd = {}
            nd['pcf_hom'] = ccount * nd_const
            nd['pcf_hom_tconst'] = ncells_per_slot * nd_const
            nd['pcf_tconst'] = ncells_per_slot * nd_normed
            
            for pcf_name in nd.keys():
                
                # pair-correlation analysis with variable BG density
                gfull, g, rbins, indices = pairCorrelationFunction_2D(pos, 
                                                                      egrid, 
                                                                      nd[pcf_name],
                                                                      rMax, dr, weight = weight)

                gfull_with_edge = np.nan *  np.ma.ones( (ccount, nbins - 1) )
                gfull_with_edge[indices] = gfull
                
                c[pcf_name][m] = gfull_with_edge
                c[pcf_name] = np.ma.masked_invalid( c[pcf_name] )


    c['rbins_pcf'] = rbins
    # ================================================================



    return c
    # ================================================================
