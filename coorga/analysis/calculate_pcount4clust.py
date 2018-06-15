#!/usr/bin/env python

import os, sys, glob

import numpy as np
import scipy.ndimage

from get_cluster_id import get_cluster_id

import tropy.analysis_tools.grid_and_interpolation as gi
from coorga.metrics.pair_correlation_function import radius_ring_edges

from camocapy.mc_tools import calculate_inner_pairnumbers

######################################################################
######################################################################

def calculate_pcount4clust( c, vname, variable_set, 
                            pcount_name = 'pcount', 
                            rMax = 500., 
                            dr = 20.,
                            verbose = False):


    '''
    Calculates Pair-correlation function for a cluster set.


    Parameters
    -----------
    c : dict
        set o cluster properties

    vname : str
        name of variable used for categorization

    variable_set : list
        list of variable bins

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



    # get cell coordinates and cutout region -------------------------
    xt, yt = c['x_mean'], c['y_mean']
    # ================================================================

    
    # get ids --------------------------------------------------------
    var = c[vname]
    full_ids = get_cluster_id(var, variable_set)
    nvar_bins = len( variable_set )

    Ncomponent = nvar_bins - 1
    # ================================================================


    # loop over all time ids -----------------------------------------
    if 'time_id' in c.keys():
        tname = 'time_id'
    else:
        tname = 'rel_time'

    time_ids = sorted( set(c[tname]) )
    
    nbins = int(rMax / dr) + 1 

    c[pcount_name] = np.nan * np.ma.ones((ncells, Ncomponent, nbins))


    for tid in time_ids:

        # masking 
        m = (c[tname] == tid) 
        ccount = m.sum()

        labels = np.arange(ccount)


        # ccount_inside_region = (m & rmask).sum()
        pos = np.array([xt[m], yt[m]]).T
        ids = full_ids[m]



        if verbose:
            print tid, m.sum()

        # get particle count matrix
        rbins, pcount = calculate_inner_pairnumbers(pos, ids, labels, rMax, 
                                                    dr = dr, 
                                                    Ncomponent = Ncomponent)

        c[pcount_name][m] = pcount.transpose(1, 0, 2)
    
    c['%s_rbins' % pcount_name] = rbins
    c['%s_variable_set' % pcount_name] = variable_set
    c['%s_variable_name' % pcount_name] = vname
    c['%s_ids' % pcount_name] = full_ids
    # ================================================================



    return c
    # ================================================================
