#!/usr/bin/env python

import os, sys
import numpy as np
import scipy.stats
import analysis_tools.grid_and_interpolation as gi

######################################################################
######################################################################

def histogram2d_kde( vlist, bins,
                     edge_based_grid = False,
                     output_standard = False):

    '''
    Make 2d histograms based on KDE (kernel density estimation) method.

    
    Parameters
    ----------
    vlist : list of numpy arrays
       list of two variables
    
    bins : list or tuple of two numpy arrays
        bin list that is used for binning the two variables

    edge_based_grid : bool, optional, default = False
        controls if edge- or midpoint-based grid is returned
    
    output_standard : bool, optional, default = False
        if True, also the standard numpy histogram results are returned


    Returns
    --------
    x1 : numpy array, 2dim
        grid of the first variable (edge- or midpoint-based)

    y1 : numpy array, 2dim
        grid of the second variable (edge- or midpoint-based)
    
    h : numpy array, 2dim
       kernel density estimates
    
    h_standard : numpy array, 2dim, optional, if output_standard = True
       numpy histogram density estimates (using normed = True)

    
    
    '''

    # prepare values .................................................
    values = np.vstack( vlist )

    # kde calculations ...............................................
    kernel = scipy.stats.gaussian_kde(values)

    # standard histogram .............................................
    h_standard, x1e, x2e = np.histogram2d(vlist[0], vlist[1], bins,
                                          normed = True)

    # the grid .......................................................
    x1g, x2g = np.meshgrid(gi.lmean(x1e), gi.lmean(x2e), indexing='ij')
    x1e, x2e = np.meshgrid(x1e, x2e, indexing='ij')

    # calculate kde based densities
    positions = np.vstack([x1g.ravel(), x2g.ravel()])        
    h = np.reshape(kernel.evaluate(positions).T, x1g.shape)

    # prepare output .................................................
    if edge_based_grid:
        x1, x2 = x1e, x2e
    else:
        x1, x2 = x1g, x2g

    if output_standard:
        return x1, x2, h, h_standard
    else:
        return x1, x2, h


######################################################################
######################################################################

if __name__ == '__main__':

            
    r = np.random.randn(100)
    r2 = r + 0.5 * np.random.randn(100)

    x, y, h = histogram2d_kde([r, r2], (30,20))
    x, y, h, hs = histogram2d_kde([r, r2], (30,20), output_standard = True)
