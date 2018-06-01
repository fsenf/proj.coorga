#!/usr/bin/env python

import numpy as np
import scipy.spatial



######################################################################
######################################################################


def nearest_neighbor_distance(x, y, t, k = 2, props = {}):
    
    '''
    The subroutine uses KDTree Algorithm to calculate 
    nearest-neighor distance for k neighbors.


    Parameters
    ----------
    x : numpy array, 1dim with shape = (ncells)
        x-coordinate

    y : numpy array, 1dim with shape = (ncells)
        y-coordinate
    
    t : numpy array, 1dim with shape = (ncells)
        t-coordinate, which is used for sequential masking

    props : dict, optional, default = {}
        set of properties used of nearest neighbors calculations, 
        e.g. what is the size of the nearest neighbor?
    

    Returns
    -------
    dist : numpy array, 2dim with shape = (ncells, number_of_neighbors)
       nearest neighbor distances field

    nn_props : dict containing numpy arrays, optional if props is not empty
       collection pof nearest neighbor properties
    '''

    tvec = set(t)
    
    dist = []

    nn_props = {}
    for pname in props.keys():
        nn_props[pname] = []

    # loop over different time instances
    for ti in sorted(tvec):
        
        mask = (t == ti)
        
        xy = np.column_stack([x[mask], y[mask]])
        kdtree =  scipy.spatial.KDTree(xy)
        
        d, index = kdtree.query(kdtree.data, k = k + 1)
        
        dist.append( d )

        for pname in props.keys():
            nn_props[pname].append( props[pname][mask][index[:,1]] )

    dist = np.row_stack(dist )

    for pname in props.keys():
        nn_props[pname] = np.hstack( nn_props[pname] )

    if len(nn_props.keys()) == 0:
        return dist[:,1:]
    else:
        return dist[:,1:], nn_props
        

######################################################################
######################################################################
