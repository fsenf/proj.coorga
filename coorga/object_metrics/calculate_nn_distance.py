#!/usr/bin/env python

import numpy as np
import scipy.spatial



######################################################################
######################################################################


def calculate_nn_distance(x, y, t, k = 2, props = {}):
    
    '''
    The subroutine uses KDTree Algorithm to calculate 
    nearest-neighor distance for k neighbors.

    USAGE
    =====
    dist = calculate_nn_distance(x, y, t, k = 2)


    INPUT
    =====
    x: x-coordinate
    y: y-coordinate
    t: t-coordinate, which is used for sequential masking

    
    OUTPUT
    ======
    dist: nearest neighbor distances field
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
