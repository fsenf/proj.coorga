#!/usr/bin/env python

import numpy as np

######################################################################
######################################################################


def get_cluster_id(v, variable_set):

    '''
    Make an ID array based on a discretization of a variable array.


    Parameters
    ----------
    v : numpy array
       variable array
    
    variable_set : list
       list of range bin edges
    
    
    Returns
    --------
    ids : numpy array, dtype = int
       variable ids
    '''


    # how many variable bins
    nbins = len( variable_set )

          
    # make ids
    ids = np.zeros_like( v ).astype( np.int )
    
    for n in range(nbins - 1):

        v1 = variable_set[n]
        v2 = variable_set[n + 1]
        
        m_id = ( v >= v1 ) & ( v < v2 )

        ids[m_id] = n

    return ids

######################################################################
######################################################################
