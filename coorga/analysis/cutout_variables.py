#!/usr/bin/env python

import numpy as np

######################################################################
######################################################################

def cutout_variables(d, vname, vrange_value):
    
    '''
    Cuts out cluster data for a certain variable value or range.
    

    Parameters
    ----------
    d : dict
        set of cell properties

    vname : str
        variable name

    vrange_value : EITHER list or tuple of float OR only one value
        variable range vrange = (v1, v2), OR variable value


    Returns
    -------
    dcut : dict
        cutout of of cell set
    '''


    # get length of range
    nvar = len( vrange_value )

    if nvar == 1:
        return cutout_variable_value(d, vname, vrange_value)
    elif nvar == 2:
        return cutout_variable_range(d, vname, vrange_value)
    
######################################################################
######################################################################

def cutout_variable_range(d, vname, vrange):
    
    '''
    Cuts out cluster data for a certain variable range.
    

    Parameters
    ----------
    d : dict
        set of cell properties

    vname : str
        variable name

    vrange : list or tuple of float
        variable range vrange = (v1, v2)


    Returns
    -------
    dcut : dict
        cutout of of cell set
    '''
    
    # get variable field
    v = d[vname]
    
    # create mask
    v1, v2 = vrange
    m = (v >= v1) & (v < v2)

    # and cut the cluster data
    dcut = {}
    for k in d.keys():
        try:
            dcut[k] = d[k][m]
        except:
            print 'No masking applies for %s' % k
    
    return dcut

######################################################################
######################################################################

def cutout_variable_value(d, vname, value):
    
    '''
    Cuts out cluster data for a certain discrete variable value.
    

    Parameters
    ----------
    d : dict
        set of cell properties

    vname : str
        variable name

    value : int, float, or str
        variable value used for masking


    Returns
    -------
    dcut : dict
        cutout of of cell set
    '''
    
    # get variable field
    v = d[vname]
    
    # create mask
    m = (v == value)

    # and cut the cluster data
    dcut = {}
    for k in d.keys():
        try:
            dcut[k] = d[k][m]
        except:
            print 'No masking applies for %s' % k
    
    return dcut
