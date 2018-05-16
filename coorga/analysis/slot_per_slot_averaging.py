#!/usr/bin/env python

import numpy as np

######################################################################
######################################################################

##LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
## automatically generated by 21-Plot_size-dependent_PCF_from_Bootstrap.ipynb
##TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


def slot_per_slot_averaging(f, d, mask = True, identifier_name = 'time_id'):

    '''
    Uses an identifier ID for categorization and averages a variable f for same IDs.

    Parameters
    ----------
    f : numpy array
       variable that is average for same IDs;
       f can be of rank n <=3, but 1st dimension must be cell dimension
    d : dict
       dictionary (typically cell propertes dict) that contains IDs entries
    mask : boolean numpy array, optional, default = True
       additional mask array to apply conditioned averaging
    identifier_name : str, default = 'time_id'
       name of the identifier variable in d
       

    Returns
    -------
    fave : numpy array
       averaged field f conditioned on different ID values
    
    '''

    gts = []
    time_ids = sorted( set(d['time_id']) )

    
    if np.ndim(f) == 2:
        ncell, ncol = f.shape
        nan = np.nan * np.zeros(ncol)

    elif np.ndim(f) == 3:
        ncell, nrow, ncol = f.shape
        nan = np.nan * np.zeros( (nrow,ncol) )
    
    elif np.ndim(f) == 1:
        nan = np.nan
        
    for tid in time_ids:

        # time masking 
        tmask = (d['time_id'] == tid) 
        
        # additional masking
        m = tmask & mask
        
        gm = f[m].mean(axis = 0)
        
        if m.sum() == 0:
            gts.append( nan )
        else:
            gts.append( gm )

    fave = np.ma.masked_invalid( np.array(gts) )
    
    return fave
