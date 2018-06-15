#!/usr/bin/env python

import numpy as np

######################################################################
######################################################################

def get_variable4cellids( cset, varname, index_time, idvec ):

    '''
    Returns variable array based on cell IDs.


    Parameters
    ----------
    cset : dict
        set of cell properties

    varname : str
        variable name for selection

    index_time : int
        time index for temporal subselection

    idvec : list
        cell IDs to be selected

    Returns
    -------
    v : numpy array
        variable array
    '''

    # (i) set time mask
    m = cset['index_time'] == index_time


    # (ii) select variable
    var = cset[varname][m]
    cid = cset['cell_id'][m]

    # (iii) get index of target ids
    ind = np.searchsorted(cid, idvec)
    
    return var[ind]
