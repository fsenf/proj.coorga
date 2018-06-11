#!/usr/bin/env python

import numpy as np
import scipy.spatial


import tropy.io_tools.hdf as hio

from coorga.inout.data_reader import read_narval_addvars

######################################################################
######################################################################


def add_variable2clustdat( cset, addvarset, variable_list, nsub = 10 ):

    '''
    Add a variable values to the cluster set given a gridded variable field.

    
    Parameters
    ----------
    cset : dict
       cluster data set

    addvarset : dict
       set of additional variables


    Returns 
    --------
    cset : dict
       updated cluster data set
    '''
    


    # First get the cell positions
    # ==============================
    cx = cset['x_mean']
    cy = cset['y_mean']
    index_time = cset['index_time']

    cell_pos = np.array( [cx, cy] ).T


    # Initialize Output Fields
    # ==============================
    for vname in variable_list:
        cset[vname] = np.nan * np.ones_like( cx )


    # Then look at the variable grid
    # ==============================
    gx = addvarset['x'][::nsub, ::nsub]
    gy = addvarset['y'][::nsub, ::nsub]
    
    grid_pos = np.array([ gx.flatten(), gy.flatten() ]).T


    # Nearest Neighbor Tree preparation
    # ==============================
    tree = scipy.spatial.KDTree( grid_pos )
    



    # Prepare the loops - time and variable-wise
    # ==============================
    tindex_set = set( index_time )



    for itime in tindex_set:

        # select the right time index
        index_mask = index_time == itime

        # get the nearest grid points
        dist, cindex = tree.query( cell_pos[index_mask, :], 1)


        for vname in variable_list:
            
            # bring variable field into shape
            v = addvarset[vname][itime, ::nsub, ::nsub].flatten()

            # make the nn interpolation
            cset[vname][index_mask] = v[cindex]

    return cset

######################################################################
######################################################################            



if __name__ == '__main__':

    
    fname = '/vols/talos/home/fabian/data/icon/narval/var4workflowtest/cluster_properties/clust_prop_icon-narval_dom01_imf_20160801_exp011.h5'
    cset = hio.read_dict_from_hdf(fname)


    addvar = {}
    vlist = ['mean_tcwv', 'mean_sst', 'mean_cape', 'max_cape']
    for vname in vlist:
        fname = '/vols/talos/home/fabian/data/icon/narval/variables/icon-narval_dom01_%s_20160801.nc' % vname
        addvarset = read_narval_addvars(fname, vname)
        addvar.update( addvarset )

    add_variable2clustdat( cset, addvar, vlist)
