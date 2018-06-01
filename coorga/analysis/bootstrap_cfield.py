#!/usr/bin/env python

# load libraries -----------------------------------------------------
import sys, os, glob

import numpy as np
import scipy.ndimage

import analysis_tools.grid_and_interpolation as gi
import analysis_tools.segmentation as seg
import io_tools.hdf as hio

from standard_config import *

sys.path.append('%s/.python' % proj_path)
import convective_orga.inout.data_reader as data_reader



######################################################################
######################################################################

def random_field_generator(c, use_poisson = False, nedge = 10, dmin = 1.5):

    '''
    Takes the cells from the cluster field c and rearranges the cell randomly. If option
    is set, the number of cells is taken from a Poisson distribution.


    
    Parameters
    ----------
    c : numpy array, dtype int 
      cluster field (integer index of cell)

    use_poisson : boolean, optional, default: False
       if True, the number of cells is drawn from Poisson distribution, 
       with repeated use of same cells possible


    Returns
    -------
    cran : numpy array
       randomly rearanged cluster field
    '''
    


    # first sort cluster and remove possibly missing ones ------------
    cs = seg.sort_clusters(c)

    cmax = cs.max()
    nrows, ncols = cs.shape
    # ================================================================


    # initialize new random cluster field ----------------------------
    cran = np.zeros_like(cs)
    # ================================================================

    
    # set up the new cluster index list ------------------------------
    if use_poisson:
        ncell = np.random.poisson(cmax)
        cindex_list = np.random.randint(1, cmax + 1, size = ncell)
    else:
        cindex_list = range(1, cmax + 1)
    # ================================================================

    
    # place cells into the new field ---------------------------------
    for cind in cindex_list:


        # (i) cutout cell ............................................
        ndist = dmin + 2
        ccut = gi.cutout_cluster(cs, cind, nedge = ndist)
        mcut = (ccut == cind)

        # outer distance
        dout =  scipy.ndimage.distance_transform_edt(mcut)
        mout = dout < dmin

        # cell dimensions
        nc_row, nc_col = ccut.shape


        # (ii) randomly draw center position
        row_edge = nc_row / 2 + 1
        col_edge = nc_col / 2 + 1

        
        # START ITERATION
        Niter = 0
        Niter_max = 10
        CELL_HAS_OVERLAP = True

        while CELL_HAS_OVERLAP:

            irow = np.random.randint(row_edge, nrows - row_edge)
            icol = np.random.randint(col_edge, ncols - col_edge)

            ir1 = irow - nc_row / 2
            ir2 = irow + nc_row / 2  +  np.mod(nc_row, 2)   # odd vs. even thing
            ic1 = icol - nc_col / 2
            ic2 = icol + nc_col / 2  +  np.mod(nc_col, 2)
        
            # (iii) check if cell is already there
            cran_cut = cran[ir1:ir2, ic1:ic2]

            if cran_cut[mout].sum() == 0:
                CELL_HAS_OVERLAP = False
            
            Niter += 1
            
            if Niter > Niter_max:
                print 'Error: Iteration did not lead to isolated cell position!'
                break 


        # (iv) place cell into the new field
        cran[ir1:ir2, ic1:ic2] = np.where( mcut, ccut, cran_cut)
    # ================================================================

    cran = seg.sort_clusters(cran)

    return cran

######################################################################
######################################################################

    

if __name__ == "__main__":


    fname = sys.argv[1]
    
    c3d = data_reader.read_cluster(fname)

    c = c3d[24]

    cran = random_field_generator(c)
