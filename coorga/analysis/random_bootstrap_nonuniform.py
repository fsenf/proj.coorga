
##LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
## automatically generated by 04-Bootstrap_non-uniform_NumberDensity.ipynb
##TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


# load libraries -----------------------------------------------------
import numpy as np
import scipy.ndimage

import analysis_tools.grid_and_interpolation as gi
import analysis_tools.segmentation as seg
import analysis_tools.statistics as stats


##############################################################################
##############################################################################


def random_field_generator_nonuniform_dist(c,  nd0, 
                                               nfixed = None,
                                               use_poisson = False, 
                                               nedge = 10, 
                                               dmin = 1.5):

    '''
    Takes the cells from the cluster field c and rearranges the cell randomly. If option
    is set, the number of cells is taken from a Poisson distribution.



    Parameters
    ------------
    c : numpy array, int
       cluster field (integer index of cell)

    nd0 : numpy array
        average number distribution (same grid as cluster field c !!!)
    
    nfixed : int, optional, default = None
        if set, nfixed specifies a constant number of cells used in bootstrapping

    use_poisson : bool, optional, default = False
       if True, the number of cells is drawn from Poisson distribution, 
       with repeated use of same cells possible

    nedge : int, optional, default = 10
       size of edge (in px) where no cells are allowed
    
    dmin : float, optional, default = 1.5
       minimum distance between two cells randomly place in the domain

    
    Returns
    --------
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
    if not nfixed == None:
        ncell = nfixed
        cindex_list = np.random.randint(1, cmax + 1, size = ncell)
    elif use_poisson:
        ncell = np.random.poisson(cmax)
        cindex_list = np.random.randint(1, cmax + 1, size = ncell)
    else:
        ncell = cmax
        cindex_list = range(1, cmax + 1)
    # ================================================================


    # draw random numbers in advance ---------------------------------
    irows = np.arange(0, nrows + 1) - 0.5
    icols = np.arange(0, ncols + 1) - 0.5
    bins = (irows, icols)
    
    fac = 100  # hope that 20 times the cell number is okay ...
    random_position = stats.draw_from_empirical_dist(bins, nd0, 
                                                     Nsamp = fac * ncell,
                                                     discrete_version = True)
    
    # get the counter axis in-front
    random_position = np.array( random_position ).T
    # ================================================================

    
    
    # place cells into the new field ---------------------------------
    ncounter = 0
    for icell, cind in enumerate(cindex_list):


        # (i) cutout cell ............................................
        ndist = dmin + 2
        ccut = gi.cutout_cluster(cs, cind, nedge = ndist)
        mcut = (ccut == cind)

        # outer distance
        #        dout =  scipy.ndimage.distance_transform_edt(mcut)
        # BUGFIX: consider outer edge --> dilatation NOT shrinking !!!
        dout =  scipy.ndimage.distance_transform_edt( mcut == False )

        mout = dout < dmin

        # cell dimensions
        nc_row, nc_col = ccut.shape


        # (ii) randomly draw center position
        row_edge = nc_row / 2 + 1
        col_edge = nc_col / 2 + 1

        
        # START ITERATION
        Niter = 0
        Niter_max = fac
        CELL_HAS_OVERLAP = True

        while CELL_HAS_OVERLAP:


            point_inside = False
            while not point_inside:
                
                # get random position
                irow, icol = random_position[ncounter]
                ncounter += 1
                
                # check if random point is inside the core domain
                irow_inside = (irow > row_edge) & (irow < nrows - row_edge)
                icol_inside = (icol > col_edge) & (icol < ncols - col_edge)
                point_inside = irow_inside & icol_inside
                
            
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
        # cran[ir1:ir2, ic1:ic2] = np.where( mcut, ccut, cran_cut)
        cran[ir1:ir2, ic1:ic2] = np.where( mcut, icell + 1, cran_cut)  #set increasing index of random cell placement
    # ================================================================

    cran = seg.sort_clusters(cran)

    return cran

######################################################################
######################################################################
