
##LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
## automatically generated by 04-Bootstrap_non-uniform_NumberDensity.ipynb
##TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


# load libraries -----------------------------------------------------
import numpy as np
import scipy.ndimage

import tropy.analysis_tools.grid_and_interpolation as gi
import tropy.analysis_tools.segmentation as seg
import tropy.analysis_tools.statistics as stats

import coorga.object_creation.cluster_analysis as cluster_analysis
from coorga.analysis.get_variable4cellids import get_variable4cellids

##############################################################################
##############################################################################


def random_field_generator_nonuniform_dist(c,  nd0, 
                                               nfixed = None,
                                               output_cell_mapping = False,
                                               use_poisson = False, 
                                               Niter_max = 100,
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

    Niter_max : int, optional, default = 100
        number of iterations to try to place a new cell

    dmin : float, optional, default = 1.5
       minimum distance between two cells randomly place in the domain

    output_cell_mapping : bool, optional, default = False
        switch if mapping between new and original cell index is provided 
    

    Returns
    --------
    cell_mapping : dict, optional, if output_cell_mapping = True
        gives the mapping between new and original cell index,
        e.g cell_mapping[ new_index ] = original_index
    
    cran : numpy array
        randomly rearanged cluster field


    '''
    


    # first sort cluster and remove possibly missing ones ------------
    
    # sorting destroys order --> bad
    cs = c # seg.sort_clusters(c)
    existing_cell_indices = np.array( list(set(cs.flatten()) - {0}) )

    nexisting = len( existing_cell_indices )

    nrows, ncols = cs.shape
    # ================================================================


    # initialize new random cluster field ----------------------------
    cran = np.zeros_like(cs)
    # ================================================================

    
    # set up the new cluster index list ------------------------------
    if not nfixed == None:
        ncell = nfixed
        irand = np.random.randint(nexisting, size = ncell)
        cindex_list = existing_cell_indices[irand]
    
    elif use_poisson:
        ncell = np.random.poisson(nexisting)
        irand = np.random.randint(nexisting, size = ncell)
        cindex_list = existing_cell_indices[irand]
 
    else:
        ncell = nexisting
        cindex_list = existing_cell_indices
    # ================================================================


    # draw random numbers in advance ---------------------------------
    irows = np.arange(0, nrows + 1) - 0.5
    icols = np.arange(0, ncols + 1) - 0.5
    bins = (irows, icols)
    
    fac = Niter_max  # hope that 100 times the cell number is okay ...
    random_position = stats.draw_from_empirical_dist(bins, nd0, 
                                                     Nsamp = fac * ncell,
                                                     discrete_version = True)
    
    # get the counter axis in-front
    random_position = np.array( random_position ).T.astype( np.int )
    # ================================================================

    
    
    # place cells into the new field ---------------------------------
    ncounter = 0
    icell = 1
    cell_mapping = {}

    for cind in cindex_list:

        # (i) cutout cell ............................................
        ndist = np.round( dmin + 2 ).astype( np.int )
        try:
            ccut = gi.cutout_cluster(cs, cind, nedge = ndist)
        except:
            # cutout fails if cell number is empty
            print 'Warning: break at cell cutout used!'
            break

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
        cran[ir1:ir2, ic1:ic2] = np.where( mcut, icell, cran_cut)  #set increasing index of random cell placement

        cell_mapping[icell] = cind

        icell += 1
    # ================================================================

    # sorting
    csort = seg.sort_clusters(cran)

    # after sorting also the index map has to be adjusted
    index_sorted = range(1, csort.max() + 1)   # new index

    index_unsorted = scipy.ndimage.measurements.mean(cran, 
                                labels = csort, 
                                index = index_sorted).astype(np.int)

    revised_mapping = []
    for i, i_sort in enumerate( index_sorted ):
        i_unsorted = index_unsorted[i]
        revised_mapping.append( cell_mapping[i_unsorted] )

    if output_cell_mapping:
        return revised_mapping, csort
    else:
        return csort

######################################################################
######################################################################


def calculate_bootstrap4clust(cx, cy, segmented, cset, nd_ref,
                                               nfixed = None,
                                               use_poisson = False,
                                               dmin = 1.5,
                                               addvarnames = []):

    '''
    Calculates a randomization of cluster field c based on a background
    number density field.


    Parameters
    ----------
    cx : numpy array, 2dim
        x-coordinate of cluster field c

    cy : numpy array, 2dim
        y-coordinate of cluster field c

    segmented : numpy array, 3dim, (ntimes, nrows, ncols)
        cluster field

    cset : dict
        set of cell properties

    nd_ref : numpy array, 2dim, shape = (nrows, ncols)
        reference number density same grid as cluster field

    nfixed : int, optional, default = None
        if set, nfixed specifies a constant number of cells used in bootstrapping

    use_poisson : bool, optional, default = False
       if True, the number of cells is drawn from Poisson distribution, 
       with repeated use of same cells possible

    dmin : float, optional, default = 1.5
       minimum distance between two cells randomly place in the domain

    addvarnames : list, optional, default = ['imf_mean']
        list of additional variables add to the bootstrap set


    Returns:
    --------
    cset_ran : dict
        cell data for randomized field

    segmented_ran : numpy, 3dim, (ntimes, nrows, ncols)
        randomized cluster field

    '''
    
    # prepare aux fields
    carea = gi.simple_pixel_area( cx, cy, xy = True )

    
    # initialize random field
    ntimes, nrows, ncols = segmented.shape
    segmented_ran = np.zeros_like( segmented )
    
    cset_ran = {}

    for addvar in addvarnames:
        cset_ran[addvar] = []
    noffset = 0

    # random shifting within time loop
    for n in range( ntimes ):
    
        c = segmented[n]
 
        # random bootstrap 
        if True: #try:
           cell_mapping, cran = random_field_generator_nonuniform_dist(c, 
                                        nd_ref,
                                        use_poisson = use_poisson,
                                        nfixed = nfixed, 
                                        dmin = dmin,
                                        output_cell_mapping = True)

        else: #except:
            cran = np.zeros_like( c )

        segmented_ran[n] = cran[:, :]

        # cluster analysis for bootstrap set ---------------------
        dset = {'clust': cran, 'x': cx, 'y': cy, 'area':carea}

        dset['rel_time'] = get_variable4cellids(cset, 'rel_time', n, [1])[0]
        dset['abs_time'] = get_variable4cellids(cset, 'abs_time', n, [1])[0]
        dset['index_time'] = n

        # get cluster properties
        cluster_analysis.cellset_analysis(dset, cset_ran, 
                                          noffset = noffset, 
                                          var_names = [],
                                          weight_names = [],
                                          do_landsea_fraction = False)

        for addvar in addvarnames:
            cset_ran[addvar] += [ get_variable4cellids(cset, 
                                                addvar, 
                                                n, 
                                                cell_mapping) ]

        noffset += cran.max()
    
    for addvar in addvarnames:
        cset_ran[addvar] = np.hstack( cset_ran[addvar] )

    return cset_ran, segmented_ran

######################################################################
######################################################################

