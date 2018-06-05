#!/usr/bin/env python

import os, sys, glob, copy

import numpy as np
import scipy.ndimage

import tropy.analysis_tools.segmentation as seg
import tropy.analysis_tools.optical_flow as oflow
import tropy.analysis_tools.grid_and_interpolation as gi
import mahotas

import tropy.io_tools.netcdf as ncio




######################################################################
######################################################################

def shift_objects(c, ind_shift, 
                  index = None, 
                  as_uv = True):



    '''
    Shifts objects in a segmented cluster field by a determined 
    index shift.


    Parameters
    ----------
    c : numpy array, 2dim
        segmented field, objects cell are indexed
    
    ind_shift : numpy array, 2dim
        index shift (forward shift)   = (row shift, column shift)

    index : list or numpy array, optional, default = None
        index of cell to be shifted

    as_uv : bool, optional, default = True
        if True shift is composed of (col_shift, row_shift)


    Returns
    --------
     cs : numpy array, 2dim
        shift cluster field
    '''



    # get field dimensions
    nrow, ncol  =  c.shape


    # get slices .....................................................
    slices = scipy.ndimage.measurements.find_objects(c)


    # set index
    if index == None:
        index = range(1, len(slices) + 1)


    ishift = ind_shift.astype(np.int)

    # loop over slices ...............................................
    cs = np.zeros_like(c)

    for i in index:
        
        if as_uv:
            dicol, dirow = ishift[i]
        else:
            dirow, dicol = ishift[i]

        try:
            # get object slice
            s = slices[i - 1]


            # shift row
            ir1, ir2 = s[0].start,  s[0].stop

            ir1 += dirow
            ir2 += dirow

            # shift col
            ic1, ic2 = s[1].start,  s[1].stop

            ic1 += dicol
            ic2 += dicol

            s_shift = ( slice(ir1, ir2), slice(ic1, ic2) )

            # check bounds
            if s_shift[0].start < 0 or s_shift[0].stop >= nrow:
                raise Exception('Object #%d crossed row bounds' % i)

            if s_shift[1].start < 0 or s_shift[1].stop >= ncol:
                raise Exception('Object #%d crossed column bounds' % i)

            cc = c[s]

            

            # mask
            m = (cc == i)  # assume that slices are sorted 

            cs[s_shift] = np.where(m, cc, cs[s_shift])
        except:
            pass


    return cs
  

######################################################################
######################################################################

def object_nowcast(f1, f2,
                   output_symmetric_flow = False,
                   output_average_flow = False,
                   output_uv = False,
                   s1 = None, 
                   cluster_method = 'connect',
                   thresh = 0.,
                   vmin = None,
                   vmax = None):




    '''
    Objects identified in field f1 (segmentation s1) are shifted
    with average flow to match objects identified in field f2.

    
    Parameters
    ----------
    f1 : numpy array, 2dim
       field at present time (to be shifted)

    f2 : numpy array, 2dim
       field at future time (taken to calculated opt. flow)

    output_symmetric_flow : bool, optional, default = False
       switch if local flow field is returned

    output_average_flow : bool, optional, default = False
       switch if local flow field is returned

    output_uv : bool, optional, default = False,
       switch if flow vector (shift per object) is returned

    s1 : numpy array, optional, default = None
       segmentation of field s1

    cluster_method : str, optional, default = 'connect'
       which method used for clustering 

    thresh : float, optional, default = 0.
       selected threshold for clustering

    vmin : float, optional, default = None
       lower minimum value for field clipping 
       if None, minimum from input fields is taken

    vmax : float, optional, default = None
       upper maximum value for field clipping 
       if None, maximum from input fields is taken


    Returns
    -------
    t1 : numpy array, 2dim
       transformed object field at present time with object shifted 

    flow : numpy array, 2dim, optional if output_symmetric_flow = True
       local flow field

    mflow : numpy array, 2dim, optional if output_average_flow = True
       cell_average flow field

    uv : numpy array, 2dim, optional if output_uv = True
       shift vector per cell
    '''

    
    
    # clip the range of the fields -----------------------------------
    if vmin == None:
        vmin = np.min([f1.min(), f2.min()])
    if vmax == None:
        vmax = np.max([f1.max(), f2.max()])


    f1c = np.clip(f1, vmin, vmax)
    f2c = np.clip(f2, vmin, vmax)
    # ================================================================



    # get symmetric flow field ---------------------------------------
    flow12 = oflow.displacement_from_opt_flow(f1c, f2c, 
                                              vmin = vmin, 
                                              vmax = vmax)

    flow21 = oflow.displacement_from_opt_flow(f2c, f1c, 
                                              vmin = vmin, 
                                              vmax = vmax)

    aflow = np.maximum( np.abs(flow21), np.abs(flow12) )
    flow = np.sign(flow12) * aflow
    # ================================================================



    # do segmentation if needed --------------------------------------
    if s1 == None:
        s1 = seg.clustering(f1c, thresh, cluster_method = cluster_method)
    # ================================================================


    # get cell-average flow ------------------------------------------
    mflow = np.zeros_like( flow )
    uv = []
    

    for i, f in enumerate([flow[:,:,0], flow[:,:,1]]):
        uvec = scipy.ndimage.measurements.mean(f, 
                                               labels = s1, 
                                               index = range(0, s1.max() + 1))

        uvec[0] = 0.

        uv.append( uvec )

        mflow[:,:,i] = uvec[s1]

    uv = np.column_stack(uv)
    # ================================================================



    # shift cluster field at previous time (like nowcasting step) ----
    t1 = shift_objects(s1, uv, as_uv = True)
    # ================================================================
 

    if output_average_flow:
        return t1, mflow
    elif output_symmetric_flow:
        return t1, flow
    elif output_uv:
        return t1, uv
    else:
        return t1

######################################################################
######################################################################

def get_area_rate(lon, lat, f1, f2, thresh,
                  cluster_method = 'connect',   # segmentation method 
                  nsub = 4,                     # subsampling factor
                  vmin = None,
                  vmax = None,
                  nedge = 20,
                  output_percentile_change = True,
                  percentiles = [50, 75, 90, 95],
                  return_vector = False, 
                  dt = 3600., **kwargs):

    '''
    Calculates area rate between two sets of time-connected objects. 

    Method:
  
    (i) The two fields f1 and f2 are segmented using threshold thresh and condition
    that f1, f2 > thresh (foreground). 

    (ii) Optimal flow transformation is calculated between the field and an average 
    shift is applied to field f1. 
    
    (iii) The shifted f1 and f2 are stack and segmented again to get time connection.


    
    Parameters
    ----------
    lon : numpy array, 2dim
       longitude field

    lat : numpy array, 2dim
       latitude field

    f1 : numpy array, 2dim
       field at present time (to be shifted)

    f2 : numpy array, 2dim
       field at future time (taken to calculated opt. flow)

    thresh : float
       selected threshold for clustering

    cluster_method : str, optional, default = 'connect'
       which method used for clustering 

    nsub : int, optional, default = 4
       subsampling factor

    vmin : float, optional, default = None
       lower minimum value for field clipping 
       if None, minimum from input fields is taken

    vmax : float, optional, default = None
       upper maximum value for field clipping 
       if None, maximum from input fields is taken

    output_percentile_change :  bool, optional, default = True
       switch if change in percentile values is also returned

    percentiles :  list,  optional, default = [50, 75, 90, 95]
       percentiles for which change is monitored

    output_vector : bool, optional, default = False
       switch if change is returned as vector (and not 2d field)

    dt : float, optional, default = 3600.
       time interval

    **kwargs: dict
       keyword argument used in clustering routine


    Returns
    -------
    da : numpy array, 2dim, shape subsampled with nsub
       area rate field ( units km * m/s )

    davec : numpy array 1dim, optional if output_vector = True
       area rate vector, sorted per cell

    dp : numpy array, 2dim, optional if output_percentile_change = True
       percentile change (not divided by dt)

    dpvec : numpy array 1dim, optional if output_vector = True AND output_percentile_change = True
       percentile change vector, sorted per cell
    '''


    # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    # Section 1: prepare input fields
    # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


    # subsampling ----------------------------------------------------
    f1c = f1[::nsub, ::nsub]
    f2c = f2[::nsub, ::nsub]
    
    clon = lon[::nsub, ::nsub]
    clat = lat[::nsub, ::nsub]
    # ================================================================
    


    # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    # Section 2: segmenation of individual fields
    # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


    # perform single field segmentation ------------------------------
    s1 = seg.clustering(f1c, thresh, cluster_method = cluster_method, **kwargs)
    s2 = seg.clustering(f2c, thresh, cluster_method = cluster_method, **kwargs)

    # s1 = mahotas.labeled.remove_bordering(s1, (nedge / nsub, nedge / nsub))
    # s2 = mahotas.labeled.remove_bordering(s2, (nedge / nsub, nedge / nsub))
    # ================================================================
    


    # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    # Section 3: optical flow 
    # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    t1, flow = object_nowcast(f1c, f2c, s1 = s1, output_symmetric_flow = True)


    # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    # Section 4: 3d segmentation (time-connected clusters)
    # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    # stack shifted, previous and actual field
    s3d = np.array([t1, s2])

    s = seg.clustering(s3d, 0,  cluster_method = cluster_method, **kwargs)
    s = mahotas.labeled.remove_bordering(s, (0, nedge / nsub, nedge / nsub))



    # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    # Section 5: calculate output fields
    # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    # area rate
    a = gi.simple_pixel_area(clon, clat)
    

    a1 = scipy.ndimage.measurements.sum(a, 
                                        labels = s[0], 
                                        index = range(0, s.max() + 1))

    a2 = scipy.ndimage.measurements.sum(a,
                                        labels = s[1], 
                                        index = range(0, s.max() + 1))

    
    davec = a2 - a1
    davec[0] = 0.
    
    da = davec[ s[1] ] * 1000. / dt   # units convection km / h to m / s



    # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    # Section 6: optionally calculate percentiles of field
    # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    if output_percentile_change:
        
        
        # PROBLEM: f2c is not transform -> inconsistent with s[0] !!!! SOLVED !!!!!
        
        ft1 = oflow.morph_trans_opt_flow(f1c, flow)
        p1 = percentiles_from_cluster(ft1, s[0], p = percentiles, 
                                      index = range(0, s[1].max() + 1))


        p2 = percentiles_from_cluster(f2c, s[1], p = percentiles,
                                      index = range(0, s[1].max() + 1))

        dpvec = np.row_stack(p2) - np.row_stack(p1)
    
        dp = dpvec[ s[1] ].transpose(2, 0, 1)

        
        if return_vector:
            return davec, dpvec
        else:
            return da, dp
    else:
        if return_vector:
            return davec
        else:
            return da




######################################################################
######################################################################

def percentiles_from_cluster(f, c, p = [25, 50, 75], index = None):

    '''
    Calculates percentiles of cells in an segmented (labeled) field.

    Functionality missing in scipy.ndimage.measurements.


    Parameters
    ----------
    f : numpy array, 2dim
        the field as basis for percentile calculation

    c : numpy array, 2dim, int
        labeled field

    p : list or numpy array, optional, default =  [25, 50, 75]
        percentiles array

    index : list or numpy array, optional, default = None
        list of cell indices to be analyzed
        if None: all cells will be analyszed


    Returns
    --------
    pc : numpy array, 1dim
        array of percentiles per cell (including background (set to zero))
    '''


    # get slices .....................................................
    slices = scipy.ndimage.measurements.find_objects(c)


    # set index
    if index == None:
        index = range(1, len(slices) + 1)

    # loop over slices ...............................................

    nperc = len(p)
    pc = [np.zeros(nperc),]

    for i in index:

        try:
            s = slices[i- 1]
            cc = c[s]
            fc = f[s]
        
            m = (cc == i)  # assume that slices are sorted 
        
            pvec =  np.percentile( fc[m], p )
        except:
            pvec = np.nan * np.ones(nperc)

        pc.append( pvec )

    return pc
            

######################################################################
######################################################################

if __name__ == '__main__':
    
    
    # input test field
    fname = '/vols/talos/home/fabian/data/icon/narval/msevi_narval_DOM01_20160801.nc'
    b = ncio.read_icon_4d_data(fname, ['bt108'], itime = None)['bt108']
    geo = ncio.read_icon_georef(fname)

    itime = 15
    ns = 2
    da, bperc = get_area_rate(geo['clon'], geo['clat'], 
                              -b[itime - 1], -b[itime], -230,
                               vmin = -260,
                              vmax = -180,
                              cluster_method = 'watershed_merge',
                              ctype = 4,
                              filter_method = 'curve',
                              marker_field = 'dist',
                              marker_method = 'iterative_shrinking',
                              numberOfIterations = 5,
                              cluster_masking = True,
                              exclude_border = True,
                              min_size = 3,
                              nsub = ns)
