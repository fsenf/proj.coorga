#!/usr/bin/env python

import os, sys, glob
import numpy  as np
import datetime
import scipy.ndimage, scipy.spatial

import io_tools.hdf as hio
import analysis_tools.grid_and_interpolation as gi
from analysis_tools.grid_and_interpolation import make_index_set, cutout_fields
import analysis_tools.segmentation as seg


from standard_config import *



######################################################################
######################################################################

def cutout_cluster(c, nc, 
                   nedge = 30):
    
    '''
    Makes a cutout of a categorial field c for an object of class / number nc.


    INPUT
    =====
    c: 2d categorial field
    nc: category number / class which is chosen
    nedge: optional, number of edge pixels included


    OUTPUT
    =====
    ccut: cutout of categorial field c
    '''


    # get index field ------------------------------------------------
    nrow, ncol = c.shape
    irow, icol = make_index_set(nrow, ncol)
    # ================================================================


    # make mask and mask the indices ---------------------------------
    mask = (c== nc)

    ir = irow[mask]
    ic = icol[mask]
    # ================================================================

    
    # determine cutout region ----------------------------------------
    ir1 = ir.min() - nedge
    ir2 = ir.max() + nedge
    
    ic1 = ic.min() - nedge
    ic2 = ic.max() + nedge

    # check bounds
    if ir1 < 0:
        ir1 = 0

    if ir2 >= nrow:
        ir2 = nrow 

    # check bounds
    if ic1 < 0:
        ic1 = 0

    if ic2 >= ncol:
        ic2 = ncol 

    region = ((ir1, ir2), (ic1, ic2))
    # ================================================================


    # do the cutout --------------------------------------------------
    ccut = cutout_fields(c, region)
    # ================================================================

    return ccut

######################################################################
######################################################################

if __name__ == '__main__':


    # get input date -------------------------------------------------
    fname = sys.argv[1]
    basename = os.path.splitext(os.path.basename(fname))[0]

    date = basename.split('_')[-1]

    # check if its is obs or sim?
    ftype = basename.split('_')[0]
    if ftype == 'msevi':
        subpath = None
    elif ftype == 'synsat':
        subpath = 'synsat_oper'

    thresh = 230.
    min_size = 40.

    t0 = datetime.datetime.strptime(date, '%Y%m%d')
    # ================================================================


    # read data ------------------------------------------------------

    # read bt108 .....................................................
    print '... read BT10.8 from %s' % fname
    b3d = hio.read_var_from_hdf(fname, 'IR_108', subpath = subpath) / 100.
    b3d = np.ma.masked_less(b3d, 100.)
     

    # get field dimensions
    ntime, nrow, ncol = b3d.shape
    
    b = b3d[10]
    # ================================================================
    

    # test clustering ------------------------------------------------
    cw = seg.connectivity_clustering(-b, -thresh)

    cw = seg.remove_small_clusters(cw, min_size)

    ccut = cutout_cluster(cw, 1)
    # ================================================================
