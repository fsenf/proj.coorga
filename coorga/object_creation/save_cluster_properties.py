#!/usr/bin/env python

import os, sys, glob, copy
import numpy  as np
import datetime

import tropy.io_tools.hdf as hio
import tropy.analysis_tools.grid_and_interpolation as gi
from tropy.standard_config import *

import coorga.inout.data_reader as data_reader
import cluster_analysis 


######################################################################
######################################################################


def main(fname, 
         expname = 'basic', 
         varname = 'bt108',
         do_output = True):

    '''
    Performed all the cluster analysis steps
    (1) input of data
    (2) analysis
    (3) output into file (optiional)


    INPUT
    =====
    fname: 
    varname: variable name of field that is segmented

    expname: optional, name of parameter set that is used for segmentation.
    
    
    OUTPUT
    ======
    segmented_field: stacked of segmented data
    cset: set of cell properties

    '''


    # data input -----------------------------------------------------
    din = data_reader.input(fname, varname)

    USE_EXISTING_CLUSTER_DATA = 'clustname' in din.keys()
    # ================================================================



    # cluster analysis part ------------------------------------------
    setup_for_later_output, b_segmented, cset = cluster_analysis.cluster_analysis(
        din, varname, expname = expname)
    # ================================================================



    # # output ---------------------------------------------------------
    if do_output:
        out = cset
        out['_settings'] = setup_for_later_output
     
     
        if USE_EXISTING_CLUSTER_DATA:
            fname = din.pop('fname', fname)
        basename = os.path.splitext(os.path.basename(fname))[0]    
        
        output_dir = din['input_dir']
        odir = '%s/cluster_properties' % output_dir
        oname = '%s/clust_prop_%s_%s.h5' % (odir, basename, expname)
        print '... save output to %s' % oname
        hio.save_dict2hdf(oname, out)
     
     
        # and segmented field ..............................................
        if  not USE_EXISTING_CLUSTER_DATA:
            out = {}
            out['b_segmented'] = b_segmented.astype(np.uint16)
            out['_settings'] = setup_for_later_output
            
            
            odir = '%s/cluster_properties' % output_dir
            oname = '%s/segmented_%s_%s.h5' % (odir, basename, expname)
            print '... save output to %s' % oname
            hio.save_dict2hdf(oname, out)
        # # ================================================================

    return b_segmented, cset
    

######################################################################
######################################################################





if __name__ == '__main__':


    # get input ------------------------------------------------------
    fname = sys.argv[1]
        

    # set the segmentation parameters by experient name
    try:
        expname = sys.argv[2]
    except:
        expname = 'basic'

    # set the variable name
    try:
        varname = sys.argv[3]
    except:
        varname = 'bt108'
        
    # ================================================================


    main( fname, expname, varname )
