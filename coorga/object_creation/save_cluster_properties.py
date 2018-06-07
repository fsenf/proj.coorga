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
    Performed all the cluster analysis steps.


    Parameters
    ----------
    fname : str
       name of input variable file

    expname : str, optional, default = 'basic'
        name of the experiment (segmentation setup)

    varname: str, optional, default = 'bt108'
       analysis variable name that is read from fname

    do_output: bool, optional, default = True
        switch that determines if output is written in predefined file
    

    Returns
    -------
    b_segmented : numpy array, 3dim
        segmented field 
    
    cset : dict
        dictionary of cell properties


    Notes
    -----
    The following workflow is applied
    (1) input of data
    (2) analysis
    (3) output into file (optiional)

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


