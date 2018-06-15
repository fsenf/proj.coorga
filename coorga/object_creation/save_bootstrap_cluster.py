#!/usr/bin/env python

import os, sys, glob, copy
import numpy  as np
import datetime

import tropy.io_tools.hdf as hio

from coorga.analysis.random_bootstrap_nonuniform import calculate_bootstrap4clust


######################################################################
######################################################################


def main(clusterfile, 
        number_density_file,
        bootstrap_options = {},
        do_output = True,
        output_segmented = True,
        output_dir = None):

    '''
    Saves Bootstrap Cluster Files.


    Parameters
    ----------
    clusterfile : str
        name of cluster file

    number_density_file : str
        name of number density file

    do_output : bool, optional, default = True
        switch that determines if output is written in predefined file
    
    output_dir : str, optional, default = None
        name of output directory


    Returns
    -------
    None
    '''


    # input clusterdata
    cset = hio.read_dict_from_hdf( clusterfile )

    # read segmented data
    segmentedfile = clusterfile.replace('clust_prop', 'segmented')
    fseg = hio.read_dict_from_hdf(segmentedfile)
    segmented = fseg['b_segmented']

    # get number density
    ndset = hio.read_dict_from_hdf( number_density_file )
    cy = ndset['y_int']
    cx = ndset['x_int']
    nd = ndset['nd_int']

   
    #perform bootstrapping
    cset_ran, segmented_ran = calculate_bootstrap4clust(cx, 
                                                        cy, 
                                                        segmented, 
                                                        cset, 
                                                        nd, 
                                                        **bootstrap_options)


    # # output ---------------------------------------------------------
    if do_output:
        out = cset_ran
     
    
        if output_dir is None:
            odir = '%s/boot_cluster_properties' % din['input_dir']
        else:
            odir = output_dir

        basename = os.path.splitext(os.path.basename( clusterfile ))[0]    
        
        oname = '%s/%s.h5' % (odir, basename)
        print '... save output to %s' % oname
        hio.save_dict2hdf(oname, out)
     
    if do_output and output_segmented:

        # and segmented field ..............................................
        out = {}
        out['b_segmented'] = segmented_ran.astype(np.uint16)
        
        basename = os.path.splitext(os.path.basename( segmentedfile ))[0]
        
        oname = '%s/%s.h5' % (odir, basename)
        print '... save output to %s' % oname
        hio.save_dict2hdf(oname, out)
        # ================================================================

    return cset_ran, segmented_ran
    

######################################################################
######################################################################


