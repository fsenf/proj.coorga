#!/usr/bin/env python

import sys, os, glob
import datetime
import numpy as np

import tropy.io_tools.hdf as hio
import tropy.analysis_tools.grid_and_interpolation as gi

from coorga.analysis.average_numberdensity import calculate_average_numberdensity

######################################################################
######################################################################

def main(collection_file, 
                do_output = True, 
                output_file = None,
                smooth_sig = 2.,
                filter_in_logspace = True,
                var_names = None,
                bins = None,
                itime_range = None,
                dx = 20. ,
                dy = 20.,):

    
    '''
    Calculates average number density.
    

    Parameters
    ----------
    collection_file : str
        filename of cluster collection

    do_output : bool, optional, default = True
        switch if output is written

    output_file : str, optional, default = None
        name of output file
        if None: name is generated based on collection filename

    smooth_sig : float, optional, default = 2.
       sigma of Gaussian filter for smoothing the results

    filter_in_logspace : bool, optional, default = True,
       switch if filter is applied in log-space

    var_names : list of two strings, optional, default = None
       list of the two variables used for binning
       if None: 'x_mean' and 'y_mean' are used
       
    bins : list or tuple of two numpy arrays, optional, default = None
       sets the bins of histogram analysis
       if None: min and max is estimated from data and (dx, dy) is used
       
    itime_range : list or tuple of two int, optional, default = None
       selects a range of relative time for the analysis
       if None: all times are used
    
    dx : float, optional, default = 20. 
       interval of x-binning if bins == None

    dy : float, optional, default = 20.
       interval of y-binning if bins == None


    Returns
    --------
    egrid : tuple of two numpy arrays
       edge-based output grid

    nd : numpy array, 2dim
       number density field.

   
    '''
    # ================================================================
        

    # input collection data ------------------------------------------
    cset = hio.read_dict_from_hdf( collection_file ) 
    # ================================================================
 

    # do number denisty calculation ----------------------------------
    egrid, nd_ref = calculate_average_numberdensity(cset,
                                    smooth_sig = smooth_sig,
                                    filter_in_logspace = filter_in_logspace,
                                    var_names = var_names,
                                    bins = bins,
                                    itime_range = itime_range,
                                    dx = dx ,
                                    dy = dy )

    
    # ================================================================
    

    # save stuff into hdf --------------------------------------------
    if do_output:
        out = {}
        out['xgrid'] = egrid[0]
        out['ygrid'] = egrid[1]
        out['number_density'] = nd_ref


        if output_file is None:
            oname = collection_file.replace('collected_cluster_props',
                                            'average_number_density')
        else:
            oname = output_file

        print '... save data to %s' % oname
        hio.save_dict2hdf(oname, out)
    # ================================================================

    return egrid, nd_ref



######################################################################
######################################################################
