#!/usr/bin/env python

import sys, os, glob
import datetime
import numpy as np

import tropy.io_tools.hdf as hio
import tropy.analysis_tools.grid_and_interpolation as gi

from coorga.analysis.average_numberdensity import calculate_average_numberdensity
import coorga.inout.data_reader as data_reader

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
                dy = 20.,
                do_interpolation = False,
                variable_filename = None,   
                variable_name = None,
                input_options = {}):

    
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

    do_interpolation : bool, optional, default = True
        if number density field is interpolated onto cluster grid

    variable_filename : str, optional, default = None
        one of the variable files to input cluster grid georef
                
    variable_name : str, optional, default = None
        variable name in variable file

    input_options : dict, optional, default = {}
        all possible input options used in the data_reader


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

    xgrid = egrid[0]
    ygrid = egrid[1]
    # ================================================================
    
    # do interpolation of reference number density on cluster grid ---
    if do_interpolation:
        if variable_filename is not None and variable_name is not None:
            varset = data_reader.input( variable_filename, 
                                        variable_name,
                                        **input_options )

            # centered base grid
            xgridc = gi.lmean( gi.lmean( xgrid, axis = 0), axis = 1)
            ygridc = gi.lmean( gi.lmean( ygrid, axis = 0), axis = 1)


            # get target grid
            cx = varset['x']
            cy = varset['y']

            # interpolation index
            ind = gi.create_interpolation_index( xgridc, ygridc, 
                                                 cx, cy, xy = True)

            # get total number of cells
            ncells = (nd_ref * dx * dy ).sum()


            # normalization of new nd field
            da = gi.simple_pixel_area( cx, cy, xy = True)
            norm = ( nd_ref[ind] * da ).sum()
            nd_int = ncells * nd_ref[ind] / norm
    # ================================================================


    # save stuff into hdf --------------------------------------------
    if do_output:
        out = {}
        out['xgrid'] = xgrid
        out['ygrid'] = ygrid
        out['number_density'] = nd_ref

        if do_interpolation:
            out['x_int'] = cx
            out['y_int'] = cy
            out['nd_int'] = nd_int


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
