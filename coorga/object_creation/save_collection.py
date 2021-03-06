#!/usr/bin/env python

import sys, os, glob
import datetime
import numpy as np
import scipy.spatial

import tropy.io_tools.hdf as hio
import tropy.analysis_tools.grid_and_interpolation as gi
from  tropy.standard_config import *

from  coorga.inout.cluster_prop_reader import read_cluster_props
from coorga.metrics.nearest_neighbor_distance import nearest_neighbor_distance as calculate_nn_distance
from cluster_analysis import create_time_id, remove_too_few_clusters
from segmentation_config import  predefined_collections

######################################################################
######################################################################

def main(flist, 
         do_output = True, 
         output_file = None,
         add_aux_varlist = []):

    '''
    Gathers a collection of cell properties for a selected 
    time period, calculates some new props and finnaly saves results.


    Parameters
    ----------
    flist : list of str
        input filelist

    do_output : bool, optional, default = True
        switch that determines if output is written in predefined file

    output_file : str
        name of output (collection) file
    
    add_aux_varlist : list of strings
        list of additional auxiliary variable names to be included


    Returns
    -------
    out: dict
        dictionary that contains cell properties


    '''
    # ================================================================
        


    # variables list -------------------------------------------------
    vlist =  add_aux_varlist + [
        'area',
        'diameter', 
        'x_mean', 
        'y_mean', 
        'rel_time', 
        'abs_time', 
        'lon_mean',
        'lat_mean',
        'dist_max', 
        'hull_dmax',
        'bt108_min',
        'shape_bt108',
        'shape_smf',
        'area_rate',
        'perc_rate',
        'shape_dist',
        'mass_lift',
        'pcount',
        'pcf',
        'pcf0',
        'pcf0_timevar',
        'rbins',
        'rbins_pcount',
        'rbins_pcf',
        'diameter_set',
        'mean_sst',
        'mean_cape',
        'max_cape',
        'mean_tcwv']
    # ================================================================


    # for rain / water column comparison ------------------------------
    for vname in ['rr', 'tcw', 'imf'] + add_aux_varlist:
        for prop in ['max', 'min', 'mean', 'p10', 'p25', 'p50', 'p75', 'p90']:
            vlist.append( '%s_%s' % (vname, prop) )
    # ================================================================


    # loop over ftype combinations -----------------------------------
    d = read_cluster_props( vlist, flist = flist)

    d['time_id'] = create_time_id(d['abs_time'], d['rel_time'])

    remove_too_few_clusters(d)
    
    d['aspect'] = d['diameter']**2 / d['hull_dmax']**2


    props = dict(diameter = d['diameter'])
    dist, nn_props = calculate_nn_distance(d['x_mean'],
                                           d['y_mean'], 
                                           d['time_id'],
                                           props = props)
    
    d['nN_distance'] = dist[:, 0]
    d['2nd_Neighbor_distance'] = dist[:, 1]
    d['nN_diameter'] = nn_props['diameter']

    R = 0.5 * d['diameter']
    d['weighted_nN_distance'] = dist[:, 0] / R
    d['weighted_2nd_Neighbor_distance'] = dist[:, 1] / R
 
    if 'imf_mean' in d.keys():
        d['mass_lift'] = d['imf_mean'] * d['area'] * 1e6 

   

    out = d
    # ================================================================



    # save stuff into hdf --------------------------------------------
    if do_output:
        if output_file is None:
            oname = 'collected_cluster_props.h5'
        else:
            oname = output_file

        print '... save data to %s' % oname
        hio.save_dict2hdf(oname, out)
    # ================================================================

    return out



######################################################################
######################################################################
