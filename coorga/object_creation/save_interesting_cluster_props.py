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
 

######################################################################
######################################################################

def main(expname, regtype, varname, date, do_output = True):

    '''
    Gathers a collection of cell properties for a selected 
    time period, calculates some new props and finnaly saves results.


    Parameters
    ----------
    expname : string
        name of the experiment (segmentation setup)

    regtype : string
        region name

    varname : string
        variable name

    date : string, typically in format %Y%m
        string that specifies the date or part of the date

    do_output: bool, optional, default = True
        switch that determines if output is written in predefined file


    Returns
    -------
    out: dict
        dictionary that contains cell properties


    '''

    # select file parameters -----------------------------------------
    if regtype == 'narval' and varname == 'bt108':
        
        ftype_combinations = [  ('msevi', 'obs'),  
                                ('synsat', 'sim'),
                                ('trans', 'trans')]

        narval_dir = '%s/icon/narval/synsat/' % local_data_path

        kws = dict(
            filepart = '_narval_DOM01_',
            fdir = '%s/cluster_properties' % narval_dir,
            date = date)

        addlist = ['fraction_of_lsm_types', ]


    elif regtype == 'narval' and varname in ['smf', 'rr', 'imf']:
        
        ftype_combinations = [  ('icon-narval', varname) ] 


        kws = dict(
            filepart = '_dom01_%s_' % varname,
            fdir = '%s/icon/narval/variables/cluster_properties'  % local_data_path,
            date = date)

        addlist = ['fraction_of_lsm_types', ]


    elif regtype == 'hdcp2':

        ftype_combinations = [('hdfd', 'obs'),]

        

        kws = dict(
            filepart = '_trop_seviri',
            time_masking = False,
            fdir = '%s/hdcp2/cluster_properties' % local_data_path,
            date = date,
            )

        addlist = []


    elif regtype == 'icon-lem':

        ftype_combinations = [('msevi', 'msevi'),
                              ('synsat1', 'synsat1'),
                              ('synsat2', 'synsat2'),
                              ('synsat3', 'synsat3')]

        

        kws = dict(
            filepart = '_3d_coarse_icon-lem-tstack_DOM',
            time_masking = False,
            fdir = '%s/icon/lem-de/cluster_properties' % local_data_path,
            date = date,
            )

        addlist = []



    # ================================================================
        


    # variables list -------------------------------------------------
    vlist =  addlist + [
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
        'diameter_set']
    # ================================================================


    # for rain / water column comparison ------------------------------
    for vname in ['rr', 'tcw', 'imf']:
        for prop in ['max', 'min', 'mean', 'p10', 'p25', 'p50', 'p75', 'p90']:
            vlist.append( '%s_%s' % (vname, prop) )
    # ================================================================


    # loop over ftype combinations -----------------------------------
    out = {}
    
    for ftype, outname in ftype_combinations:
        d = read_cluster_props( vlist,
                                ftype = ftype,
                                expname = expname, 
                                **kws)


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
        

        out[outname] = d
    # ================================================================



    # save stuff into hdf --------------------------------------------
    if do_output:
        fdir = kws['fdir']
        oname = '%s/collected_cluster_props_%s_%s_%s.h5' % (fdir,  varname, date, expname)
        print '... save data to %s' % oname
        hio.save_dict2hdf(oname, out)
    # ================================================================

    return out



######################################################################
######################################################################

if __name__ == '__main__':

    # possible inputs ------------------------------------------------
    try:
        expname = sys.argv[1]
    except:
        expname = 'basic'


    try: 
        regtype = sys.argv[2]
    except:
        regtype = 'narval'


    try: 
        varname = sys.argv[3]
    except:
        varname = 'imf'


    try:
        date = sys.argv[4]
    except:
        date = '201608'
    # ================================================================


    main(expname, regtype, varname, date)

