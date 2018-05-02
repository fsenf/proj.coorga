#!/usr/bin/env python

import sys, os, glob
import datetime
import numpy as np
import io_tools.hdf as hio
import analysis_tools.grid_and_interpolation as gi
import scipy.spatial

from standard_config import *


sys.path.append('%s/.python' % proj_path)

from  convective_orga.inout.cluster_prop_reader import read_cluster_props


######################################################################
######################################################################

def remove_to_less_clusters(dset):

    tname = 'time_id'

    t = dset[tname]
    tvec = set(t)

    # loop over different time instances
    for ti in sorted(tvec):
        
        mask = (dset[tname] == ti) 
        nonmask = np.logical_not(mask)

        if mask.sum() < 3:

            for vname in dset.keys():
                dset[vname] = dset[vname][nonmask]
    return
 
######################################################################
######################################################################



def calculate_nn_distance(x, y, t, k = 2, props = {}):
    
    '''
    The subroutine uses KDTree Algorithm to calculate 
    nearest-neighor distance for k neighbors.

    USAGE
    =====
    dist = calculate_nn_distance(x, y, t, k = 2)


    INPUT
    =====
    x: x-coordinate
    y: y-coordinate
    t: t-coordinate, which is used for sequential masking

    
    OUTPUT
    ======
    dist: nearest neighbor distances field
    '''

    tvec = set(t)
    
    dist = []

    nn_props = {}
    for pname in props.keys():
        nn_props[pname] = []

    # loop over different time instances
    for ti in sorted(tvec):
        
        mask = (t == ti)
        
        xy = np.column_stack([x[mask], y[mask]])
        kdtree =  scipy.spatial.KDTree(xy)
        
        d, index = kdtree.query(kdtree.data, k = k + 1)
        
        dist.append( d )

        for pname in props.keys():
            nn_props[pname].append( props[pname][mask][index[:,1]] )

    dist = np.row_stack(dist )

    for pname in props.keys():
        nn_props[pname] = np.hstack( nn_props[pname] )

    if len(nn_props.keys()) == 0:
        return dist[:,1:]
    else:
        return dist[:,1:], nn_props
        


######################################################################
######################################################################

def create_time_id(abs_time, rel_time):


    '''
    assume 
    * abs_time as array of days since epoche
    * rel_time as array of hours
    '''

    id_list = []
    
    for i, tr in enumerate(rel_time):

        init_day = abs_time[i] - rel_time[i]/24.
        init_time = datetime.datetime(1970,1,1) + datetime.timedelta(days = init_day)
        init_str = init_time.strftime('%Y%m%d_%H%M')
        
        hour = '%s' %  str(np.int(tr)).zfill(2)
        
        id_list.append( '%s_%s' % (init_str, hour) )

    return np.array( id_list )



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


    # for rain / water colum comparison ------------------------------
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

        remove_to_less_clusters(d)
        
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
    fdir = kws['fdir']
    oname = '%s/collected_cluster_props_%s_%s_%s.h5' % (fdir,  varname, date, expname)
    print '... save data to %s' % oname
    hio.save_dict2hdf(oname, out)
    # ================================================================
