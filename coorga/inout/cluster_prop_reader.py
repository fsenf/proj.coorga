#!/usr/bin/env python

# Copyright (c) 2015 - 2016

# TROPOS,
# Permoserstr. 15
# 04318 Leipzig, Germany. 

# Author:
# ====== 
# Fabian Senf <senf@tropos.de>


# This program is free software; you can redistribute it and/or 
# modify it under the terms of the GNU General Public License as 
# published by the Free Software Foundation; either version 3 of 
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details.

'''
Input of cluster data.
'''


import sys, os, glob
import datetime
import numpy as np

import tropy.io_tools.hdf as hio
from tropy.standard_config import *

narval_dir = '%s/icon/narval' % local_data_path

######################################################################
######################################################################

def params_from_clustname(clustfile):

    '''
    Reads parameters from standardized cluster filename.

    Assumes that date is in the filename and the only int part and the following
    format: *{vname}_{date}_{expname}*


    Parameters
    -----------
    clusterfile : str
        clusterfile name


    Returns
    --------
    vname : str
        variable name

    date : int
        date 

    expname : str
        name of segmentation setup
    '''

    # dissect into pieces

    basename = os.path.splitext( os.path.basename( clustfile) )[0]
    name_vec = basename.split('_')

    # get the date
    for i, entry in enumerate( name_vec ):
        try:
            date = int(entry)
            idate = i
        except:
            pass
    
    # get the others
    vname = name_vec[idate - 1]
    expname = name_vec[idate + 1]

    return vname, date, expname

######################################################################
######################################################################



def read_cluster_props( vlist,
                        ftype = 'msevi',
                        expname = 'basic',
                        filepart = '_narval_DOM01_',
                        fdir = '%s/cluster_properties' % narval_dir,
                        subpath = None,
                        time_masking = True,
                        date = '', 
                        as_array = True):

    '''
    Reads cluster properties from a stack of files.

    
    Parameters
    ----------
    vlist : list
        list of variables (cell properties)
    
    
    ftype : str, optional, default = 'msevi'
        file type (or mode) that is used for input, e.g. 'msevi' vs. 'synsat'

    expname : str, optional, default = 'basic'
        name of the segmentation setup

    filepart : str, optional, default = '_narval_DOM01_'
        part of the filenam edirectly after the ftype string

    fdir : str, optional, default = '%s/cluster_properties' % narval_dir
        file directory
    
    subpath : str, optional, default = None
        subpath for grouped variables, e.g. /main_group/sub_group1/...'
     
  
    time_masking  : str, optional, default = True
        switch if masking based on time variable should be done

    date : str, optional, default = ''
        date or part of date string

    as_array : str, optional, default = True
        if return is provided as numpy array



    Returns
    --------
    dset : dict
       dictionary of cell properties
    '''



    if type(vlist) == type(''):
        vlist = [vlist,]


    fname = 'clust_prop_%s%s*%s*_%s.h5' % (ftype, filepart, date, expname)
    flist = glob.glob('%s/%s' % (fdir, fname))

    dset = {}
    for f in sorted(flist):

        print '... read data from %s' % f
        fname = os.path.basename(f)
        dset[fname] = {}

        for vname in vlist:

            try:
                dset[fname][vname] = hio.read_var_from_hdf(f, vname, subpath = subpath)
            except:
                print 'Error: %s not in %s' % (vname, f)

    dout = {}
    if as_array:
        for f in sorted(dset.keys()):

            if dset[f].has_key('rel_time') and time_masking:
                tr = dset[f]['rel_time']
                m = (tr < 120)
            else:
                vname = dset[f].keys()[0]
                m = np.ones_like(dset[f][vname]).astype(np.bool)
                
            for vname in dset[f].keys():
                
                try:
                    d = dset[f][vname][m]
                
                    try:
                        dout[vname].append( d )
                    except:
                        dout[vname] = [d,]

                except:
                    # if no masking applies e.g. for rbins
                    dout[vname] = dset[f][vname]


        for vname in dout.keys():

            if np.ndim( dout[vname][0] ) == 1:
                stack = np.hstack

            elif np.ndim( dout[vname][0] ) >= 2:
                stack = np.row_stack


            dout[vname] = stack( dout[vname] )

    else:
        dout = dset

    return dout

######################################################################
######################################################################

def read_collected_clusterdata(date = '201608',
                               expname = 'exp010',
                               vname = 'smf',
                               fdir = '%s/cluster_properties' % narval_dir,
                               bootstrap = False, 
                               bootstrap_extension_suffix = '_dext500'):
    
    '''
    Calculates average number density.
    

    Parameters
    ----------
    date : str, optional, default = '201608'
       date or part of date string

    expname : str, optional, default = 'exp010'
       name for segmentation setup

    vname : str, optional, default = 'smf'
       name of variable on which cell analysis is based

    fdir : str, optional, default = '%s/cluster_properties' % narval_dir
       file directory name

    bootstrap : bool, optional, default = False
       switch if bootrstrap file is used

    bootstrap_extension_suffix : str, optional, default = '_dext500'
       possible suffix of bootstrap file (could contain ID)

    
    Returns
    --------
    dset : dict of numpy arrays
       full set of cell properties
    '''
    

    # reading the cluster data

    if bootstrap:
        fname = '%s/bootstrap_collected_cluster_props_%s_%s_%s_%s.h5' % (fdir, 
                                                                         vname, 
                                                                         date, 
                                                                         expname, 
                                                                         bootstrap_extension_suffix
        )        
    else:
        fname = '%s/collected_cluster_props_%s_%s_%s.h5' % (fdir, 
                                                            vname, 
                                                            date, 
                                                            expname)

    print '...read data from ', fname

    dset = hio.read_dict_from_hdf(fname)

      
    return dset[vname]


#####################################################################
#####################################################################


if __name__ == '__main__':

    cprop = read_cluster_props(['x_mean', 'y_mean', 'fraction_of_lsm_types'], 
                               date = '201608')
