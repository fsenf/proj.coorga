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


def read_cluster_props( vlist,
                        ftype = 'msevi',
                        expname = 'basic',
                        filepart = '_narval_DOM01_',
                        fdir = '%s/cluster_properties' % narval_dir,
                        subpath = None,
                        time_masking = True,
                        date = '', 
                        as_array = True):

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


if __name__ == '__main__':

    cprop = read_cluster_props(['x_mean', 'y_mean', 'fraction_of_lsm_types'], 
                               date = '201608')
