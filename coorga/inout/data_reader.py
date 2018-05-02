#!/usr/bin/env python

import os, sys, glob
import numpy  as np
import datetime

import io_tools.hdf as hio
import io_tools.netcdf as ncio

import analysis_tools.grid_and_interpolation as gi

from standard_config import *


######################################################################
######################################################################


def read_narval_data(fname):

    '''
    Reads the time stack of Narval data, either meteoat or synsat.

    
    USAGE
    ======
    dset = read_narval_data(fname)


    INPUT
    =====
    fname: filename of data file
    
    
    OUTPUT
    ======
    dset: dataset dictionary containing georef and bt108 data.
    '''

    # read land sea data ---------------------------------------------
    narval_dir = '%s/icon/narval' % local_data_path
    lsm_name = '%s/aux/narval_landsea_coast_mask.h5' % narval_dir
    
    print '... read land-sea-mask from %s' % lsm_name
    dset = hio.read_dict_from_hdf(lsm_name)

    lsm = dset['mask50']    
    # ================================================================
    

    
    # read bt108 -----------------------------------------------------
    print '... read BT10.8 from %s' % fname
    basename, file_ext = os.path.splitext(os.path.basename(fname))

    date = basename.split('_')[-1]
    t0 = datetime.datetime.strptime(date, '%Y%m%d')

    # check if its is obs or sim?
    ftype = basename.split('_')[0]
    if ftype in ['msevi', 'trans']:
        subpath = None
    elif ftype == 'synsat':
        subpath = 'synsat_oper'
        
    # read bt108 from hdf
    if file_ext == '.h5':
        b3d = hio.read_var_from_hdf(fname, 'IR_108', subpath = subpath) / 100. 
    elif file_ext == '.nc':
        vname = 'bt108'
        b3d = ncio.read_icon_4d_data(fname, vname, itime = None)[vname]
        
    b3d = np.ma.masked_invalid( b3d )
    b3d = np.ma.masked_less(b3d, 100.)

    ntime, nrow, ncol = b3d.shape
    # ================================================================


    # prepare time vector --------------------------------------------
    rel_time = np.arange(1, ntime + 1)


    day_shift = t0 - datetime.datetime(1970, 1,1)
    day_shift = day_shift.total_seconds() / (24.*3600)

    abs_time = day_shift + rel_time / 24.
    # ================================================================
    
    
    # get georef .....................................................
    gfile = '%s/aux/target_grid_geo_reference_narval.h5' % narval_dir
    geo = hio.read_dict_from_hdf(gfile)
    lon, lat = geo['lon'], geo['lat']

    x,y = gi.ll2xy(lon, lat)
    area = np.abs( gi.simple_pixel_area(lon, lat) )
    # ================================================================



    # prepare output .................................................
    dset = {}
    vnames = ['x', 'y', 'lon', 'lat', 'lsm', 'area', 'rel_time', 'abs_time']
    vvec   = [ x ,  y ,  lon ,  lat ,  lsm ,  area ,  rel_time ,  abs_time ]
    for i, vname in enumerate(vnames):
        dset[vname] = vvec[i]

    dset['bt108'] = b3d
    dset['input_dir'] = narval_dir
    # ================================================================


    return dset

######################################################################
######################################################################


def read_narval_addvars(fname, vname):

    '''
    Reads the time stack of Narval data, either meteoat or synsat.

    
    USAGE
    ======
    dset = read_narval_data(fname)


    INPUT
    =====
    fname: filename of data file
    
    
    OUTPUT
    ======
    dset: dataset dictionary containing georef and bt108 data.
    '''

    # read land sea data ---------------------------------------------
    narval_dir = '%s/icon/narval' % local_data_path
    lsm_name = '%s/aux/narval_landsea_coast_mask.h5' % narval_dir
    
    print '... read land-sea-mask from %s' % lsm_name
    dset = hio.read_dict_from_hdf(lsm_name)

    lsm = dset['mask50']    
    # ================================================================
    

    
    # read bt108 -----------------------------------------------------
    print '... read %s from %s' % (vname, fname)
    basename, file_ext = os.path.splitext(os.path.basename(fname))

    date = basename.split('_')[-1]
    t0 = datetime.datetime.strptime(date, '%Y%m%d')

    b3d = ncio.read_icon_4d_data(fname, vname, itime = None)[vname]
    b3d = np.ma.masked_invalid( b3d )

    ntime, nrow, ncol = b3d.shape
    # ================================================================


    # prepare time vector --------------------------------------------
    rel_time = np.arange(1, ntime + 1)


    day_shift = t0 - datetime.datetime(1970, 1,1)
    day_shift = day_shift.total_seconds() / (24.*3600)

    abs_time = day_shift + rel_time / 24.
    # ================================================================
    
    
    # get georef .....................................................
    gfile = '%s/aux/target_grid_geo_reference_narval.h5' % narval_dir
    geo = hio.read_dict_from_hdf(gfile)
    lon, lat = geo['lon'], geo['lat']

    x,y = gi.ll2xy(lon, lat)
    area = np.abs( gi.simple_pixel_area(lon, lat) )
    # ================================================================



    # prepare output .................................................
    dset = {}
    dset[vname] = b3d

    addnames = ['x', 'y', 'lon', 'lat', 'lsm', 'area', 'rel_time', 'abs_time']
    vvec   = [ x ,  y ,  lon ,  lat ,  lsm ,  area ,  rel_time ,  abs_time ]
    for i, aname in enumerate(addnames):
        dset[aname] = vvec[i]

    dset['input_dir'] = os.path.dirname(fname)
    # ================================================================


    return dset

######################################################################
######################################################################



def read_hdcp2_data(fname):


    # data fields ----------------------------------------------------
    vlist = ['tb108', 'lon', 'lat', 'time']
    dset = ncio.read_icon_4d_data(fname, vlist,
                                    itime = None)

    b3d = dset.pop('tb108')
    b3d = np.ma.masked_less(b3d, 100.)
    # ================================================================

    
    # geo ref --------------------------------------------------------
    lon, lat = dset['lon'], dset['lat']
    
    x,y = gi.ll2xy(lon, lat)
    area = np.abs( gi.simple_pixel_area(lon, lat) )
    # ================================================================


    # time conversions -----------------------------------------------
    abs_time = dset['time'] / (3600. * 24 )
    rel_time = np.mod(abs_time, 1) * 24.
    # ================================================================
    

    # prepare output .................................................
    vnames = ['x', 'y', 'area', 'rel_time', 'abs_time']
    vvec   = [ x ,  y ,  area ,  rel_time ,  abs_time ]
    for i, vname in enumerate(vnames):
        dset[vname] = vvec[i]

    dset['bt108'] = b3d
    dset['lsm'] = np.ones_like(x)
    dset['input_dir'] = os.path.dirname( fname )
    # ================================================================




    return dset


######################################################################
######################################################################

def read_icon_lem_data(fname):


    # data fields ----------------------------------------------------
    vlist = ['bt108', 'lon', 'lat', 'time']
    dset = ncio.read_icon_4d_data(fname, vlist,
                                    itime = None)

    b3d = dset.pop('bt108')
    b3d = np.ma.masked_less(b3d, 100.)
    # ================================================================

    
    # geo ref --------------------------------------------------------
    lon, lat = dset['lon'], dset['lat']
    
    x,y = gi.ll2xy(lon, lat)
    area = np.abs( gi.simple_pixel_area(lon, lat) )
    # ================================================================


    # time conversions -----------------------------------------------
    rel_time = 24 * (dset['time'] - dset['time'][0])

    t0 = datetime.datetime(1970, 1, 1)
    abs_time = []
    for t in dset['time']:
        day = str( int(t) )
        subday = np.mod(t, 1)

        tobj = datetime.datetime.strptime(day, '%Y%m%d')
        tobj += datetime.timedelta(days = subday)

        dt = ( tobj - t0).total_seconds()
        abs_time.append( dt / (24. * 3600.) )

    abs_time = np.array( abs_time )

    # ================================================================
    

    # prepare output .................................................
    vnames = ['x', 'y', 'area', 'rel_time', 'abs_time']
    vvec   = [ x ,  y ,  area ,  rel_time ,  abs_time ]
    for i, vname in enumerate(vnames):
        dset[vname] = vvec[i]

    dset['bt108'] = b3d
    dset['lsm'] = np.ones_like(x)
    dset['input_dir'] = os.path.dirname( fname )
    # ================================================================




    return dset


######################################################################
######################################################################


def read_cluster(fname, cname = 'iclust'):

    c = ncio.read_icon_4d_data(fname, cname, itime = None)[cname]

    return c.astype(np.int)


######################################################################
######################################################################



def read_vstack_from_hdflist( flist, vname,
                              itime0 = None,
                              nsub = None,
                              subpath = None):




    v = []
    for f in flist:
        
        try:
            vs = hio.read_var_from_hdf(f, vname, subpath = subpath)
            
            # do subsampling ?
            if nsub != None:
                vs = vs[:, ::nsub, ::nsub]
                
            if itime0 != None:
                vs = vs[itime0:]

            
            v.append( vs )
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            print 'ERROR for file %s' % f

    v = np.vstack( v )

    return v

######################################################################
######################################################################



if __name__ == '__main__':

    fname = '%s/hdcp2/hdfd_trop_seviri00_l1_tb108_v00_20150704000000.nc' % local_data_path
    dset = read_hdcp2_data(fname)
    print dset['rel_time'], dset['abs_time']


    fname = '%s/icon/lem-de/msevi_3d_coarse_icon-lem-tstack_DOM_20150704.nc' % local_data_path
    dset = read_icon_lem_data(fname)


