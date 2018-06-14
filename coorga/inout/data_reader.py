#!/usr/bin/env python

import os, sys, glob
import numpy  as np
import datetime
import netCDF4

import tropy.io_tools.hdf as hio
import tropy.io_tools.netcdf as ncio
import tropy.analysis_tools.grid_and_interpolation as gi
from tropy.standard_config import *

######################################################################
######################################################################

def input(fname, varname, **input_options):
    
    '''
    Input Interface as data reader wrapper.

    
    Parameters
    ----------
    fname : str 
       name of data file

    varname: str
       name of variable (if fname is cluster file than 
       varname should be the correspsong basic field)


    Returns
    --------
    dset : dict
       dictionary that collects a set of input data
    '''



    # ----------------------------------------------------------------
    # two chances: 
    # (i) fname is standard basic data input filename
    #       -> then, everything works like usual
    #
    # (ii) fname is the cluster id filename 
    #       -> then, sat filename has to be constructed
    #       -> cluster data have to be read


    if check_if_clusterfile(fname):

        # get data filename
        clustname = fname[:]
        fname, expname = data_fname_from_cluster_fname(
            clustname, fileext = 'nc') 

        print '... read cluster data'
        cset = dict( iclust = read_cluster( clustname ) )
        cset['clustname'] = clustname
        cset['fname'] = fname
    # ================================================================


    # read data ------------------------------------------------------
    if varname == 'bt108':
        if 'narval' in fname:
            dset = read_narval_data(fname)

        elif 'hdfd' in fname:
            dset = read_hdcp2_data(fname)

        elif 'icon-lem' in fname:
            dset = read_icon_lem_data(fname)


    else:
        if 'narval' in fname:
            dset = read_narval_addvars(fname, varname, **input_options)
    # ================================================================
            


    # combined data with cluster info if there ...
    try:
        dset.update( cset )
    except:
        # no cluster file ...
        pass

    return dset


######################################################################
######################################################################

def check_if_clusterfile(fname):

    '''
    Checks if file is netcdf and contains a cluster field.

    
    Parameters
    ----------
    fname : str
       filename


    Returns
    --------
    decision : bool
       bool variable and says if file is clusterfile
    
    '''

    f = netCDF4.Dataset(fname, 'r')
    v = f.variables
    vnames = v.keys()
    f.close()

    if 'iclust' in vnames or 'b_segmented' in vnames:
        return True
    else:
        return False


######################################################################
######################################################################


def data_fname_from_cluster_fname(clustname, fileext = 'h5'):


    '''
    Connects the data filename to the cluster filename.

    
    Parameters
    ----------
    clustname : str
       cluster file name


    fileext : str, optional, default = 'h5'
       file extention of date file


    Returns
    --------
    newname : str
       name of the data file that belongs to the clusterfile
      
    expname : str
       name of the segmentation setup

    '''

    # directory reconstruction
    # ========================
    dirname = os.path.dirname( clustname )

    subdirectories = dirname.split( '/' )

    # remove the last subdirectory
    newdir = '/'.join( subdirectories[:-1] )
 



    # filename reconstruction
    # =======================

    # the standard filename contains 'segmented' in the first place
    #                           and the expname in the last place
    
    basename = os.path.splitext(os.path.basename(clustname))[0]
    baselist = basename.split('_')

    expname = baselist[-1]
    newbase = '_'.join( baselist[1:-1] )

    newname = '/%s/%s.%s' % ( newdir, newbase, fileext )
    
    return newname, expname


######################################################################
######################################################################


def read_narval_data(fname):

    '''
    Reads the time stack of Narval data, either meteoat or synsat.

    
    Parameters
    ----------
    fname : str
       filename of data file
    
    
    Returns
    --------
    dset : dict
       dataset dictionary containing georef and bt108 data.
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
    index_time = np.arange(ntime)

    day_shift = t0 - datetime.datetime(1970, 1,1)
    day_shift = day_shift.total_seconds() / (24.*3600)

    abs_time = day_shift + rel_time / 24.
    # ================================================================
    
    
    # get georef .....................................................
    gfile = '%s/aux/target_grid_geo_reference_narval.h5' % narval_dir
    geo = hio.read_dict_from_hdf(gfile)
    lon, lat = geo['lon'], geo['lat']

    # centered sinusoidal
    x, y = gi.ll2xyc( lon, lat )
    area = np.abs( gi.simple_pixel_area(lon, lat) )
    # ================================================================



    # prepare output .................................................
    dset = {}
    vnames = ['x', 'y', 'lon', 'lat', 'lsm', 'area', 'rel_time', 'abs_time', 'index_time']
    vvec   = [ x ,  y ,  lon ,  lat ,  lsm ,  area ,  rel_time ,  abs_time ,  index_time]
    for i, vname in enumerate(vnames):
        dset[vname] = vvec[i]

    dset['bt108'] = b3d
    dset['input_dir'] = narval_dir
    # ================================================================


    return dset

######################################################################
######################################################################

def generic_addvar_reader(fname, vname):

    '''
    Reads the time stack of additional variables.


    Parameters
    ----------
    fname : str 
        filename of data file
    
    vname : str
        variable name 
        (variable should be contained in file)


    Returns
    --------
    dset : dict
        dataset dictionary containing georef and add data.
    '''

    dset = ncio.read_icon_4d_data(fname, [vname,], itime = None)

    geo = ncio.read_icon_georef(fname)
    geo['x'], geo['y'] = gi.ll2xyc( geo['lon'], geo['lat'] )

    dset.update( geo )

    return dset
   
######################################################################
######################################################################



def read_narval_addvars(fname, vname, 
                            domain_center = None,
                            region_slice = None):

    '''
    Reads the time stack of Narval data, either meteoat or synsat.


    Parameters
    ----------
    fname : str 
        filename of data file
    
    vname : str
        variable name 
        (variable should be contained in file)

    domain_center : tuple of floats, optional, default = None
        setting the projection center to (clon, clat)
        if None: not used

    region_slice : tuple of floats, optional, default = None
        cutout of fields for form  ((irow1, irow2), (icol1, icol2))
        if None: not used


    Returns
    --------
    dset : dict
        dataset dictionary containing georef and bt108 data.
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

    b3d = ncio.read_icon_4d_data(fname, [vname,], itime = None)[vname]
    b3d = np.ma.masked_invalid( b3d )

    ntime, nrow, ncol = b3d.shape
    # ================================================================


    # prepare time vector --------------------------------------------
    rel_time = np.arange(1, ntime + 1)
    index_time = np.arange(ntime)



    day_shift = t0 - datetime.datetime(1970, 1,1)
    day_shift = day_shift.total_seconds() / (24.*3600)

    abs_time = day_shift + rel_time / 24.
    # ================================================================
    
    
    # get georef .....................................................
    gfile = '%s/aux/target_grid_geo_reference_narval.h5' % narval_dir
    geo = hio.read_dict_from_hdf(gfile)
    lon, lat = geo['lon'], geo['lat']

    if domain_center is not None:
        mlon, mlat = domain_center
    else:
        mlon, mlat = None, None

    x,y = gi.ll2xyc(lon, lat, mlon = mlon, mlat = mlat)

    area = np.abs( gi.simple_pixel_area(lon, lat) )
    # ================================================================



    # prepare output .................................................
    dset = {}
    dset[vname] = b3d

    addnames = ['x', 'y', 'lon', 'lat', 'lsm', 'area', 'rel_time', 'abs_time', 'index_time']
    vvec     = [ x ,  y ,  lon ,  lat ,  lsm ,  area ,  rel_time ,  abs_time ,  index_time]
    for i, aname in enumerate(addnames):
        dset[aname] = vvec[i]

    dset['input_dir'] = os.path.dirname(fname)
    # ================================================================


    # do cutout if wanted --------------------------------------------
    field_names =  [ 'x', 'y', 'lon', 'lat', 'lsm', 'area' , vname ]

    if region_slice is not None:
        for name in field_names:
            dset[name] = gi.cutout_fields( dset[name], region_slice, 
                                           vaxis = 0)
    # ================================================================

    return dset

######################################################################
######################################################################



def read_hdcp2_data(fname):


    '''
    Reads BT10.8 data from files generated by the HDCP2 O module.


    Parameters
    ----------
    fname : str
       file name


    Returns
    --------
    dset : dict
        dictionary of datasets
    '''


    # data fields ----------------------------------------------------
    vlist = ['tb108', 'lon', 'lat', 'time']
    dset = ncio.read_icon_4d_data(fname, vlist,
                                    itime = None)

    b3d = dset.pop('tb108')
    b3d = np.ma.masked_less(b3d, 100.)
    # ================================================================

    
    # geo ref --------------------------------------------------------
    lon, lat = dset['lon'], dset['lat']
    
    x,y = gi.ll2xyc(lon, lat, lon0 = 10, lat0 = 50)
    area = np.abs( gi.simple_pixel_area(lon, lat) )
    # ================================================================


    # time conversions -----------------------------------------------
    abs_time = dset['time'] / (3600. * 24 )
    rel_time = np.mod(abs_time, 1) * 24.

    ntime = len(rel_time)
    index_time = np.arange(ntime)
    # ================================================================
    

    # prepare output .................................................
    vnames = ['x', 'y', 'area', 'rel_time', 'abs_time', 'index_time']
    vvec   = [ x ,  y ,  area ,  rel_time ,  abs_time ,  index_time]
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

    '''
    Reads BT10.8 data from files generated for ICON-LEM runs.


    Parameters
    ----------
    fname : str
       file name


    Returns
    --------
    dset : dict
        dictionary of datasets
    '''


    # data fields ----------------------------------------------------
    vlist = ['bt108', 'lon', 'lat', 'time']
    dset = ncio.read_icon_4d_data(fname, vlist,
                                    itime = None)

    b3d = dset.pop('bt108')
    b3d = np.ma.masked_less(b3d, 100.)
    # ================================================================

    
    # geo ref --------------------------------------------------------
    lon, lat = dset['lon'], dset['lat']
    
    x,y = gi.ll2xyc(lon, lat, lon0 = 10, lat0 = 50)
    area = np.abs( gi.simple_pixel_area(lon, lat) )
    # ================================================================


    # time conversions -----------------------------------------------
    rel_time = 24 * (dset['time'] - dset['time'][0])

    ntime = len(rel_time)
    index_time = np.arange(ntime)
 
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
    vnames = ['x', 'y', 'area', 'rel_time', 'abs_time', 'index_time']
    vvec   = [ x ,  y ,  area ,  rel_time ,  abs_time ,  index_time]
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

    '''
    Reads Cluster Data from netcdf files.


    Parameters
    ----------
    fname : str
       file name

    cname : str, optional, default = 'iclust'
       name of cluster index variable


    Returns
    --------
    c : numpy array, int
       cluster index field
    
    '''


    c = ncio.read_icon_4d_data(fname, cname, itime = None)[cname]

    return c.astype(np.int)


######################################################################
######################################################################



def read_vstack_from_hdflist( flist, vname,
                              itime0 = None,
                              nsub = None,
                              subpath = None):



    '''
    Reads a variable stack from hdf list.


    Parameters
    ----------
    flist : list
       list of hdf file names

    vname : str
       variable name 

    itime0 : int, optional, default = None
       first time index to count out data temporally

    nsub : int, optional, default = None
       subsampling factor

    subpath : str, optional, default = None
       subpath for grouped variables, e.g. /main_group/sub_group1/...'
    
    

    Returns
    --------
    v : numpy array
        stack of variable
    '''


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


    fname = '%s/icon/lem-de/synsat/msevi_3d_coarse_icon-lem-tstack_DOM_20150704-shifted-double_day.nc' % local_data_path
    dset = read_icon_lem_data(fname)


