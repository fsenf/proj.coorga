#!/usr/bin/env python

import os, sys, glob
import numpy  as np
import datetime
import scipy.ndimage, scipy.spatial

import io_tools.hdf as hio
import analysis_tools.grid_and_interpolation as gi
import analysis_tools.segmentation as seg
import mahotas

from standard_config import *
import data_reader
import netCDF4

sys.path.append('../tools/')

import cluster_analysis 
        


######################################################################
######################################################################



def predefined_setups(name):


    # BASIS SETUP ----------------------------------------------------
    if name == 'basic':
        return dict( name = name,
                     cluster_method = 'watershed_merge', 
                     min_size = 40,
                     base_field = 'BT108', 
                     thresh = 230,
                     dradius = 3,
                     ctype = 4,
                     marker_field = 'dist',
                     filter_method = 'curve',
                     marker_method = 'iterative_shrinking',
                     numberOfIterations = 5,
                     cluster_masking = True,
                     siggauss = 1.,
                     exclude_border = True
                 )
    # ================================================================


    elif name == 'exp001':
        b = predefined_setups('basic')
        b.update(dict( cluster_method = 'watershed' ))
        b.update(dict( dradius = 10 ))
        b.update(dict( name = name ))
        b.update(dict( marker_field = 'field'))
        return b


    elif name == 'exp002':
        b = predefined_setups('basic')
        b.update(dict( cluster_method = 'connect' ))
        b.update(dict( name = name ))
        return b

    elif name == 'exp003':
        b = predefined_setups('basic')
        b.update(dict( thresh = 210 ) )
        b.update(dict( name = name ))
        return b

    elif name == 'exp004':
        b = predefined_setups('exp002')
        b.update(dict( thresh = 210 ) )
        b.update(dict( name = name ))
        return b

    elif name == 'exp005':
        b = predefined_setups('exp002')
        b.update(dict( min_size = 1 ))
        b.update(dict( name = name ))
        return b

    elif name == 'basicde':
        b = predefined_setups('basic')
        b.update( dict( thresh = 240 ) )
        return b

    return 


######################################################################
######################################################################

def read_satdata(fname):
    
    if 'narval' in fname:
        return data_reader.read_narval_data(fname)

    elif 'hdfd' in fname:
        return data_reader.read_hdcp2_data(fname)

    elif 'icon-lem' in fname:
        return data_reader.read_icon_lem_data(fname)



######################################################################
######################################################################

def check_if_clusterfile(fname):

    f = netCDF4.Dataset(fname, 'r')
    v = f.variables
    vnames = v.keys()
    f.close()

    if 'iclust' in vnames:
        return True
    else:
        return False


######################################################################
######################################################################

def data_fname_from_cluster_fname(clustname, fileext = 'h5'):

    basename = os.path.splitext(os.path.basename(clustname))[0]
    narval_dir = '%s/icon/narval' % local_data_path
 
    if 'narval' in basename:
        baselist = basename.split('_')
        inarv = baselist.index('narval')

        newbase = '_'.join( baselist[inarv - 1 : inarv + 3] )
        clustersuff = '_'.join( baselist[inarv + 3 : ] )
        newname = '%s/%s.%s' % (narval_dir, newbase, fileext)

    return newname, clustersuff
    

######################################################################
######################################################################




if __name__ == '__main__':


    # get input ------------------------------------------------------
    fname = sys.argv[1]
    basename = os.path.splitext(os.path.basename(fname))[0]
        
    try:
        expname = sys.argv[2]
    except:
        expname = 'basic'


    sett = predefined_setups(expname)
    thresh = sett.pop('thresh')
    min_size = sett.get('min_size')
    # ================================================================


    # ----------------------------------------------------------------
    # two chances: 
    # (i) fname is standard sat.-data input filename
    #       -> then, everything works like usual
    #
    # (ii) fname is the cluster id filename 
    #       -> then, sat filename has to be constructed
    #       -> cluster data have to be read

    if check_if_clusterfile(fname):
        USE_EXISTING_CLUSTER_DATA = True

        # get data filename
        clustname = fname[:]
        fname, expname = data_fname_from_cluster_fname(clustname, fileext = 'nc') 

        print '... read cluster data'
        c3d = data_reader.read_cluster(clustname)
    else:
        USE_EXISTING_CLUSTER_DATA = False 

    basename = os.path.splitext(os.path.basename(fname))[0]
    # ================================================================


    # read the satellite data ----------------------------------------
    din = read_satdata(fname)
    output_dir = din['input_dir']

    dset = din.copy()
    # ================================================================




    # cluster properties ---------------------------------------------
    # get field dimensions
    b3d = din['bt108']
    ntime, nrow, ncol = b3d.shape

    b_segmented = np.zeros_like( b3d )
    noffset = 0
    cset = {}

    for itime in range(ntime):
        
        print '... calculate cluster properties for %d' % itime

        # get local time field
        b = b3d[itime]

        # set masks
        b.data[b.mask] = thresh / 2.
        
        # get cluster field
        if  USE_EXISTING_CLUSTER_DATA:
            c = c3d[itime]
            c = seg.remove_small_clusters(c, min_size = min_size)

        else:
            c = seg.clustering(-b, -thresh, **sett)


        # remove edge connections
        dist_to_edge = 10.
        c = mahotas.labeled.remove_bordering(c, dist_to_edge)

        b_segmented[itime] = c
        # c = seg.remove_small_clusters(c, min_size = min_size)
        
        # get distance field
        dist = scipy.ndimage.distance_transform_edt(c)


        dset['clust'] = c
        dset['bt108'] = b

        dset['bt108_trans'] = np.where(b < thresh, thresh - b, 0)
        dset['dist'] = dist
        dset['rel_time'] = din['rel_time'][itime]
        dset['abs_time'] = din['abs_time'][itime]

        

        cluster_analysis.cluster_analysis(dset, cset, noffset = noffset) 
        noffset += c.max()

    for cprop in cset.keys():
        cset[cprop] = np.array( cset[cprop] )
    # ================================================================


    # # output ---------------------------------------------------------
    out = cset

    if not USE_EXISTING_CLUSTER_DATA:
        out['_settings'] = predefined_setups(expname)
    
    
    odir = '%s/cluster_properties' % output_dir
    oname = '%s/clust_prop_%s_%s.h5' % (odir, basename, expname)
    print '... save output to %s' % oname
    hio.save_dict2hdf(oname, out)


    # and segmented field ..............................................
    if  not USE_EXISTING_CLUSTER_DATA:
        out = {}
        out['b_segmented'] = b_segmented.astype(np.uint16)
        out['_settings'] = predefined_setups(expname)
        
        
        odir = '%s/cluster_properties' % output_dir
        oname = '%s/segmented_%s_%s.h5' % (odir, basename, expname)
        print '... save output to %s' % oname
        hio.save_dict2hdf(oname, out)
    # # ================================================================


