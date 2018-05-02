#!/usr/bin/env python

import os, sys, glob, copy
import numpy  as np
import datetime
import scipy.ndimage, scipy.spatial

import io_tools.hdf as hio
import analysis_tools.grid_and_interpolation as gi
import analysis_tools.segmentation as seg
import mahotas

from standard_config import *
import netCDF4

sys.path.append('../../2016-05_convective_orga/tools/')
sys.path.append('../../2016-05_convective_orga/inout/')

import cluster_analysis 
import data_reader
        
sys.path.append('../tools')
from special_threshold_calculations import special_threshold_calculations



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
        

    # set the segmentation parameters by experient name
    try:
        expname = sys.argv[2]
    except:
        expname = 'basic'

    # set the variable name
    try:
        varname = sys.argv[3]
    except:
        varname = 'bt108'
    
    # ALSO ALLOW FOR AUX DATA
    do_aux_input = True
    aux_vname = 'tcw'

    sett = predefined_setups(expname.split('-')[0])
    setup_for_later_output = copy.copy(sett)

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

    # THIS IS THE MAIN CHANGE IN THE PROGRAM!!!
    din = data_reader.read_narval_addvars(fname, varname)

    output_dir = din['input_dir']

    dset = din.copy()

    if do_aux_input:
        aux_fname = fname.replace(varname, aux_vname)
        daux = data_reader.read_narval_addvars(aux_fname, aux_vname)
    # ================================================================




    # cluster properties ---------------------------------------------
    # get field dimensions

    # check if field should be inverted
    if len(expname.split('-')) == 2 and expname.split('-')[1] == 'inv':
        b3d = -din[varname]
    else:
        b3d = din[varname]


    ntime, nrow, ncol = b3d.shape

    b_segmented = np.zeros_like( b3d )
    noffset = 0
    cset = {}

    for itime in range(ntime):
        
        print '... calculate cluster properties for %d' % itime

        # get local time field
        b = b3d[itime]

        if thresh == 'relative50_for_mass_flux':
            thresh = special_threshold_calculations(din['lon'], din['lat'], b, 
                                                    method = 'relative50_for_mass_flux')

        # set masks
        b.data[b.mask] = thresh / 2.
        
        # get cluster field
        if  USE_EXISTING_CLUSTER_DATA:
            c = c3d[itime]
            c = seg.remove_small_clusters(c, min_size = min_size)

        else:
            c = seg.clustering(b, thresh, **sett)


        # remove edge connections
        dist_to_edge = 10.
        c = mahotas.labeled.remove_bordering(c, dist_to_edge)

        b_segmented[itime] = c
        # c = seg.remove_small_clusters(c, min_size = min_size)
        
        # get distance field
        dist = scipy.ndimage.distance_transform_edt(c)


        dset['clust'] = c
        dset[varname] = b
        dset['%s_trans' % varname] = np.where(b < thresh, 0, b - thresh)

        dset['dist'] = dist
        dset['rel_time'] = din['rel_time'][itime]
        dset['abs_time'] = din['abs_time'][itime]

        if do_aux_input:
            dset[aux_vname] = daux[aux_vname][itime]
            var_names =  [varname, 'dist', aux_vname]
        else:
            var_names =  [varname, 'dist']

        cluster_analysis.cluster_analysis(dset, cset, 
                                          noffset = noffset, 
                                          var_names = var_names, 
                                          weight_names = ['%s_trans' % varname, 'dist'],) 
        noffset += c.max()

    for cprop in cset.keys():
        cset[cprop] = np.array( cset[cprop] )
    # ================================================================


    # # output ---------------------------------------------------------
    out = cset

    if not USE_EXISTING_CLUSTER_DATA:
        out['_settings'] = setup_for_later_output
    
    
    odir = '%s/cluster_properties' % output_dir
    oname = '%s/clust_prop_%s_%s.h5' % (odir, basename, expname)
    print '... save output to %s' % oname
    hio.save_dict2hdf(oname, out)


    # and segmented field ..............................................
    if  not USE_EXISTING_CLUSTER_DATA:
        out = {}
        out['b_segmented'] = b_segmented.astype(np.uint16)
        out['_settings'] = setup_for_later_output
        
        
        odir = '%s/cluster_properties' % output_dir
        oname = '%s/segmented_%s_%s.h5' % (odir, basename, expname)
        print '... save output to %s' % oname
        hio.save_dict2hdf(oname, out)
    # # ================================================================


