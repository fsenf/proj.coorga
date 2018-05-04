#!/usr/bin/env python

import copy
import numpy  as np
import datetime
import scipy.spatial
import mahotas

# tropy libs
import tropy.analysis_tools.grid_and_interpolation as gi
import tropy.analysis_tools.segmentation as seg
from tropy.standard_config import *

# local libs
from segmentation_config import predefined_setups


######################################################################
######################################################################

def remove_too_few_clusters(dset):

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


def single_cell_analysis(dset, 
                  var_names = ['bt108', 'dist'],
                  weight_names = ['bt108_trans', 'dist'],
                  do_landsea_fraction = True):

    '''
    Calculates a set of properties for a masked cell.
    
    INPUT
    =====
    dset: input data set, including the cell mask

    
    OUTPUT
    ======
    cset: set of cell properties
    '''



    # get cell mask --------------------------------------------------
    mask = dset['mask']

    # and cell coordinates
    xc =  dset['x'][mask]
    yc =  dset['y'][mask]

    # area and land-sea-coast mask
    ac = dset['area'][mask]

    # and variables
    vc = {}
    for  vname in var_names:
        vc[vname] = dset[vname][mask]
    # ================================================================


    # fill the cell set with cell averages ---------------------------
    cset = {}

    cset['x_mean'] = xc.mean()
    cset['y_mean'] = yc.mean()
    
    cset['lon_mean'], cset['lat_mean'] = gi.xy2ll(cset['x_mean'], cset['y_mean'])

    for vname in var_names:
        cset['%s_mean' % vname] =   vc[vname].mean()

        cset['%s_min' % vname] =   vc[vname].min()
        cset['%s_max' % vname] =   vc[vname].max()
        
        for p in [10,25,50,75,90]:
            cset['%s_p%d' % (vname, p)] =   np.percentile(vc[vname], p)
    # ================================================================


    # times ----------------------------------------------------------
    cset['abs_time'] = dset['abs_time']
    cset['rel_time'] = dset['rel_time'] 
    cset['loc_time'] = dset['rel_time'] +  cset['lon_mean'] * 24. / 360.
    # ================================================================



    # integral cell values -------------------------------------------
    npx = len( xc )
    cset['npx'] = npx
    cset['area'] = ac.sum()
    cset['diameter'] = 2 * np.sqrt( cset['area'] / np.pi)
    # ================================================================


    # land sea coast type contributions ------------------------------
    if do_landsea_fraction:
        lsmc = dset['lsm'][mask]
        h, xe = np.histogram(lsmc, np.arange(0.5,6.5, 1))
        cset['fraction_of_lsm_types'] = h.astype(np.float) / npx
    # ================================================================


    # convex hull ----------------------------------------------------
    try:
        p = np.column_stack( [xc, yc] )
        hull = scipy.spatial.ConvexHull(p)
        xh, yh = xc[hull.vertices], yc[hull.vertices]

        # exclude convex hull to have equal shapes per cell  
        # cset['hull_xv'] = xc[hull.vertices]
        # cset['hull_yv'] = yc[hull.vertices]
        
        # get maximum distance in convex hull 
           
        # make a mesh for fast matrix-based distance calculation
        xx, yy = np.meshgrid( xh, yh )
     
        # get squared direction deviations
        dxq = (xx - xx.transpose())**2
        dyq = (yy - yy.transpose())**2
     
     
        dist_matrix = np.sqrt( dxq + dyq ) 
        
        cset['hull_dmax'] = dist_matrix.max()
    except:
        print 'ERROR in convex hull calculation'

        # exclude convex hull to have equal shapes per cell 
        # cset['hull_xv'] = [ 0 , ]
        # cset['hull_yv'] = [ 0 , ]
        cset['hull_dmax'] = 0.
    # ================================================================


    # calculate weighted averages ------------------------------------
    for wname in weight_names:
        
        weight_base = dset[wname][mask]
        
        w = weight_base / weight_base.sum()

        cset['x_weighted_ave_%s' % wname] = (w * xc).mean()
        cset['y_weighted_ave_%s' % wname] = (w * yc).mean()
 
    # ================================================================


    # volume and shape for different variables -----------------------
    for vname in var_names:
        
        vmin = cset['%s_min' % vname] 
        vmax = cset['%s_max' % vname] 
        vmean = cset['%s_mean' % vname] 

        shape = (vmean - vmin) / (vmax - vmin)
        vol  =  npx * shape

        cset['shape_%s' % vname]  =  shape
        cset['volume_%s' % vname] =  vol

    # ================================================================


    return cset
    
    

 
######################################################################
######################################################################


def cellset_analysis(dset, cset, noffset = 0, **kwargs):


    '''
    Performs subsequent single cell analysis for a set of cells.
    
    INPUT
    =====
    dset: input data set, including the cell mask

    
    OUTPUT
    ======
    cset: set of cell properties
    '''



    c = dset['clust']

    for nc in range(1, c.max() + 1):

        # get number of pixels in the cluster
        numpix =  (c == nc).sum()

        # skip analysis if cluster has no pixels
        if numpix <= 0:
            continue

        # make cell mask
        dset['mask'] = (c == nc)

        # do cell analysis
        cname = 'cell_%s' % str(nc + noffset).zfill(8)

        cell = single_cell_analysis(dset, **kwargs)
        cell['cell_name'] = cname
        cell['cell_id'] = nc

        for cprop in cell.keys():

            try:
                cset[cprop].append( cell[cprop] ) 
            except:
                cset[cprop] = [ cell[cprop], ]

    return cset
        


######################################################################
######################################################################


def cluster_analysis(din, varname,
                     expname = 'basic',
                     dist_to_edge = 11,
                     aux_names = [],
                     verbose = True):
    
    '''
    Performs segmentation and cell analysis 
    for a temporal data stack.

    
    INPUT
    =====
    din: input data set, including georeference.
    varname: variable name of field that is segmented

    expname: optional, name of parameter set that is used for segmentation.
    
    
    OUTPUT
    ======
    segmented_field: stacked of segmented data
    cset: set of cell properties

    '''




    #LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    # (1) get segmentation setup
    #LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    sett = predefined_setups(expname.split('-')[0])
    setup_for_later_output = copy.copy(sett)

    thresh = sett.pop('thresh')
    min_size = sett.get('min_size')



    #LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    # (2) preparation of input field
    #LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    # check if field should be inverted
    if len(expname.split('-')) == 2 and expname.split('-')[1] == 'inv':
        field = -din[varname]

    elif varname == 'bt108':
        field = -din[varname]
        thresh = -thresh 

        if verbose:
            print '... invert field'
            print
            
    else:
        field = din[varname]


    # get field dimensions
    ntime, nrow, ncol = field.shape


    segmented_field = np.zeros_like( field )
    noffset = 0

    cluster_set = {}
    dset = din.copy()


    #LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    # (3) time loop over segmentation and cell analysis
    #LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    for itime in [0, 1]: #range(ntime):
        
        if verbose:
            print '... calculate cluster properties for %d' % itime


        # get local time field
        f = field[itime]

        if thresh == 'relative50_for_mass_flux':
            thresh = special_threshold_calculations(din['lon'], din['lat'], b, 
                                                    method = 'relative50_for_mass_flux')

        # set masks
        f.data[f.mask] = thresh / 2.
        

        # get categorial field through segmnetation
        # ======================================

        if 'iclust' in din.keys():
            # field is already clustered
            c = din['iclust'][itime]
        else:
            c = seg.clustering(f, thresh, **sett)


        # remove edge connections
        c = mahotas.labeled.remove_bordering(c, dist_to_edge)
        c = seg.remove_small_clusters(c, min_size = min_size)


        segmented_field[itime] = c

        
        # cluster analysis
        # ==================

        # prepare input for cluster analysis
        # -----------------------------------

        dset['clust'] = c
        dset[varname] = f

        # set also transformed field
        trans_name = '%s_trans' % varname
        df = np.abs( f - thresh )
        dset[trans_name] = np.where( f < thresh, 0, df )

        # get distance field
        dist = scipy.ndimage.distance_transform_edt(c)
        dset['dist'] = dist
        
        dset['rel_time'] = din['rel_time'][itime]
        dset['abs_time'] = din['abs_time'][itime]
        


        var_names =  [varname, trans_name, 'dist']

        for aux_name in aux_names:            
            dset[aux_name] = din[aux_name][itime]
            var_names += [ aux_name, ]



        # perform analysis
        # ----------------
        cellset_analysis(dset, cluster_set, 
                         noffset = noffset, 
                         var_names = var_names, 
                         weight_names = ['%s_trans' % varname, 'dist'],) 

        noffset += c.max()

    for cprop in cluster_set.keys():
        cluster_set[cprop] = np.array( cluster_set[cprop] )
    # ================================================================



    return setup_for_later_output, segmented_field, cluster_set


######################################################################
######################################################################
