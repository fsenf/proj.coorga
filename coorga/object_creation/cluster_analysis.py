#!/usr/bin/env python

import numpy  as np
import scipy.spatial

import analysis_tools.grid_and_interpolation as gi

from standard_config import *


######################################################################
######################################################################


def cell_analysis(dset, 
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


def cluster_analysis(dset, cset, noffset = 0, **kwargs):


    '''
    Performs subsequent cell analysis.
    
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

        cell = cell_analysis(dset, **kwargs)
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
