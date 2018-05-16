#!/usr/bin/env python

import os, sys
import numpy as np
import skimage.transform, skimage.feature

import io_tools.netcdf as ncio
import analysis_tools.grid_and_interpolation as gi
import analysis_tools.segmentation as seg
from plotting_tools.bmaps import narval_map, Basemap

import scipy.ndimage

######################################################################
######################################################################

def register_box(f, f2, tracer_field = None):

    '''
    Shift f to match f2. An tracer field can be assigned on which basis
    the shift is determined.


    INPUT
    =====
    f: field that is shifted to match the other field
    f2: target field that is used as template to match field f
    tracer_field: optional; tuple(t, t2); field on which shift is determined


    OUTPUT
    ======
    f_trans: shifted field f to match the other field f2
    '''

    # get the tracer field, default is input field
    if tracer_field == None:
        t, t2 = f, f2
    else:
        t, t2 = tracer_field 


    # get translation shift
    shift, error, phase = skimage.feature.register_translation(f2, f)


    # setup morph transformation, shift must be reverted !!!
    shift = -shift[::-1]
    tform = skimage.transform.SimilarityTransform(translation = shift)

    
    # apply transformation (to normalized field and care for edges)
    fn = (f - f.min()) / (f.max() - f.min())

    Nan = -999.
    fn_trans = skimage.transform.warp(fn, tform, cval = Nan)
    fn_trans = np.ma.masked_equal(fn_trans, Nan)


    f_trans =  (f.max() - f.min()) * fn_trans + f.min()

    return f_trans


######################################################################
######################################################################

def cutout_field_for_subdivision(f, bsize = (100, 100)):

    '''
    Simple routine that removes edges of a field in a way 
    that the remaining subdomain field shape is dividable by 
    given box size. 
    
    The subdomain is centered in the input field. 

    
    INPUT
    =====
    f: input field
    bsize: optional, box size tuple


    OUTPUT
    ======
    fc: cutted field
    '''
    
    nrow, ncol = f.shape

    # subsampled dimensions
    nr_sub = (nrow / bsize[0]) 
    nc_sub = (ncol / bsize[1])

    # find the matching dimensions
    nr_matching = nr_sub * bsize[0]
    nc_matching = nc_sub * bsize[1]

    
    # cut edges to get a matching field
    ir_off = (nrow - nr_matching) / 2 + 1
    ic_off = (ncol - nc_matching) / 2 + 1

    ir1 = ir_off
    ir2 = nr_matching + ir_off
    ic1 = ic_off
    ic2 = nc_matching + ic_off

    fc = f[ir1:ir2, ic1:ic2]
                                
    return fc

######################################################################
######################################################################
                                

def subdivide_field(f, bsize = (100, 100)):

    '''
    A field is cutout and reshaped into given, smaler size subdomains.


    INPUT
    =====
    f: input field
    bsize: optional, box size tuple


    OUTPUT
    ======
    fdiv: subdivided field, box dimension are at the end of array
    '''
    
    nrow, ncol = f.shape

    # subsampled dimensions
    nr_sub = (nrow / bsize[0]) 
    nc_sub = (ncol / bsize[1])

    fc =  cutout_field_for_subdivision(f, bsize = bsize)
    fdiv = fc.reshape(nr_sub, bsize[0], nc_sub, bsize[1]).transpose(0,2,1,3)

    return fdiv
    
######################################################################
######################################################################

def Eulerian_change(f, f2, bsize = (100, 100)):

    '''
    Box-average Eulerian change of field.
    '''

    
    df = f2 - f

    df_sub =  subdivide_field(df, bsize = bsize)

    return df_sub.mean(axis = -1).mean(axis = -1)


######################################################################
######################################################################


def Lagrangian_change(f, f2, bsize = (100, 100), tracer_field = None):

    '''
    Lagrangian change is calculated as following:

    (i) f is interpreted as function F(x,t) and f2 as F(x, t + 1)
    (ii) f is shifted to match f2, -> F_t(x,t) =  F(x + dx, t)
    (iii) dF = F(x, t + 1) - F(x + dx, t)


    INPUT
    =====
    f: field to be shifted
    f2: target field to switch the other field is shifted
    bsize: optional, box size tuple


    OUTPUT
    ======
    df: mean Lagrangian change
    '''
 

    # (1) subdivde fields
    fsub = subdivide_field(f, bsize = bsize)
    fsub2 = subdivide_field(f2, bsize = bsize)

    if tracer_field != None:
        t, t2 = tracer_field
        tsub = subdivide_field(t, bsize = bsize)
        tsub2 = subdivide_field(t2, bsize = bsize)


    # (2) loop over outer dimensions
    nrow, ncol, nbox1, nbox2 = fsub.shape

    df = np.zeros((nrow, ncol))

    
    for irow in range(nrow):
        for icol in range(ncol):

            # take cutouts
            fc  =  fsub[irow, icol]
            fc2 = fsub2[irow, icol]

            if tracer_field != None:
                tc = tsub[irow, icol]
                tc2 = tsub2[irow, icol]
                tracer_sub = (tc, tc2)
            else:
                tracer_sub = None

            # get Lagrangian transform 
            fct = register_box(fc, fc2, tracer_field = tracer_sub)

            df[irow, icol] = (fc2 - fct).mean()

    df = np.ma.masked_invalid(df)
    df[df.mask] = 0

    return df

######################################################################
######################################################################

def semi_Lagrangian_change4tstack(f3d, **kws):


    '''
    Calculates Lagrangian change for 3d field (1st dimension time).

    It is calculated in a semi-Lagrangian way, i.e. each time step
    the difference between a fixed box and its spatially-shifted 
    counterpart is computed.

    INPUT
    =====
    f3d: 3d field (1st dimension time)
    symmetric: optional, if time trend is calculated as average
                         between forward and backward difference

    OUTPUT
    ======
    df: Lagrangian change of a field 
    '''



    return timechange(f3d, method = 'Lagrangian', **kws)


######################################################################
######################################################################

def timechange(f3d, 
               method = 'Lagrangian',
               tracer_field = None,
               symmetric = True,
               bsize = (100,100)):

    '''
    Calculates Lagrangian or Eulerian change for 3d field (1st dimension time).

    The Lagrangian change is calculated in a semi-Lagrangian way, i.e. each time step
    the difference between a fixed box and its spatially-shifted 
    counterpart is computed.

    INPUT
    =====
    f3d: 3d field (1st dimension time)
    symmetric: optional, if time trend is calculated as average
                         between forward and backward difference

    OUTPUT
    ======
    df: Lagrangian change of a field 
    '''


    # get dimensions ................................................
    ntime, nrow, ncol = f3d.shape

    df = []

    for itime in range(ntime - 1):
        f1 = f3d[itime]
        f2 = f3d[itime + 1]

        if tracer_field == None:
            tracer_sub = None
        else:
            t1 = tracer_field[itime]
            t2 = tracer_field[itime + 1]
            tracer_sub = (t1, t2)
       
        print '... calculate change for %d' % itime

        if method == 'Lagrangian':
            df.append( Lagrangian_change(f1, f2, 
                                         bsize = bsize,
                                         tracer_field = tracer_sub) ) 

        elif method == 'Eulerian':
            df.append( Eulerian_change(f1, f2, 
                                         bsize = bsize) )


    df = np.dstack( df ).transpose(2,0,1)

    if symmetric:

        # left and right values
        dfl = df[0:1]
        dfr = df[-1:]

        # center values
        dfc = 0.5 * ( df[1:] + df[:-1] )

        df = np.row_stack([dfl, dfc, dfr])

    return df

######################################################################
######################################################################

def map_prop_onto_clusterfield(prop, c):

    cs = seg.sort_clusters(c)

    # get pixel sums per object
    pmean = scipy.ndimage.measurements.mean(prop, 
                                            labels = cs,  
                                            index = range(cs.max()+1))

    pmean[0] = 0.




    return pmean[cs]


######################################################################
######################################################################

    


if __name__ == '__main__':


    try:
        itime = int(sys.argv[1])
    except:
        itime = 12
    
    date = '20160810'

    fname = '/vols/talos/home/fabian/data/icon/narval/msevi_narval_DOM01_%s.nc' % date
    
    b = ncio.read_icon_4d_data(fname, 'bt108', itime = itime)['bt108']
    b2 = ncio.read_icon_4d_data(fname, 'bt108', itime = itime + 1)['bt108']
    geo = ncio.read_icon_georef(fname)
    lon, lat = geo['clon'], geo['clat']

    bc = gi.cutout_field4box(b, (750,750), 200)
    bc2 = gi.cutout_field4box(b2, (750,750), 200)

    b = np.clip(b, 108, 230)
    b2 = np.clip(b2, 108, 230)

    bct = register_box(bc, bc2)

    bsize = (100,100)
    
    lon100 =  subdivide_field(lon, bsize ).mean(axis=-1).mean(axis=-1)
    lat100 =  subdivide_field(lat, bsize ).mean(axis=-1).mean(axis=-1)

    #lon500 =  subdivide_field(lon, (200, 200) ).mean(axis=-1).mean(axis=-1)
    #lat500 =  subdivide_field(lat, (200, 200) ).mean(axis=-1).mean(axis=-1)

    db100 = Lagrangian_change(b, b2, bsize)
    #db500 = Lagrangian_change(b, b2, (200, 200))


    # map onto cluster field
    fname = '/vols/talos/home/fabian/data/icon/narval/cluster_properties/segmented_msevi_narval_DOM01_%s_basic.nc' % date
    c = ncio.read_icon_4d_data(fname, 'iclust', itime = itime + 1)['iclust'].astype(np.int16)

    
    cc =  cutout_field_for_subdivision(c, bsize)
    clon =  cutout_field_for_subdivision(lon, bsize)
    clat =  cutout_field_for_subdivision(lat, bsize)

    db100 = db100.repeat(bsize[0], axis = 0).repeat(bsize[1], axis = 1)
    df = map_prop_onto_clusterfield(db100.data, cc)

    df = np.ma.masked_where(cc==0, df)


    # do some plotting for testing ...
    import pylab as pl

    # subsampling ....................................................
    ns = 2
    b = df[::ns, ::ns]
    
    lon = clon[::ns, ::ns]
    lat = clat[::ns, ::ns]
    # ================================================================

    
    # geo map --------------------------------------------------------
    mp = Basemap(projection='cyl',
#                 projection='mill',
                 llcrnrlat = -9,
                 urcrnrlat = 19,
                 llcrnrlon = -67.5,
                 urcrnrlon = 14.5,
                 resolution='i')
    x,y = lon, lat
    # ================================================================
 

    fig = pl.figure( figsize = (18,6))
    mp.drawcoastlines(linewidth = 1, color='gold')
    mp.warpimage()
    mp.pcolormesh(x,y,b, vmin = -5, vmax = 5, cmap = pl.cm.RdBu)
    
    pl.savefig('../pics/test_BTchange/msevi_narval_DOM01_%s_basic_%s.jpg' % (date, 
                                                                             str(itime).zfill(2)))

