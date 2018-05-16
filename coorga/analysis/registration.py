#!/usr/bin/env python

# load libraries -----------------------------------------------------
import sys, os, glob, copy
import matplotlib
matplotlib.use('AGG')

import numpy as np
import datetime
import scipy.ndimage
import skimage.feature 

import io_tools.hdf as hio
import io_tools.netcdf as ncio

import io_tools.radolan as rado
import analysis_tools.grid_and_interpolation as gi
import analysis_tools.parallax
import analysis_tools.segmentation as seg
from analysis_tools.make_hrv_upscaling import make_hrv_upscaling 

from standard_config import *

    
######################################################################
######################################################################

def Lagrangian_change(f1, f2,
                      bsize = (100, 100),
                      tracer_field = None,
                      return_index = False):

    
    '''
    Lagrangian change is calculated as following:

    (i) f1 is interpreted as function F(x,t) and f2 as F(x, t + 1)
    (ii) f1 is shifted to match f2, -> F_t(x,t) =  F(x + dx, t)
    (iii) dF = F(x, t + 1) - F(x + dx, t)


    INPUT
    =====
    f1: field to be shifted
    f2: target field to switch the other field is shifted
    bsize: optional, box size tuple


    OUTPUT
    ======
    df: mean Lagrangian change
    '''

    if bsize[0] != bsize[1]:
        raise Exception('Registration is not able not handle different bsizes')

    if np.mod(bsize[0], 2) == 0:
        raise Exception('Registration is not able not handle even bsize')
        
    # get field shape ................................................
    nrow, ncol = f1.shape



    # determine subsampling parameters ...............................
    nrow_sub = nrow / bsize[0]
    ncol_sub = ncol / bsize[1]

    nrow_fit = nrow_sub * bsize[0]
    ncol_fit = ncol_sub * bsize[1]

    ir_offset = (nrow - nrow_fit) / 2
    ic_offset = (ncol - ncol_fit) / 2

    ir0 = ir_offset + bsize[0] / 2
    ic0 = ic_offset + bsize[1] / 2
    

    row_range = bsize[0] * np.arange(nrow_sub) + ir0
    col_range = bsize[1] * np.arange(ncol_sub) + ic0
    

    df = np.zeros((nrow_sub, ncol_sub))
    ir_shift =  np.zeros((nrow_sub, ncol_sub))
    ic_shift =  np.zeros((nrow_sub, ncol_sub))

    # get 3d field with time at first place
    f = np.dstack([f1, f2]).transpose(2,0,1)
    
    for i, irow in enumerate(row_range):
        print '... calculate df for row (%d / %d )' % (irow, nrow)
        for j, icol in enumerate(col_range):

            ind0 = (irow, icol)
            tubesize = bsize[0]
            
            ind, freg = cutout_registered_box(f, ind0, tubesize,
                                                  Niter = 3,
                                              #   Niter = 20,
                                                  itime = 0)

            df[i, j] = (freg[1] - freg[0]).mean()
            ir_shift[i, j] = ind[0,-1] - irow
            ic_shift[i, j] = ind[1,-1] - icol

    if return_index:
        return ir_shift, ic_shift, df
    else:
        return df


######################################################################
######################################################################

def semi_Lagrangian_change4tstack(f3d,
                                    tracer_field = None,
                                    symmetric = True,
                                    bsize = (100,100)):

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
        df.append( Lagrangian_change(f1, f2, 
                                        bsize = bsize,
                                        tracer_field = tracer_sub) ) 

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


def cutout_registered_box(f, ind0, bsize,
                            cmethod = 'register',
                            nstep = 1,
                            Niter = 3,
                            itime = None,
                            fmin = -1e23, 
                            regist_method = 'median'):
    
    '''
    This routine makes a Lagranian cutout of a smaller 3d tube from a 3d field. 
    It uses image registration method to find shift with maximum cross-correlation.


    INPUT
    =====
    f: 3d field (time axis at 1st dimension, time array centered around considered time)
    ind0: center index (conserved for central time index)
    bsize: box size

    Niter: optional, number of iterations for successive registrations
    fmin: optional, minimum threshold for data (assumes mountains not valleys in data)
    regist_method: optional, method to choose for registration, 
                   e.g., average vs. individual index shift


    OUTPUT
    ======
    ind: index array used for registration, function of time
    ftube: 3d tube of registered field cutout
    
    '''


    
    # field properties
    ntime = f.shape[0]

    if itime == None:
        itime = ntime / 2
    

    # set inital index as time array  ................................
    ind = np.array(ind0).repeat(ntime).reshape(2,-1)

    # initial cutout
    ftube = gi.tube_cutout4box(f, ind, bsize)


    eps = 1e23
    # iterate over shift .............................................
    for ni in range(Niter):
        
        
        # mask
        freg =  np.where((ftube.data > fmin) & (ftube.mask == False) , ftube.data, fmin)
        
        # determine shift from image registration
        ishift = get_tube_shift_from_registration(freg,
                                                    cmethod = cmethod,
                                                    nstep = nstep,
                                                    cumulative = False, )
        #                                                  edge_filter = 'Hamming')
        


        # choose the shift to apply ..................................
        if regist_method == 'median':
            # calculate median shift
            med_shift = np.median(ishift, axis = 0) * np.ones_like(ishift)
            applied_shift = med_shift

        elif regist_method == 'centered_median':

            # calculate median shift just from the center part
            cmed_shift = np.median(
                ishift[ntime/4.:3*ntime/4.], 
                                          axis = 0) * np.ones_like(ishift)
            applied_shift = cmed_shift
                
        elif regist_method == 'smoothed_index':
            applied_shift = scipy.ndimage.median_filter(ishift, 3)

        else:
            print '...take unfiltered index shift', ishift
            applied_shift = ishift


        # sump it up and center it around central time
        cum_shift = np.round( applied_shift.cumsum(axis=0) ).transpose()


        cum_shift[0] -= cum_shift[0, itime]
        cum_shift[1] -= cum_shift[1, itime]

#        ind -= cum_shift

        ind_old = copy.copy(ind)
        ind = ind - cum_shift

        eps = np.abs(ind - ind_old).sum()
        # print eps, ind
#        if eps < 0.1:
#            break
        
        # and cutout
        ftube = gi.tube_cutout4box(f, ind, bsize)



    return ind, ftube

######################################################################
######################################################################


def get_tube_shift_from_registration(f3d, 
                                        cmethod = 'register',
                                        cumulative = True, 
                                        edge_filter = None,
                                        nstep = 1):

    '''
    Performs fourier shift calculation via image registration.

    INPUT
    =====
    f3d: 3d field, 1st dimension is time


    OUTPUT
    ======
    ishift: cumulated shift index
    '''


    ntime, nrows, ncols = f3d.shape

    if edge_filter == 'Hamming':
        g3d = whamm(f3d)
    else:
        g3d = f3d

    s = []
    for n in range(ntime - nstep):

        f1 = g3d[n]
        f2 = g3d[n + nstep]

        if cmethod == 'register':
            shift, error, diffphase = skimage.feature.register_translation(f1, f2, 10)
        
        elif cmethod == 'crosscorr':
            nshift = 5
            shift = shift_max_crosscorr(f1, f2, nshift = nshift)
            shift = np.array(shift)


        s.append(shift / nstep)

    # append edge
    s.insert(0, [0,0])

    # append edge
    for i in range(nstep -1):
        s.append( s[-1] )


    s = np.array(s)
        

    
    if cumulative:
        return np.round( s.cumsum(axis=0) ).transpose()
    else:
        return s


######################################################################
######################################################################


def whamm(f3d):
    '''
    Simple Hamming window implementation for testing purposes.
    '''
    
    ntime, nrow, ncol = f3d.shape

    f3d_hamm = f3d * hamming_filter(ncol)

    f3d_trans = f3d_hamm.transpose(0,2,1) * hamming_filter(nrow)

    return f3d_trans.transpose(0,2,1)

######################################################################
######################################################################

def hamming_filter(n, no_edge = True):

    if no_edge:
        return np.hamming(n)
    else:
     
        nhalf = n / 2
     
        f = np.ones(nhalf)    
        nedge = n / 3
     
        f[:nedge] = np.hamming(2 * nedge)[:nedge]
        
     
        
        if np.mod(n, 2) == 1:
            return np.hstack([f,1, f[::-1]])
     
        else:
            return np.hstack([f, f[::-1]])



        
######################################################################
######################################################################

def cross_correlation_matrix(f1, f2, nshift = 3):

    '''
    Calculates cross-correlation matrix between two fields,
    target field f1 and search field f2 which is sequentially 
    shifted.

    INPUT
    =====
    f1: target field
    f2: search field

    OUTPUT
    ======
    cc: cross-correlation matrix

    '''



    # get dimension of fields
    nrow, ncol = f1.shape

    if nrow != ncol:
        raise Exception('Assume equal size of dimension!')


    # initialized cc marix
    cc = np.zeros((2*nshift + 1, 2*nshift + 1))


    # set center index and boxsize
    icenter = (nrow / 2, ncol / 2)
    bsize = nrow - nshift

    # cutout base /target field
    c1 = gi.cutout_field4box(f1, icenter, bsize)

    # loop over shifts
    for ir_shift in range(-nshift, nshift + 1):
        for ic_shift in range(-nshift, nshift + 1):
 
            # indices of cc matrix
            i = ir_shift + nshift
            j = ic_shift + nshift

            # shift center index
            icenter = (nrow / 2 + ir_shift, ncol / 2 + ic_shift)

            # cutout search field
            c2 = gi.cutout_field4box(f2, icenter, bsize)


            # calculate cross correlation
            cc[i, j] = np.corrcoef(c1.flatten(), c2.flatten())[0,1]

    return cc

######################################################################
######################################################################


def shift_max_crosscorr(f1, f2, nshift = 3):

    '''
    Calculates index shift for maximum cross-correlation.


    INPUT
    =====
    f1: target field
    f2: search field

    OUTPUT
    ======
    ind: index shift

    '''

    # calculate cross correlation matrix
    cc = cross_correlation_matrix(f1, f2, nshift = nshift)

    #  get maximum index relative to (0,0)
    imax = gi.i2iset(cc, cc.argmax())

    # get dimension of fields
    nrow, ncol = cc.shape

    # set center index 
    icenter = (nrow / 2, ncol / 2)
 
    # maximum index relative to center index
    return icenter[0] - imax[0], icenter[1] - imax[1]


######################################################################
######################################################################



if __name__ == '__main__':

    f = np.zeros((295, 295))
    f = np.zeros((302, 302))
    print Lagrangian_change(f,f, bsize = (51,51))
