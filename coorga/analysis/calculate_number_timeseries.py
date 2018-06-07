
##LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
## automatically generated by 25-Autocorrelation_of_PCF.ipynb
##TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

# TODO 
# make a routine without IO
# use atlantic masking capabilities


import tropy.analysis_tools.grid_and_interpolation as gi
import numpy as np

######################################################################
######################################################################

def trans_coords2equator(x, y):

    clon, clat = gi.xy2ll(1.*x, 1.*y)

    return gi.ll2xy(clon, clat, lon0 = 0, lat0 = 0)

######################################################################
######################################################################

def get_number_timeseries(clustfile, vname = 'smf'):

    dset = hio.read_dict_from_hdf(clustfile)
    d = dset[vname]

    xs, ys, tid, trel = d['x_mean'], d['y_mean'], d['time_id'], d['rel_time']
    clon, clat = gi.xy2ll(xs, ys)

    # determine typical land-sea-type
    lsm = (d['fraction_of_lsm_types'] * np.arange(5)).sum(axis = 1)


    # transform coordinates
    x, y = trans_coords2equator(xs,ys)

    
    # prepare cutout
    xrange = -5500., -2000.
    yrange = 0., 2000.


    rmask = (x > xrange[0]) & (x < xrange[1])& (y > yrange[0]) & (y < yrange[1])


    # loop over time
    time_ids = sorted( set(d['time_id']) )

    nvec = []
    for tid in time_ids:

        # masking 
        tmask = (d['time_id'] == tid)
    
        ncount = (rmask & tmask).sum()

        nvec.append( ncount )
    
    nvec = np.ma.masked_invalid( np.array( nvec ) )
    
    return nvec
    

######################################################################
######################################################################
