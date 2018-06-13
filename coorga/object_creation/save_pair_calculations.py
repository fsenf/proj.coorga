#!/usr/bin/env python

import sys, os, glob
import datetime
import numpy as np

import tropy.io_tools.hdf as hio
import tropy.analysis_tools.grid_and_interpolation as gi

from coorga.analysis.calculate_pcf4clust import calculate_pcf4clust
from coorga.analysis.calculate_pcount4clust import calculate_pcount4clust


######################################################################
######################################################################

def main(collection_file, number_density_file,
                         pcount_variables = None,
                         weight = 1., 
                         rMax = 500., 
                         dr = 20.,
                         no_pcf = False):
    
    '''
    Saves pair correlations.
    

    Parameters
    ----------
    collection_file : str
        filename of cluster collection

    number_density_file : str
        filename of number density file

    pcount_variables : dict, optional, default = None
        set of pcount variable names with bins

    weight : float, optional, default = 1.
        weight between equi-distant and equal-area range binning
        0 = equal-area
        1 = equi-distant

    rMax : float, optional, default = 500.
        largest range bin

    dr : float, optional, default = 20.
        range ring interval
 

    Returns
    --------
    None
    '''
    # ================================================================
        

    # input collection data ------------------------------------------
    cset = hio.read_dict_from_hdf( collection_file ) 
    # ================================================================
 
    if not no_pcf:
        # input number denisty data ---------------------------------------
        ndset = hio.read_dict_from_hdf( number_density_file )
    
        egrid = (ndset['xgrid'], ndset['ygrid'])
        nd_ref = ndset['number_density']
        # ================================================================
        
    
        # pcf calculations ----------------------------------------------- 
        calculate_pcf4clust( cset, egrid, nd_ref, 
                             weight = weight, 
                             rMax = rMax, 
                             dr = dr, 
                             verbose = True)
        # ================================================================

   
    # pcount calculation ---------------------------------------------
    if pcount_variables is not None:
        varnames = pcount_variables.keys()
        nvar = len( varnames )

        for vname in varnames:
            vbins = pcount_variables[vname]
            
            calculate_pcount4clust( cset, vname, vbins, 
                            pcount_name = 'pcount_%s' % vname, 
                            rMax = rMax, 
                            dr = dr,
                            verbose = True)
    # ================================================================


    # update cluster data --------------------------------------------
    print '... update collection file %s' % collection_file
    hio.update_dict_in_hdf(collection_file, cset)
    # ================================================================

    return 



######################################################################
######################################################################
