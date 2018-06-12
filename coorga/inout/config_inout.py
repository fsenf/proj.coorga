#!/usr/bin/env python

import os
import json

from coorga.object_creation.segmentation_config import predefined_collections

######################################################################
######################################################################

def write_cluster_collection_config(*args, **kwargs):

    '''
    Writes a predefined cluster collection config to disk.

    
    Parameters
    -----------
    args :  list of arguments
        arguments forwarded to predefined_collections function

    fname : str, optional, default = None
        name of output config file
        if None: filename is automatically generated from arguments
    '''



    # make fname -----------------------------------------------------
    fname = kwargs.get('fname', None)

    if fname is None:
        arg_str = '_'.join( args )
        fname = 'collection_config_%s.json' % arg_str


    # get collection -------------------------------------------------
    coll = predefined_collections(*args)


    # do output ------------------------------------------------------
    print 'output collection config into ', fname
    dict2json( fname, coll)    

    return

######################################################################
######################################################################

def json2dict( fname ):

    '''
    Reads a Json File.


    Parameters
    ----------
    fname : str
        input filename


    Returns
    --------
    d : dict
        content of json as dict

    '''

    with open(fname, "r") as read_file:
        d = json.load(read_file)

    return d

######################################################################
######################################################################


def read_cluster_collection_config( fname ):

    '''
    Reads a Config for making a Clustering Collection.


    Parameters
    ----------
    fname : str
        input filename


    Returns
    --------
    coll : dict
        setup for a cluster collection

    '''

    return json2dict

######################################################################
######################################################################


def dict2json( fname, d ):

    '''
    Reads a Config for making a Clustering Collection.


    Parameters
    ----------
    fname : str
        input filename
    d : dict
        data for output

    Returns
    --------
    None
    '''

    with open(fname, "w") as write_file:
        json.dump(d, write_file, indent = 4)


    return

######################################################################
######################################################################




if __name__ == '__main__':

    # testing
    write_cluster_collection_config('narval', 'imf', '201608')
    
    print read_cluster_collection_config('collection_config_narval_imf_201608.json')    

