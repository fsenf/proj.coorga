#!/usr/bin/env python


from  tropy.standard_config import *


######################################################################
######################################################################



def predefined_setups(name):

    '''
    A set of predefined segmentation parameter sets.

    Parameters
    ----------
    name : str
        name of a parameter set

    Returns
    ----------
    setup : dict
       predefined set of segmentation parameters as dictionary
    '''

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

    elif name == 'exp010':
        b = predefined_setups('basic')
        b.update( dict( thresh = 50 ) )
        b.update( dict( min_size = 36 ) )
        return b

    elif name == 'exp011':
        b = predefined_setups('basic')
        b.update( dict( thresh = 'relative50_for_mass_flux' ) )
        b.update( dict( min_size = 36 ) )
        b.update(dict( name = name ))
        return b

    elif name == 'exp012':
        b = predefined_setups('basic')
        b.update( dict( thresh = 'full_domain_relative50_for_mass_flux' ) )
        b.update( dict( min_size = 36 ) )
        b.update(dict( name = name ))
        return b


    elif name == 'exp022':
        b = predefined_setups('exp011')
        b.update( dict( min_size = 4 ) )
        b.update(dict( name = name ))
        return b

    elif name == 'exp023':
        b = predefined_setups('exp022')
        b.update( dict(  cluster_method = 'connect' ) )
        b.update(dict( name = name ))
        return b


    elif name == 'exp123':
        b = predefined_setups('basic')
        b.update( dict( thresh = 0.5 ) )
        b.update( dict( min_size = 16 ) )
        return b


    # setups made after review of PCF paper 2019
    elif name == 'exp030':
        
        s = dict( min_size = 36, 
                  base_field = 'max_w',
                  thresh = 0.7,
                  name = name )
        
        b = predefined_setups('basic')
        b.update( s )
        return b

    elif name == 'exp031':
        
        s = dict( min_size = 36, 
                  base_field = 'max_mflux',
                  thresh = 0.5,
                  name = name )
        
        b = predefined_setups('basic')
        b.update( s )
        return b

    elif name == 'exp032':
        
        s = dict( min_size = 36, 
                  base_field = 'max_cflux',
                  thresh = 6e-4,
                  name = name )
        
        b = predefined_setups('basic')
        b.update( s )
        return b

    elif name == 'exp033':
        
        s = dict( min_size = 36, 
                  base_field = 'imf',
                  thresh = 800.,
                  name = name )
        
        b = predefined_setups('basic')
        b.update( s )
        return b

    elif name == 'exp034':
        
        s = dict( min_size = 36, 
                  base_field = 'tcwp',
                  thresh = 3.,
                  name = name )
        
        b = predefined_setups('basic')
        b.update( s )
        return b


    return 


######################################################################
######################################################################


def available_expnames():

    '''
    Print available Experiment names.
    '''

    elist = ['basic', 'exp001', 'exp002', 'exp003', 'exp004', 'exp005', 
             'basicde', 'exp010', 'exp011', 'exp022', 'exp023', 'exp123',
             'epx030', 'exp031', 'exp032', 'exp033', 'exp034']

    return elist

######################################################################
######################################################################


def predefined_collections(regtype, varname, date):
    
    '''
    A set of predefined collections for generating temporal stacks 
    of clusterfiles.

    Parameters
    ----------
    regtype : string
        region name

    varname : string
        variable name

    date : string, typically in format %Y%m
        string that specifies the date or part of the date


    Returns
    --------
    ftype_combinations : list
        list of datatype - modes - combined together

    kws : dict
        keywords for filelist generation, used in the read_cluster_props function
    
    addlist : list
        names of specific variable to be added
    '''


    # select file parameters -----------------------------------------
    if regtype == 'narval' and varname == 'bt108':
        
        ftype_combinations = [  ('msevi', 'obs'),  
                                ('synsat', 'sim'),
                                ('trans', 'trans')]

        narval_dir = '%s/icon/narval/synsat/' % local_data_path

        kws = dict(
            filepart = '_narval_DOM01_',
            fdir = '%s/cluster_properties' % narval_dir,
            date = date)

        addlist = ['fraction_of_lsm_types', ]


    elif regtype == 'narval' and varname in ['smf', 'rr', 'imf']:
        
        ftype_combinations = [  ('icon-narval', varname) ] 


        kws = dict(
            filepart = '_dom01_%s_' % varname,
            fdir = '%s/icon/narval/variables/cluster_properties'  % local_data_path,
            # fdir = '%s/icon/narval/var4workflowtest/cluster_properties'  % local_data_path,
            date = date)

        addlist = ['fraction_of_lsm_types', ]


    elif regtype == 'hdcp2':

        ftype_combinations = [('hdfd', 'obs'),]

        

        kws = dict(
            filepart = '_trop_seviri',
            time_masking = False,
            fdir = '%s/hdcp2/cluster_properties' % local_data_path,
            date = date,
            )

        addlist = []


    elif regtype == 'icon-lem':

        ftype_combinations = [('msevi', 'msevi'),
                              ('synsat1', 'synsat1'),
                              ('synsat2', 'synsat2'),
                              ('synsat3', 'synsat3')]

        

        kws = dict(
            filepart = '_3d_coarse_icon-lem-tstack_DOM',
            time_masking = False,
            fdir = '%s/icon/lem-de/cluster_properties' % local_data_path,
            date = date,
            )

        addlist = []



    return ftype_combinations, kws, addlist

######################################################################
######################################################################
