#!/usr/bin/env python

######################################################################
######################################################################



def predefined_setups(name):

    '''
    A set of predefined segmentation parameter sets.

    INPUT
    =====
    name: name of a parameter set

    OUTPUT
    ======
    setup: a set of segmentation parameters as dictionary
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

    return 


######################################################################
######################################################################
