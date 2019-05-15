
'''
This module collects methods for threshold calculation.
'''

import numpy as np

import tropy.analysis_tools.statistics as stats

from coorga.inout.atlantic_masking import atlantic_masking

####################################################################
####################################################################

def special_threshold_calculations(lon, lat, f, method = None):
    
    '''
    Calculates a special threshold for a field.


    Parameters
    ----------
    lon : numpy array, 2-dim (nrow, ncol)
        longitude field

    lat : numpy array, 2-dim (nrow, ncol)
        latitude field

    f : numpy array, 2-dim (nrow, ncol)
        addtional field on which threshold calculation is based 
    
    method : optional, str
        string that identifies a use method for threshold calculation
    
        method = 'relative50_for_mass_flux' takes threshold where 50% of the domain (Atlantic)
        integrated mass flux is enclosed


    Returns
    -------
    thresh : numpy float value
         calculated threshold value


    Notes
    -----
    Special threshold calculations are needed in some segmentation applications.
    
    '''
    
    if method ==  'relative50_for_mass_flux':
        
        # peform region masking
        fm = atlantic_masking(lon, lat, f, mask_value = 0.)

        # calculate cummulative data sum
        fp = np.where(fm > 0., fm, 0.)
        csf = stats.cumsum_data_fraction(fp)
        
        thresh = fp[csf < 50.].min()
 
    elif 'full_domain_relative' in method:

        method_format = 'full_domain_relative%d_for_mass_flux'
        
        method_list = method.split('_')

        # extract number
        last_two_digits = float( method_list[2][-2:] )

        # assume that number is used as threshold
        thresh = last_two_digits

        # no region masking
        #fm = atlantic_masking(lon, lat, f, mask_value = 0.)
        fm = f

        # calculate cummulative data sum
        fp = np.where(fm > 0., fm, 0.)
        csf = stats.cumsum_data_fraction(fp)
        
        thresh = fp[csf < thresh].min()
        
    return thresh
