"""
Calculate the Gauss Quadrature points & weights
@author: abhishek
Edited by group 4
"""
# notes no error catchment for non int npts, or non quadrangular

import numpy as np

def GaussQuad(npts=4) :
    
    """
    Calculates the Gauss Quadrature points & weights
    Inputs 
        ntps  : number of quadrature integration points
    Outputs 
        point : vector quadrature points (npts)
        weit  : vector quadrature weights (npts)
    """
    point, weit = np.zeros(4),np.zeros(4) # Initialize vector size
    ##### Point and Weit values, given #####
    point[0], point[1] = 0.861136311590453,  0.339981043583856
    point[2], point[3] = -point[1], -point[0]
    weit[0],  weit[1]  = 0.347854845137454, 0.652145154862526
    weit[2],  weit[3]  = weit[1], weit[0]
    return point, weit
