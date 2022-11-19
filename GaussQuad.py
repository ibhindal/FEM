"""
Calculate the Gauss Quadrature points & weights
@author: abhishek
Edited by group 4
"""
# notes no error catchment for non int npts, or outside of range 2 - 6

import numpy as np
from numpy import *

def GaussQuad(npts) :
    """
    Calculates the Gauss Quadrature points & weights
    Inputs 
        ntps  : number of quadrature integration points
    Outputs 
        point : vector quadrature points (npts)
        weit  : vector quadrature weights (npts)
    """
    point = weit = zeros(npts)          #initialize vector size

    if npts == 2 :
        point[0] = 1.0/sqrt(3.0)
        point[1] = -point[0]
        weit[0] = 1.000000000000000
        weit[1] = weit[0]
        return point, weit              #individual returns means later conditional if statements are not evaluated

    if npts == 3 :
        point[0], point[1] = 0.774596669241483, 0.000000000000000
        point[2] = -point[0]
        weit[0],  weit[1] = 0.555555555555556, 0.888888888888889
        weit[2] = weit[0]
        return point, weit

    if npts == 4 :
        point[0], point[1] = 0.861136311590453,  0.339981043583856
        point[2],  point[3] = -point[1], -point[0]
        weit[0],  weit[1]  = 0.347854845137454, 0.652145154862526
        weit[2], weit[3] = weit[1], weit[0]
        return point, weit

    if npts == 5 :
        point[0] = 0.0000000000000000
        point[1], point[3] = 1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0)), 1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0))
        point[2], point[4] = -point[1], -point[3]
        weit[0], weit[1], weit[3] = 128.0/225.0, (322.0+13.0*sqrt(70.0))/900.0, (322.0-13.0*sqrt(70.0))/900.0
        weit[2], weit[4] = weit[1], weit[3]
        return point, weit

    if npts == 6 :
        point[0], point[2], point[4] = 0.2386191861, 0.6612093865, 0.9324695142
        point[1], point[3], point[5] = -point[0], -point[2], -point[4]
        weit[0], weit[2], weit[4] = 0.4679139346, 0.3607615730, 0.1713244924
        weit[1], weit[3], weit[5] = weit[0], weit[2], weit[4]
        return point, weit

    return point, weit          #error case returns vectors populated with 0's



""" print(GaussQuad(2))
print(GaussQuad(3))
print(GaussQuad(4))
print(GaussQuad(5))
print(GaussQuad(6)) """