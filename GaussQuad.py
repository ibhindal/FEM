"""
Calculate the Gauss Quadrature points & weights
@author: abhishek
"""
import numpy as np
from numpy import *
def GaussQuad(npts) :

# Calculates the Gauss Quadrature points & weights

# Inputs :
# ntps : number of quadrature integration points
# Outputs :
# point : vector quadrature points (npts)
# weit : vector quadrature weights (npts)

    if npts == 2 :
        point, weit = zeros(2), zeros(2)
        point[0] = 1.0/sqrt(3.0)
        point[1] = -point[0]
        weit[0] = 1.000000000000000
        weit[1] = weit[0]
    if npts == 3 :
        point, weit = zeros(3), zeros(3)
        point[0] = 0.774596669241483
        point[1] = 0.000000000000000
        point[2] = -point[0]
        weit[0] = 0.555555555555556
        weit[1] = 0.888888888888889
        weit[2] = weit[0]
    if npts == 4 :
        point, weit = zeros(4), zeros(4)
        point[0] = 0.861136311590453
        point[1] = 0.339981043583856
        point[2] = -point[1]
        point[3] = -point[0]
        weit[0] = 0.347854845137454
        weit[1] = 0.652145154862526
        weit[2] = weit[1]
        weit[3] = weit[0]
    if npts == 5 :
        point, weit = zeros(5), zeros(5)

        point[0] = 0.0000000000000000
        point[1] = 1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0))
        point[2] = -point[1]
        point[3] = 1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0))
        point[4] = -point[3]
        weit[0] = 128.0/225.0
        weit[1] = (322.0+13.0*sqrt(70.0))/900.0
        weit[2] = weit[1]
        weit[3] = (322.0-13.0*sqrt(70.0))/900.0
        weit[4] = weit[3]
    if npts == 6 :
        point, weit = zeros(6), zeros(6)
        point[0] = 0.2386191861
        point[1] = -point[0]
        point[2] = 0.6612093865
        point[3] = -point[2]
        point[4] = 0.9324695142
        point[5] = -point[4]
        weit[0] = 0.4679139346
        weit[1] = weit[0]
        weit[2] = 0.3607615730
        weit[3] = weit[2]
        weit[4] = 0.1713244924
        weit[5] = weit[4]
    return point, weit

