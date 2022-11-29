"""
Function to calculate area
@author: abhishek
"""
import numpy as np
from numpy import *
import shapeFun
from JacobianMat import ajacob
def areaFun(qpt,qwt,xyel) :
    """
Calculates the area of an element
Inputs :
qpt : a quadrature integration points
qwt : a quadrature integration weights
Outputs :
A : area of the quadrangular element
"""
    nne =xyel.shape[0]
    nquad = qpt.shape[0]
    A = 0.0
    for ii in range(nquad) :
        for jj in range(nquad) :
            xi = qpt[ii]
            eta = qpt[jj]
            sf, dsfdx, dsfde = shapeFun.elem_shape(xi,eta,nne)
            jacobmat, detj = ajacob(dsfdx, dsfde, xyel)
            A = A + detj*qwt[ii]*qwt[jj]
    return A, xi, eta