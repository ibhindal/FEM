import numpy as np
import scipy as sp
from numpy import *
from scipy import *

"""
Python subroutine to calculate shape functions and their derivatives for a
quadrangular element, for both linear (4 noded) and quadradric (8 noded) element
"""
def elem_shape(xi,eta,nne) :
    """
    calculates the shape function at specified points in the local coordinate specified by xi and eta
    Inputs :
        xi, eta : local coordinates in the quadrangular element
        nne : number of nodes per element
    Returns :
        sn : vector of shape functions (1 by nne)
        dndx : vector derivatives of shape functions w.r.t xi (1 by nne)
        dnde : vector derivatives of shape functions w.r.t eta (1 by nne)
    """
    sn,  dndx, dnde = zeros(nne)

# case for 4 nodes per element (linear)
    if nne == 4 :
        sn[0]   = (1.0-xi)*(1.0-eta)/4.0
        sn[1]   = (1.0+xi)*(1.0-eta)/4.0
        sn[2]   = (1.0+xi)*(1.0+eta)/4.0
        sn[3]   = (1.0-xi)*(1.0+eta)/4.0

        dndx[0] = -(1.0-eta)/4.0
        dnde[0] = -(1.0-xi )/4.0
        dndx[1] =  (1.0-eta)/4.0
        dnde[1] = -(1.0+xi )/4.0
        dndx[2] =  (1.0+eta)/4.0
        dnde[2] =  (1.0+xi )/4.0
        dndx[3] = -(1.0+eta)/4.0
        dnde[3] =  (1.0-xi )/4.0
        return sn, dndx, dnde

# case for 8 nodes per element (quadratic) cyclic numbering
    if nne == 8 :
        sn[0]   = -(1.0-xi)*   (1.0-eta)*(1.0+xi+eta)/4.0
        sn[1]   =  (1.0-xi**2)*(1.0-eta)/2.0
        sn[2]   =  (1.0+xi)*   (1.0-eta)*(-1.0+xi-eta)/4.0
        sn[3]   =  (1.0+xi)*   (1.0-eta**2)/2.0
        sn[4]   =  (1.0+xi)*   (1.0+eta)*(-1.0+xi+eta)/4.0
        sn[5]   =  (1.0-xi**2)*(1.0+eta)/2.0
        sn[6]   =  (1.0-xi)*   (1.0+eta)*(-1.0-xi+eta)/4.0
        sn[7]   =  (1.0-xi)*   (1.0-eta**2)/2.0

        dndx[0] = -((eta-1.0)*(2.0*xi+eta))/4.0
        dnde[0] = -((xi-1.0)*(xi+2.0*eta))/4.0
        dndx[1] =  (eta-1.0)*xi
        dnde[1] =  ((xi-1.0)*(xi+1.0))/2.0
        dndx[2] = -((eta-1.0)*(2.0*xi-eta))/4.0
        dnde[2] = -((xi+1.0)*(xi-2.0*eta))/4.0
        dndx[3] = -((eta-1.0)*(eta+1.0))/2.0
        dnde[3] = -eta*(xi+1.0)
        dndx[4] =  ((eta+1.0)*(2.0*xi+eta))/4.0
        dnde[4] =  ((xi+1.0)*(xi+2.0*eta))/4.0
        dndx[5] = -(eta+1.0)*xi
        dnde[5] = -((xi-1.0)*(xi+1.0))/2.0
        dndx[6] =  ((eta+1.0)*(2.0*xi-eta))/4.0
        dnde[6] =  ((xi-1.0)*(xi-2.0*eta))/4.0
        dndx[7] =  ((eta-1.0)*(eta+1.0))/2.0
        dnde[7] =  eta*(xi-1.0)
    return sn, dndx, dnde