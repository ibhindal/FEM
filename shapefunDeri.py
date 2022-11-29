import numpy as np
import scipy



def shapefunDeri(xi, eta) : 
    """Calculate the derivative of the shape function for a 2D linear quadrangular 
    element with respect to local parametric coordinates xi and eta 
    Inputs : 
    xi, eta : local parametric coordinates 
    Returns : 
    dndx : vector of shape function derivates w.r.t xi dnde 
     : vector of shape function derivates w.r.t eta
     """ 
    sn = np.zeros(4)
    sn[0] = (1.0-xi)*(1.0-eta)/4.0
    sn[1] = (1.0+xi)*(1.0-eta)/4.0
    sn[2] = (1.0+xi)*(1.0+eta)/4.0
    sn[3] = (1.0-xi)*(1.0+eta)/4.0
    dndx, dnde = np.zeros(4), np.zeros(4)
    dndx[0] = -(1.0-eta)/4.0
    dnde[0] = -(1.0-xi)/4.0
    dndx[1] = (1.0-eta)/4.0
    dnde[1] = -(1.0+xi)/4.0
    dndx[2] = (1.0+eta)/4.0
    dnde[2] = (1.0+xi)/4.0
    dndx[3] = -(1.0+eta)/4.0
    dnde[3] = (1.0-xi)/4.0
    
    return sn, dndx, dnde