import numpy as np

def shapefunDeri(xi, eta) : 
    """Calculate the derivative of the shape function for a 2D linear quadrangular 
    element with respect to local parametric coordinates xi and eta 
    Inputs : 
        xi, eta : local parametric coordinates 
    Returns : 
        dndx : vector of shape function derivates w.r.t xi dnde (nne x 1)
        dnde : vector of shape function derivates w.r.t eta (nne x 1)
        sn : Shape function (nne x 1)
     """ 
    sn = np.zeros(4)
    dndx, dnde = np.zeros(4), np.zeros(4)
    sn = np.array([(1 - xi) * (1 - eta) / 4, (1 + xi) * (1 - eta) / 4, (1 + xi) * (1 + eta) / 4, (1 - xi) * (1 + eta) / 4])
    dndx = np.array([-(1 - eta) / 4, (1 - eta) / 4, (1 + eta) / 4, -(1 + eta) / 4])
    dnde = np.array([-(1 - xi) / 4, -(1 + xi) / 4, (1 + xi) / 4, (1 - xi) / 4])
    
    return sn, dndx, dnde




