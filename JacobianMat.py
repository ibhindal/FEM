
# Python subroutine to evaluate Jacobian for a 2D quadrangular element 5

import numpy as np

def ajacob(dndx,dnde,nodecoor) : 
    
    """ 
    calculates Jacobian at the local coordinate point (xi,eta) 
    Inputs : 
        dndx : vector derivatives of shape functions w.r.t xi (nne by 1) 
        dnde : vector derivatives of shape functions w.r.t eta (nne by 1) 
        nodecoor : node coordinate matrix (nne X 2) 
    
    Returns : 
        ai   : the Jacobian matrix (2 by 2) 
        detj : the determinant of the Jacobian matrix
        I    : the inverse of the jacobian
    """ 
    ai = np.zeros([2,2])                 # Initialize 2x2 Jacobian matrix
    ai[0,0] = np.dot(dndx,nodecoor[:,0]) # Sum of dndx with xi 
    ai[0,1] = np.dot(dndx,nodecoor[:,1]) # Sum of dndx with yi 
    ai[1,0] = np.dot(dnde,nodecoor[:,0]) # Sum of dnde with xi 
    ai[1,1] = np.dot(dnde,nodecoor[:,1]) # Sum of dnde with yi 
    
    detj = ai[0,0] * ai[1,1] - ai[1,0] * ai[0,1] # Determinant of the Jacobian 
     
    I = np.zeros # Initialize inverse of Jacobian
    
    if detj == 0 : 
        print("Singular matrix !!!") 
    else : 
        I = np.linalg.inv(ai) # Calculate inverse of the Jacobian
        
    return ai, detj, I