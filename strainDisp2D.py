import numpy as np

def strainDisp2D(SF,nodeCoor,Jacob) : 
    
    """ 
    Calculates the strain displacement matrix of a 2D element 
    Inputs : 
        SF       : Shape function matrix 
        nodeCoor : Node coordinate matrix (nne X 2) 
        Jacob    : Jacobian data structure (Jacobian, determinant, inverse Jacobian)  
    Returns :  
        B        : Strain Displacement matrix 
    """ 
    
    nne = nodeCoor.shape[0]                                 # Number of elements
    sf, dsfdx, dsfde = SF['sf'], SF['dndx'], SF['dnde']     # Assigning variables to corresponding SF dictionary items
    I = Jacob['invJ']                                       # Assigning inverse Jacobian variable to corresponding Jacob dictionary item 
    r1 = np.c_[dsfdx,np.zeros(nne)].flatten()               # ?
    r2 = np.c_[dsfde,np.zeros(nne)].flatten()               # ?
    R = I.dot(np.c_[r1, r2].T)                              # ?
    dudx, dudy = R[:2]                                      # ? 
    r1 = np.c_[np.zeros(nne),dsfdx].flatten()               # ?
    r2 = np.c_[np.zeros(nne),dsfde].flatten()               # ?
    R = I.dot(np.c_[r1, r2].T)                              # ?
    dvdx, dvdy = R[:2]                                      # ?
    B = np.c_[dudx, dvdy, dudy +dvdx].T                     # ?
    return B
