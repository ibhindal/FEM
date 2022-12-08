import numpy as np

def strainDisp2D(SF,nodeCoor,Jacob) : 
    """ 
    calculates the strain displacement matrix of a 2D element 
    Inputs : 
        SF       : Shape fun matrix (ShapeFun, ) 
        nodeCoor : node coordinate matrix (nne X 2) 
        Jacob    : Jacobian data structure (Jacobian, determinant, inverse Jacobian)  
    Returns :  
        B : Strain Displacement matrix 
    """ 
    nne = nodeCoor.shape[0] 
    sf, dsfdx, dsfde = SF[:3]
    I = Jacob['invJ'] 
    r1 = np.c_[dsfdx,np.zeros(nne)].flatten() 
    r2 = np.c_[dsfde,np.zeros(nne)].flatten() 
    R = I.dot(np.c_[r1, r2].T)
    dudx, dvdy = R[:2]
    r1 = np.c_[np.zeros(nne),dsfdx].flatten()
    r2 = np.c_[np.zeros(nne),dsfde].flatten() 
    R = I.dot(np.c_[r1, r2].T) 
    dvdx, dvdy = R[:2]
    B = np.c_[dudx, dvdy, dudy +dvdx].T 
    return B
