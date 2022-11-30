
# Python subroutine to evaluate Jacobian for a 2D quadrangular element
# 5
# @author: abhishek
import numpy as np

def ajacob(dndx,dnde,nodecoor) : 
    """ 
    calculates Jacobian at the local coordinate point (xi,eta) 
    Inputs : 
    dndx : vector derivatives of shape functions w.r.t xi (1 by nne) 
    dnde : vector derivatives of shape functions w.r.t eta (1 by nne) 
    nodecoor : node coordinate matrix (nne X 2) 
    
    Returns : 
    ai : the Jacobian matrix (2 by 2) 
    detj : the determinant of the Jacobian matrix
    """ 
    ai = np.zeros([2,2]) 
    ai[0,0] = np.dot(dndx,nodecoor[:,0]) #sum of dndx with xi 
    ai[0,1] = np.dot(dndx,nodecoor[:,1]) #sum of dndx with yi 
    ai[1,0] = np.dot(dnde,nodecoor[:,0]) #sum of dnde with xi 
    ai[1,1] = np.dot(dnde,nodecoor[:,1]) #sum of dnde with yi 
    
    # determinant of the Jacobian 
    detj = ai[0,0]*ai[1,1]- ai[1,0]*ai[0,1] 
    
    # calculate the inverse of the jacobian 
    I= np.zeros
    
    if detj == 0 : 
        print("Singular matrix !!!") 
    else : 
        #print "Jacobian for element ", elemc, "is : ", ai 
        # #print "determinant of Jacobian is :", detj 
         
        I = np.linalg.inv(ai) 
        
        
    return ai, detj, I