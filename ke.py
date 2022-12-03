import numpy as np
from GaussQuad import GaussQuad
from shapefunDeri import shapefunDeri
from JacobianMat import ajacob
from B_function import B_function

def kecalc(npts,E,v,xyel):
    """
    Function for calculating the element siffness matrix
    
    Input 
        npts: number of points of the element
        E: Youngs modulus
        v: poissions ratio
        xyel: coordinates of element points
    Output
        ke:  
        D: tensor stiffness matrix
     """
      
    point, weit = GaussQuad(npts)
    dmat        = E/((1-v)**2) * np.array([[1, v, 0],[v, 1, 0],[0, 0,((1-v)/2)]])
    D           = dmat

    #D = {}
    #for i in range(npts):          #for each material make a D matrix calculation
    #    D[i] = dmat

    ke = np.zeros([8,8])
    count = 0

    for ii in range(npts) :
        for jj in range(npts):
            count += 1
            #print(count)
            xi  = point[ii]
            eta = point[jj]
            wti = weit[ii]
            wtj = weit[jj]
            sn, dndx, dnde = shapefunDeri(xi, eta)
            ai, detj, I    = ajacob(dndx,dnde,xyel)
            B  = B_function(sn, dndx, dnde ,xyel)
            ke = ke + B.T.dot(D[ii]).dot(B) * detj * wti * wtj
        print('KeCalc - Completed!')
    return ke, D