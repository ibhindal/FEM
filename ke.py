import numpy as np
from GaussQuad import GaussQuad
from shapefunDeri import shapefunDeri
from JacobianMat import ajacob
from B_function import B_function

def kecalc(npts,E,v,xyel):
    """
    Function for calculating the element siffness matrix for an element

    Input 
        npts: number of points of the element
        E: Youngs modulus
        v: poissions ratio
        xyel: coordinates of element points (4x2 matrix)
    Output
        ke: elastic stiffness matrix
        D: stiffness matrix
     """
      
    point, weit = GaussQuad(npts)
    dmat        = E/((1-v)**2) * np.array([[1, v, 0],[v, 1, 0],[0, 0,((1-v)/2)]])

    ke = np.zeros([8,8])
    Stress = np.zeros(4)
    count = 0

    for i in range(npts) :                                      # for number of points in element
        for j in range(npts):                                   # for number of points in element
            count += 1                                           
            xi  = point[i]
            eta = point[j]
            wti = weit[i]
            wtj = weit[j]
            sn, dndx, dnde = shapefunDeri(xi, eta)              # the shape function from point i to point j
            ai, detj, I    = ajacob(dndx,dnde,xyel)             # the jacobian from point i to point j
            B  = B_function(sn, dndx, dnde ,xyel)               # the B matrix from point i to point j
            ke = ke + np.dot(np.dot(B.T, dmat), B) * detj * wti * wtj   # elastic stiffness matrix, square matrix
            
            DdotB = np.dot(dmat, B)
        # print('KeCalc - Completed!')
    return ke, DdotB
