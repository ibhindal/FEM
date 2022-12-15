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
        v: Poissions ratio
        xyel: Coordinates of element points (4x2 matrix)
    Output
        ke: Elastic stiffness matrix
        DdotB: Strain Displacement Matrix dot product with Elasticity Tensor
     """
      
    point, weit = GaussQuad(npts)                                                 # Retrieve Gauss Quadrature points and weits from GaussQuad subroutine
    dmat        = E/((1-v)**2) * np.array([[1, v, 0],[v, 1, 0],[0, 0,((1-v)/2)]]) # Elasticity tensor, based on plane stress assumption

    ke = np.zeros([8,8])    # Element stiffness matrix initialization
    DdotB = np.zeros([3,8]) # Strain Displacement Matrix dot product with Elasticity Tensor initialization
    count = 0

    for i in range(npts) :                                                        # For number of points in element
        for j in range(npts):                                                     # For number of points in element
            count += 1                                           
            xi  = point[i]                                                        # Assigning point values to Xi
            eta = point[j]                                                        # Assigning point values to Eta
            wti = weit[i]                                                         
            wtj = weit[j]
            sn, dndx, dnde = shapefunDeri(xi, eta)                                # The Shape Function from point i to point j
            ai, detj, I    = ajacob(dndx,dnde,xyel)                               # The Jacobian and its determinant & inverse, from point i to point j
            B  = B_function(sn, dndx, dnde ,xyel)                                 # Strain displacement matrix from point i to point j. 
            ke = ke + np.dot(np.dot(B.T, dmat), B) * detj * wti * wtj             # Elastic stiffness matrix, square matrix
            
            DdotB = np.dot(dmat, B)                                               # Strain Displacement Matrix dot product with Elasticity Tensor
    return ke, DdotB 
