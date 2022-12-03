import numpy as np
import shapefunDeri
import JacobianMat
import strainDisp2D

def B_function(N, dndx, dnde ,nodeXY):
   """
   b function is the strain displacement matrix requuires xi, eta and nodal coords
   Inputs :
        N      : 
        dndx   :
        dnde   :
        nodeXY :
   Returns:
        B : the strain diplacement matrix 
   """
    #N, dndx, dnde = shapefunDeri.shapefunDeri(xi, eta) 

    Jacob, SF = {}, {} #Jacobi
    SF['sf'], SF['dndx'], SF['dnde'] = N, dndx, dnde 
    Jacob['J'], Jacob['detJ'], Jacob['invJ'] = JacobianMat.ajacob(SF['dndx'], 
                                                                  SF['dnde'], nodeXY) 
    B = strainDisp2D.strainDisp2D(SF,nodeXY,Jacob)
    #print(B)
    return B