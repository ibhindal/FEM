import numpy as np
import shapeFun
import shapefunDeri
import JacobianMat
import strainDisp2D




def B_function(xi, eta):
    xi, eta = 1.,-1. # user-defined local coordinate 
    N = shapeFun.elem_shape(xi, eta) 
    dndx, dnde = shapefunDeri.shapefunDeri(xi, eta) 
    nodeXY = np.array([[0.0, 0.0], 
                        [100.0, 0.0], 
                        [100.0, 100.0], 
                        [0.0, 100.0]]) 
    uxy = np.array([[0.0, 0.0], 
                    [0.1, 0.0], 
                    [0.1,-0.03], 
                    [0.0,-0.03]]) 

    Jacob, SF = {}, {} 
    SF['sf'], SF['dndx'], SF['dnde'] = N, dndx, dnde 
    Jacob['J'], Jacob['detJ'], Jacob['invJ'] = JacobianMat.ajacob(SF['dndx'], 
                                                                  SF['dnde'], nodeXY) 
    B = strainDisp2D.strainDisp2D(SF,nodeXY,Jacob)
    #print(B)
    return B