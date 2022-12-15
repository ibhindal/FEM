import numpy as np
import shapefunDeri
import JacobianMat
import strainDisp2D

def B_function(N, dndx, dnde ,nodeXY):
    
   """
   Prepares the shape function, jacobian and nodal coordinates 
   in the right format to be passed on to strainDisp2D
   Inputs :
        N      : Shape function matrix
        dndx   : Derivative of          w.r.t xi
        dnde   : Derivative of          w.r.t eta
        nodeXY : x y co-ordinates of the element
   Returns:
        B : Strain diplacement matrix, passed on from strainDisp2D
   """
   
   Jacob, SF = {}, {}                                                                              # Jacobian, Shape function dictionaries
   SF['sf'],   SF['dndx'],    SF['dnde']    = N, dndx, dnde                                        # Assigning SF dictionary items to function inputs inputs
   Jacob['J'], Jacob['detJ'], Jacob['invJ'] = JacobianMat.ajacob(SF['dndx'], SF['dnde'], nodeXY)   # Assigning Jacob dictionary items to ajacob function returns
   B = strainDisp2D.strainDisp2D(SF,nodeXY,Jacob)                                                  # Strain displacement matix
   return B