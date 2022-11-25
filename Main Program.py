#Main Program



import os
import numpy as np
import meshio
import areaFun
import B_function
import GaussQuad
import JacobianMat
import shapefunDeri
import strainDisp2D
import mainmesh

#############################################


filenm = 'CW_P1_Prothesis.txt'



mesh = meshio.read(filenm) 


p,w = GaussQuad.GaussQuad(2)

qpt=p 
qwt=w 
nne = xyel.shape[0]
nquad = qpt.shape[0]
ke = np.zeros([8,8])
for ii in range(nquad) :
    for jj in range(nquad) :
        xi = qpt[ii]
        eta = qpt[jj]
        sf, dsfdx, dsfde = shapeFun.elem_shape(xi,eta,nne)
        jacobmat, detj = ajacob(dsfdx, dsfde, xyel)
            

shapefundx,shapefunde = []
jacob=
count = 1
for i in range(2) :
    for j in range(2) :
        count += 1
        shapefund= shapefunDeri.shapefunDeri(p)
    
        jacob=JacobianMat.jacobmat(shapefundx,shapefunde,count)
    
        areaFun.areaFun()

        B_function.B_function()

  

    strainDisp2D.strainDisp2D()



