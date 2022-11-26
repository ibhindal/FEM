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
#import pandas as pd

#############################################


filenm = 'CW_P1_Prothesis.m'

mesh = meshio.read(filenm) 
#mesh.msh[quads]
"""
For Matlab:
flnm = "ke.mat"
fl_list = sp.io.loadmat(flnm)
Kel = fl_list['ke']
"""


p,w = GaussQuad.GaussQuad(2)

qpt=p 
qwt=w 

nquad = qpt.shape[0]
ke = np.zeros([8,8])
for ii in range(nquad) :
    for jj in range(nquad) :
        xi = qpt[ii]
        eta = qpt[jj]
        dsfdx, dsfde = shapefunDeri.shapefunDeri(xi,eta)
        jacobmat, detj = JacobianMat.ajacob(dsfdx, dsfde, xyel)
            

shapefundx,shapefunde = []
jacob = []
count = 1


for i in range(2) :
    for j in range(2) :
        count += 1
        shapefund= shapefunDeri.shapefunDeri(p)
    
        jacob=JacobianMat.jacobmat(shapefundx,shapefunde,count)
    
        bfun=B_function.B_function(xi,eta)

    #ke=ke +np.inv(Bfun) * Ce * bfun
  

    strainDisp2D.strainDisp2D()



