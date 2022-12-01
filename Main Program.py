#Main Program

import os
import numpy as np
import areaFun
import GaussQuad
import JacobianMat
import shapefunDeri
import strainDisp2D
import assembleSys
import scipy.io as spio
#import pandas as pd
from ke import kecalc
import matplotlib.pyplot as plt


#############################################


mat = spio.loadmat('data.mat', squeeze_me=(True))

#Assigning variables to the data imported from MATLAB
globalNodeCoor = mat['nodecoor']
elemconnect = mat['elemconnect']

#Node coordinates of a 4 noded element. Assumes all elements in the mesh are quadrangular.
nnodes = globalNodeCoor.shape[0] #Number of nodes 
nelem = int(nnodes/4) #Number of elements
elemNodeCoor = np.zeros((nelem,4,2))

for j in range(nelem):
    for i in range(4):
        for k in range(2):
            #display(elemNodeCoor[0,i,k])
            #display(globalNodeCoor[i,k])
            elemNodeCoor[j,i,k] = globalNodeCoor[i,k]



#mesh.msh[quads]
con_mat=np.array(con_mat)
nodeCoor=[]
ndof =[]

dofpn = 2                               # dof per node
npe =4                                  # nodes per element     
ndof = nnodes*dofpn                     # total number of dof
npts=4                              
nelmmat=4                               # number of elements per material
#element info


for i in range(nelem):
    D=elemconnect[i][4]                         #list of materials for each element
    #con_mat = np.array(elemconnect[i][0,1,2,3] )# connectivity matrix
    nodeCoor = globalNodeCoor[i]                # node coordinate matrix

p,w = GaussQuad.GaussQuad(2)

qpt=p 
qwt=w 

xyel= np.zeros([4,2])
nquad = qpt.shape[0]
ke = np.zeros([8,8])




for ii in range(nquad) :
    for jj in range(nquad) :
        xi = qpt[ii]
        eta = qpt[jj]
        sn,dsfdx, dsfde = shapefunDeri.shapefunDeri(xi,eta)
        jacobmat, detj, I = JacobianMat.ajacob(dsfdx, dsfde, xyel)
            

shapefundx,shapefunde = [],[]
jacob = []
count = -1
Kg=np.zeros(ndof,ndof) #global stiffness matrix
ElemDistMat= np.zeros([8,ndof]) #Element distribution matrix

for d in range(4):
   for e in range(nelmmat):
    
    ElemDistMat= np.zeros([8,ndof]) #Element distrribution matrix
    ke=kecalc(npts,D,xyel)
    con_matrix =con_mat[e,:]
    Kg = assembleSys(Kg,ke,con_matrix) 


plt.plot(Kg)          
              

#ke=kecalc(npts,D,xyel)

'''''
E = 2.1e11 
nu = 0.3 
print(str((E/(1-nu**2))))
D = E/(1-nu**2)*np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]])
strainVec = np.dot(bfun,uxy.flatten()) 
stressVec = np.dot(D,strainVec) 
'''