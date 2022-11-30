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
import assembleSys
#import pandas as pd
from ke import kecalc


#############################################


filenm = 'plateHole_linear.vtk'

mesh = meshio.read(filenm) 
#mesh.msh[quads]
nelem=[],[]
npe=[],[]
con_mat= []
nodeCoor=[]
nnodes =[],[]
ndof =[]

dofpn = 2                               # dof per node
nelem = mesh.cells[2].data.shape[0]     # number of elements
npe = mesh.cells[2].data.shape[1]       # nodes per element
con_mat = mesh.cells[2]                 # connectivity matrix
nodeCoor = mesh.points                  # node coordinate matrix
nnodes = mesh.points.shape[0]           # number of nodes
ndof = nnodes*dofpn                     # total number of dof
npts=4
nelmmat=4                               # number of elements per material
#element info
D=1 #change this to bens function

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
              

#ke=kecalc(npts,D,xyel)

'''''
E = 2.1e11 
nu = 0.3 
print(str((E/(1-nu**2))))
D = E/(1-nu**2)*np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]])
strainVec = np.dot(bfun,uxy.flatten()) 
stressVec = np.dot(D,strainVec) 
'''