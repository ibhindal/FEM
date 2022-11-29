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


filenm = 'plateHole_linear.vtk'

mesh = meshio.read(filenm) 
#mesh.msh[quads]
nelem=[]
npe=[]
con_mat= []
nodeCoor=[]
nnodes =[]
ndof =[]

dofpn = 3                               # dof per node
nelem = mesh.cells[2].data.shape[0]     # number of elements
npe = mesh.cells[2].data.shape[1]       # nodes per element
con_mat = mesh.cells[2]                 # connectivity matrix
nodeCoor = mesh.points                  # node coordinate matrix
nnodes = mesh.points.shape[0]           # number of nodes
ndof = nnodes*dofpn                     # total number of dof

"""
For Matlab:
flnm = "ke.mat"
fl_list = sp.io.loadmat(flnm)
Kel = fl_list['ke']
"""

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


for i in range(2) :
    for j in range(2) :
        count += 1
        shapefund= shapefunDeri.shapefunDeri(p[0], p[1]) 
        jacob=JacobianMat.ajacob(shapefund[0][count],shapefund[1][count],nodeCoor)


        bfun=B_function.B_function(p[0], p[1])

    #ke=ke +np.inv(Bfun) * Ce * bfun
  

    #strainDisp2D.strainDisp2D()



E = 2e5 
nu = 0.3 
print(str((E/(1-nu**2))))
D = E/(1-nu**2)*np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]])
strainVec = np.dot(B,uxy.flatten()) 
stressVec = np.dot(D,strainVec) 
