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




p,w = GaussQuad.GaussQuad(2)

qpt=p 
qwt=w 
xyel= np.zeros([4,2])

nquad = qpt.shape[0]
ke = np.zeros([8,8])

E_head=2.1e11
nu_head=0.3
E_stem=1.14e11
nu_stem=0.3
E_cortical=1.6e10
nu_cortical=0.3
E_trebecular=1e9
nu_trebecular=0.3
E_marrow=3e8
nu_marrow=0.45

Mat_prop=[[E_head,nu_head],[E_stem,nu_stem],[E_cortical,nu_cortical],[E_trebecular,nu_trebecular],[E_marrow,nu_marrow]]

def dmat(val):
    E=val[0]
    v= val[1]
    dmat=E/((1-v)**2) * np.array([[1, v, 0],
                                [v, 1, 0], 
                                [0, 0,((1-v)/2)]])
    return dmat




De = {}
for i, val in enumerate(Mat_prop):
    De[i] = dmat(val)




for ii in range(nquad) :
    for jj in range(nquad) :
        xi = qpt[ii]
        eta = qpt[jj]
        sn,dsfdx, dsfde = shapefunDeri.shapefunDeri(xi,eta)
        jacobmat, detj, I = JacobianMat.ajacob(dsfdx, dsfde, xyel)
            

shapefundx,shapefunde = [],[]
jacob = []
count = -1

for d in range(4):
    for i in range(2) :
        for j in range(2) :
            count += 1
            shapefund= shapefunDeri.shapefunDeri(p[0], p[1]) 
            dev=shapefund[1][2]
            print(dev)
            #jacob=JacobianMat.ajacob(shapefund[1][count],shapefund[2][count],nodeCoor)
            jacob=JacobianMat.ajacob(shapefund[1],shapefund[2],nodeCoor)
            bfun=B_function.B_function(p[0], p[1],nodeCoor)
            ke=kecalc(npts,D,xyel)

  

    #strainDisp2D.strainDisp2D()
Kg=np.zeros(ndof,ndof) #global stiffness matrix
#ke=kecalc(npts,D,xyel)

'''''
E = 2.1e11 
nu = 0.3 
print(str((E/(1-nu**2))))
D = E/(1-nu**2)*np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]])
strainVec = np.dot(bfun,uxy.flatten()) 
stressVec = np.dot(D,strainVec) 
'''