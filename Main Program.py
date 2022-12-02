#Main Program

import os
import numpy as np
import areaFun
import GaussQuad
import JacobianMat
import shapefunDeri
import strainDisp2D
import assembleSys
import scipy as sp
import scipy.io as spio
#import pandas as pd
from ke import kecalc
import meshio
import matplotlib.pyplot as plt
from assembleSys import assembleSys
from DirichletBC import DirichletBC


#############################################


mat = spio.loadmat('data.mat', squeeze_me=(True))

# Assigning variables to the data imported from MATLAB
globalNodeCoor = mat['nodecoor']
elemconnect = mat['elemconnect']

# Node coordinates of a 4 noded element. Assumes all elements in the mesh have 4 nodes.
nnodes = globalNodeCoor.shape[0] # Total number of nodes in mesh  
nelem = int(nnodes/4) # Total number of elements in mesh
elemNodeCoor = np.zeros((nelem,4,2)) # Array storing xy coordinates of the 4 nodes in an element

# Loops through the globalNodeCoor array and fills elemNodeCoor with the xy coordinates of the nodes
x = -1
for j in range(nelem):
    for i in range(4):
        x += 1
        for k in range(2):
            elemNodeCoor[j,i,k] = globalNodeCoor[x,k]


#mesh.msh[quads]
con_mat=np.zeros((nelem,4))
nodeCoor=[]
ndof =[]

dofpn = 2                               # dof per node
npe =4                                  # nodes per element     
ndof = nnodes*dofpn                     # total number of dof
npts=4                              
nelmmat= 1                             # number of elements per material
#element info
D=np.zeros(nelem)

for i in range(nelem):
    D[i]=elemconnect[i][4]                         #list of materials for each element
    for j in range(4) :
        con_mat[i][j] = np.array(elemconnect[i][j] )# connectivity matrix
    nodeCoor = globalNodeCoor[i]                    # node coordinate matrix

p,w = GaussQuad.GaussQuad(2) #getting values from GaussQuad
qpt=p 
qwt=w 

xyel= np.zeros([4,2]) # update this variable
nquad = qpt.shape[0]  # a matric for the coordinates from gaus quad
ke = np.zeros([8,8]) #defining the stiffness matrix

Kg=np.zeros((ndof,ndof)) #global stiffness matrix

for c in range(nquad):# what we need to do is extract each line from elem connect and input each number as an index into nodecoor then assign it to the
    for w in range(4):
        a = elemconnect[c][w]
        bx,by = globalNodeCoor[a-1] 
    
        xyel[w,0],xyel[w,1] = bx, by # probably would be an idea to stick this in an array instead of putting the whole program in a for loop. i mean either  works...

    for ii in range(nquad) :
        for jj in range(nquad) :
            xi = qpt[ii]                                                    # Stolen from abishek to calculate the
            eta = qpt[jj]                                                   # Shapefunction and jacobian
            sn,dsfdx, dsfde = shapefunDeri.shapefunDeri(xi,eta)
            jacobmat, detj, I = JacobianMat.ajacob(dsfdx, dsfde, xyel)
                

    shapefundx,shapefunde = [],[]
    jacob = []
    count = -1
    

    for d in range(4): #this is wrong. need code for all the elements of each material type.- for all elements where D=1, when thats done all where D=2 ...
        for e in range(nelmmat):
        
            ElemDistMat= np.zeros([8,ndof]) #Element distrribution matrix
            
            ke,dee=kecalc(npts,d,xyel) 
            con_matrix =con_mat[e,:]
            Kg = assembleSys(Kg,ke,con_matrix)   


plt.plot(Kg)          

##### Dirichlet BC (built-in edge y=0) #######
nodesbc = np.where(nodeCoor[:,1] == 0)[0]   # find the nodes on edge y=0
dofbc = np.c_[3*nodesbc, 3*nodesbc+1, 3*nodesbc+2].flatten()
K_bc = DirichletBC(Kg,dofbc)    # system matrix after boundary conditions

# Plot the sparsity of the system stiffness matrix
plt.spy(K_bc, markersize=0.1)
plt.show()
flnmfig = "sparsity_K_beforeBC.png"
plt.savefig(flnmfig)


'''
 
print(str((E/(1-nu**2))))
D = E/(1-nu**2)*np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]])
strainVec = np.dot(bfun,uxy.flatten()) 
stressVec = np.dot(D,strainVec) 
'''


################################################################################
######## CALCULATING THE EIGENVALUES AND VECTORS OF THE SYSTEM MATRIX ##########
################################################################################
nmodes = 6      # number of modes to calculate~did he mean nodes or is modes correct?
eigVal, eigVec = sp.sparse.linalg.eigsh(K_bc, k=nmodes, which='SM')



# write the mesh in vtk format for visualization in Paraview (for e.g.)
meshvtk = meshio.Mesh(nodeCoor, con_matrix)

for ii in range(nmodes) :
    nm = "eVec%d" %(ii+1)
    meshvtk.point_data[nm] = np.c_[np.zeros_like(eigVec[::3,ii]),
                np.zeros_like(eigVec[::3,ii]), eigVec[::3,ii]]


#meshvtk.cell_data = con_matrix.cell_data
meshio.write("linear_Mesh.vtk", meshvtk)
