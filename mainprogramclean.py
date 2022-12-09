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

#Load in mesh
Meshfilename = 'data.mat'
mat = spio.loadmat(Meshfilename, squeeze_me=(True)) 

# Assigning variables to the data imported from MATLAB
globalNodeCoor = mat['nodecoor']
elemconnect = mat['elemconnect'] - 1

nnodes = globalNodeCoor.shape[0]        # Total number of nodes in mesh  
nelem = elemconnect.shape[0]            # Total number of elements in mesh
elemNodeCoor = np.zeros((nelem,4,2))    # Node coordinates of a 4 noded element. Assumes all elements in the mesh have 4 nodes.
                                        # Array storing xy coordinates of the 4 nodes in an element


#initialising variables and constants 
dofpn    = 2                             # dof per node, x co-or and y co-or
npe      = npts = 4                      # nodes per element  # number of points per element
ndof     = nnodes*dofpn                  # total number of dof
nodeCoor = []                            # node coordinate, holding value for the current node                             

# Young's Modulus and Poission's Ratio for different materials 
E_head        = 2.1e11
nu_head       = 0.3
E_stem        = 1.14e11
nu_stem       = 0.3
E_cortical    = 1.6e10
nu_cortical   = 0.3
E_trebecular  = 1.0e9
nu_trebecular = 0.3
E_marrow      = 3.0e8
nu_marrow     = 0.45
Mat_Prop_Dict = {"Head" : [E_head, nu_head], "Stem" : [E_stem, nu_stem], "Cortical":[E_cortical,nu_cortical], "Trebecular":[E_trebecular,nu_trebecular], "Marrow":[E_marrow,nu_marrow]}
Material      = {0 : "Head", 1 : "Stem", 2 : "Cortical", 3 : "Trebecular", 4 : "Marrow"}
nelmmat       = []                     # number of elements per material

# initializes sizes
con_mat        = np.zeros((nelem,4))                # connectivity matrix, matrix to associate nodes to belong to an element
MaterialforElm = np.zeros(nelem)                    # materials for each element, follows key of Matrial dictionary


#initilizes variables            
shapefundx,shapefunde = [],[]
jacob = []
count = -1
D     = np.zeros(nelem)                             # D matrix, in tensor form, a (1 x nelem) matrix 



#         Kg = assembleSys(Kg,ke,con_matrix)          #geometric (initial stress) stiffness matrix
Kg = sp.sparse.csr_matrix((ndof,ndof))
for i in range(nelem) :                         # for each element                  
    #ElemDistMat = np.zeros([8,ndof])           # Element distribution matrix
    MatNo = elemconnect[i,4]
    E = Mat_Prop_Dict[Material[MatNo]][0]       # Youngs modulus of current material
    v = Mat_Prop_Dict[Material[MatNo]][1]       # Poissions ration of current material
    nodeXY = globalNodeCoor[elemconnect[i,0:4].T,:] # finds the node coordinates on the fly
    ke = kecalc(npts,E,v,nodeXY)    # calculates the element siffness matrix
                                                # ke: elastic stiffness matrix 
                                                # D : the D matrix in tensor form, stress = D x strain
    # con_matrix = con_mat[i,:]
    Kg = assembleSys(Kg,ke,elemconnect[i,0:4])          #geometric (initial stress) stiffness matrix

print('Stiffness Matrices - Completed!')

plt.spy(Kg, markersize=0.1)
plt.show()


##### Dirichlet BC (built-in edge y=0) ####### 
#nodesbc = np.where(nodeCoor[:,1] == 0)[0]   # find the nodes on edge y=0
nodesbc = np.where(xyels[:,:,1] == 0)[0]   # find the nodes on edge y=0
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



meshio.write("linear_Mesh.vtk", meshvtk)
