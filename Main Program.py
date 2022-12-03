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
elemconnect = mat['elemconnect']

nnodes = globalNodeCoor.shape[0]        # Total number of nodes in mesh  
nelem = int(nnodes/4)                   # Total number of elements in mesh
elemNodeCoor = np.zeros((nelem,4,2))    # Node coordinates of a 4 noded element. Assumes all elements in the mesh have 4 nodes.
                                        # Array storing xy coordinates of the 4 nodes in an element

# Loops through the globalNodeCoor array and populates elemNodeCoor with the xy coordinates of each node
for j in range(nelem):                  #for each element
    for i in range(4):                  #for each node in each element
        for k in range(2):              #for x and y coordinates 
            elemNodeCoor[j,i,k] = globalNodeCoor[j*4 + i,k]

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
E_trebecular  = 1e9
nu_trebecular = 0.3
E_marrow      = 3e8
nu_marrow     = 0.45
Mat_Prop_Dict = {"Head" : [E_head, nu_head], "Stem" : [E_stem, nu_stem], "Cortical":[E_cortical,nu_cortical], "Trebecular":[E_trebecular,nu_trebecular], "Marrow":[E_marrow,nu_marrow]}
Material      = {0 : "Head", 1 : "Stem", 2 : "Cortical", 3 : "Trebecular", 4 : "Marrow"}
nelmmat       = []                     # number of elements per material

# initializes sizes
con_mat        = np.zeros((nelem,4))                # connectivity matrix, matrix to associate nodes to belong to an element
MaterialforElm = np.zeros(nelem)                    # materials for each element, follows key of Matrial dictionary

# Connectivity Matrix and Element Material Matrix population
for i in range(nelem):                              # for each element
    MaterialforElm[i] = elemconnect[i][4]           
    for j in range(4) :                             # for each node in the element
        con_mat[i][j] = np.array(elemconnect[i][j]) # connectivity matrix #check this
    #nodeCoor = globalNodeCoor[i]                    

#Gauss Quadrature points & weights of a Corner (npts=4)
p,w = GaussQuad.GaussQuad(4)                        # the Gauss Quadrature points & weights
qpt = p                                             # point, vector quadrature points (npts)
qwt = w                                             # weight, vector quadrature weights (npts)

"""update this variable. Replace with xyel = elemNodeCoor[j,:,:]? Where j is the element in reference"""
nquad = qpt.shape[0]                                # Shape of the points returned from GuassQuad, 1x2 for npts=2, 1x4 for npts=4
ke    = np.zeros([8,8])                             # local stiffness matrix
Kg    = np.zeros((ndof,ndof))                       # global stiffness matrix
xyel  = np.zeros([4,2])                             # x y coordinnate matrix for the current element, node coordinate matrix (nne X 2)
#xyels = np.zeros([nelem,4,2])                      # x y coordinates for all elements
xyels=[]

for c in range(nquad):                              # for each points' connection
    for w in range(4):                              # for each point
        a = elemconnect[c][w]                       # get the connectivity of that point
        bx, by = globalNodeCoor[a]                  # get the co-ordinates of the point
        xyel[w,0],xyel[w,1] = bx, by                # store the co-ordinates of the point
    #xyels[nelem, :, :] = xyel                      # array of xyel
    xyels.append(xyel)

for ii in range(nquad) :                            # for each points' connection
    for jj in range(nquad) :                        # for each points' connection
        xi = qpt[ii]                                # xi, eta : local parametric coordinates
        eta = qpt[jj]                            
        sn,dsfdx, dsfde = shapefunDeri.shapefunDeri(xi,eta) #Calculate the derivative of the shape function for a 2D linear quadrangular element with respect to local parametric coordinates xi and eta
                                                            #  sn   : 
                                                            #  dsfdx : vector of shape function derivates w.r.t xi dsfde 
                                                            #  dsfde : vector of shape function derivates w.r.t eta 
        jacobmat, detj, I = JacobianMat.ajacob(dsfdx, dsfde, xyel) # calculates Jacobian at the local coordinate point (xi,eta)
                                                                    # jacobmat : the Jacobian matrix (2 by 2) 
                                                                    # detj     : the determinant of the Jacobian matrix

#initilizes variables            
shapefundx,shapefunde = [],[]
jacob = []
count = -1
D     = np.zeros(nelem)                             # D matrix, in tensor form, a (1 x nelem) matrix   

for MatNo in range(5):                              # for each material, removed nelmmat as it is equal to 1
    for e in range(nelem) :                         # for each element                  
        #ElemDistMat = np.zeros([8,ndof])           # Element distribution matrix
        E = Mat_Prop_Dict[Material[MatNo]][0]       # Youngs modulus of current material
        v = Mat_Prop_Dict[Material[MatNo]][1]       # Poissions ration of current material
        ke, D[e] = kecalc(npts,E,v,xyels[e])        # calculates the element siffness matrix
                                                    # ke: 
                                                    # D :
        con_matrix = con_mat[e,:]
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
stressVec = np.dot(D,strainVec) #jjj
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
