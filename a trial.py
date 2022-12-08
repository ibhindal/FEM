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

# Assign variables to the data imported from MATLAB
global_node_coor = mat['nodecoor']
elemconnect = mat['elemconnect']

nnodes = global_node_coor.shape[0]        # Total number of nodes in mesh  
nelem = elemconnect.shape[0]              # Total number of elements in mesh
elem_node_coor = np.zeros((nelem,4,2))    # Node coordinates of a 4 noded element. Assumes all elements in the mesh have 4 nodes.
                                          # Array storing xy coordinates of the 4 nodes in an element

# Loop through the global_node_coor array and populate elem_node_coor with the xy coordinates of each node
for j in range(nelem):                                                    # for each element
    for i in range(4):                                                    # for each node in each element
        for k in range(2):                                                # for x and y coordinates 
            node = elemconnect[j,i]                                       # Get the elements i(th) connection node
            node_x, node_y = global_node_coor[node, :]                    # Get the global co-ordinates of the node
            elem_node_coor[j,i,0],elem_node_coor[j,i,1] = node_x, node_y  # Store the nodes co-ordinates with respect to the element

#initialise variables and constants 
dofpn    = 2                             # dof per node, x co-or and y co-or
npe      = npts = 4                      # nodes per element  # number of points per element
ndof     = nnodes*dofpn                  # total number of dof
node_coor = []                           # node coordinate, holding value for the current node                             

# Young's Modulus and Poisson's Ratio for different materials 
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
mat_prop_dict = {"Head" : [E_head, nu_head], "Stem" : [E_stem, nu_stem], "Cortical":[E_cortical,nu_cortical], "Trebecular":[E_trebecular,nu_trebecular], "Marrow":[E_marrow,nu_marrow]}
material      = {0 : "Head", 1 : "Stem", 2 : "Cortical", 3 : "Trebecular", 4 : "Marrow"}
nelmmat       = []                      # number of elements per material

# initializes sizes
con_mat        = np.zeros((nelem,4))                # connectivity matrix, matrix to associate nodes to belong to an element
MaterialforElm = np.zeros(nelem)                    # materials for each element, follows key of Matrial dictionary

# Connectivity Matrix and Element Material Matrix population                           
con_mat[i] = elemconnect[i][:4]                 # for each element
for j in range(4) :                             # for each node in the element
    con_mat[i][j] = np.array(elemconnect[i][j]) # connectivity matrix #check this
    #nodeCoor = globalNodeCoor[i]                    

#Gauss Quadrature points & weights of a Corner (npts=4)
# Gauss Quadrature points & weights of a Corner (npts=4)
p, w = GaussQuad.GaussQuad(4)  # the Gauss Quadrature points & weights
qpt = p  # point, vector quadrature points (npts)
qwt = w  # weight, vector quadrature weights (npts)
nquad = qpt.shape[0]  # number of quadrature points

# Shape functions and their derivatives
sn, dsfdx, dsfde = shapefunDeri(qpt, qpt)  # calculate shape functions and their derivatives

# Calculate Jacobian matrix and its determinant and inverse
jacobmat, detj, I = JacobianMat(dsfdx, dsfde, xyel)

# Calculate the element stiffness matrix
ke = kecalc(sn, dsfdx, dsfde, detj, I, E, v)

# Initialize global stiffness matrix
Kg = np.zeros((ndof, ndof))

# Assemble the global stiffness matrix
for i in range(nelem):
    # Connectivity matrix for element i
    con_mat = elemconnect[i][:4]

    # Node coordinates for element i
    xyel = np.array([globalNodeCoor[j] for j in con_mat])

    # Calculate element stiffness matrix and material property tensor
    ke, D = kecalc(sn, dsfdx, dsfde, detj, I, E[i], v[i])

    # Assemble global stiffness matrix
    Kg += assemble(con_mat, ke)
    # calculates the element siffness matrix
                                                    # ke: elastic stiffness matrix 
                                                    # D : the D matrix in tensor form, stress = D x strain
    #Kg = assembleSys(Kg,ke,con_matrix)          #geometric (initial stress) stiffness matrix


plt.plot(Kg)          

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


#meshvtk.cell_data = con_matrix.cell_data
meshio.write("linear_Mesh.vtk", meshvtk)
