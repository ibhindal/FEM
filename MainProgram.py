#Main Program
import os
import numpy as np
import assembleSys
import scipy as sp
import scipy.io as spio
from ke import kecalc
import matplotlib.pyplot as plt
from assembleSys import assembleSys
from DirichletBC import DirichletBC

#Load in mesh
Meshfilename = 'data.mat'
mat = spio.loadmat(Meshfilename, squeeze_me=(True)) 
print("file Loaded")

# Assigning variables to the data imported from MATLAB
globalNodeCoor = mat['nodecoor']
elemconnect    = mat['elemconnect'] - 1

nnodes       = globalNodeCoor.shape[0]  # Total number of nodes in mesh  
nelem        = elemconnect.shape[0]     # Total number of elements in mesh
elemNodeCoor = np.zeros((nelem,4,2))    # Node coordinates of a 4 noded element. Assumes all elements in the mesh have 4 nodes.
print("data imported")                  # Array storing xy coordinates of the 4 nodes in an element

#initialising variables and constants 
dofpn    = 2                            # degrees of freedom per node, x co-or and y co-or
npts     = 4                            # number of points per element
ndof     = nnodes*dofpn                 # total number of dof
nodeCoor = []                           # node coordinate, holding value for the current node                             

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
print("Matrial properties determined")

# initializes variables
con_mat        = np.zeros((nelem,4))                # connectivity matrix, matrix to associate nodes to belong to an element
MaterialforElm = np.zeros(nelem)                    # materials for each element, follows key of Matrial dictionary           
shapefundx = shapefunde = jacob = []                # shape function (dx), shape function (de), Jacobian (will store the jacibian, inverse jacobian and determinate)
D     = np.zeros(nelem)                             # D matrix, in vector form, a (1 x nelem) matrix          
Kg    = sp.sparse.csr_matrix((ndof,ndof))           # geometric (initial stress) stiffness matrix

for i in range(nelem) :                             # for each element                  
    MatNo = elemconnect[i,4]    
    E = Mat_Prop_Dict[Material[MatNo]][0]           # Youngs modulus of current material
    v = Mat_Prop_Dict[Material[MatNo]][1]           # Poissions ratio of current material
    nodeXY = globalNodeCoor[elemconnect[i,0:4].T,:] # finds the node coordinates of current element
    ke = kecalc(npts,E,v,nodeXY)                    # calculates the element siffness matrix
                                                        # ke: elastic stiffness matrix 
                                                        # D : the D matrix in tensor form, stress = D x strain
    Kg = assembleSys(Kg,ke,elemconnect[i,0:4])      # geometric (initial stress) stiffness matrix

plt.spy(Kg, markersize=0.1)                         # plots the mattrix showing sparcity?
plt.show()
print('Stiffness Matrices calculated')

##### Dirichlet BC (built-in edge y=(min_y)) ####### 
min_y    = np.min(globalNodeCoor[:,1])
print(min_y)
nodesbc = np.where(globalNodeCoor[:,1] == min_y)                        # find the nodes on the bottom edge y=(min_y)
dofbc   = np.c_[2*nodesbc[0], 2*nodesbc[0]+1].flatten()                 #
K_bc    = DirichletBC(Kg,dofbc)                                         # system matrix after boundary conditions

# Plot the sparsity of the system stiffness matrix
plt.spy(K_bc, markersize=0.1)
plt.show()
flnmfig = "sparsity_K_beforeBC.png"
plt.savefig(flnmfig)

a = float('-inf')
for i in range(nelem): # for each element
    if (3 == elemconnect[i,4]): # check element is the correct material
        if (True):  # check element is higher than the last 
            a = i # update a to have the new highest found element
topnodeTr = 2*a+1                                   # replace 171 with the index of the top node of the trebecular bone #Isaac: i think *2 for force fx and fy, +1 is ?
topnodeHead = 2*np.max(elemconnect[:,1]) + 1        # replace 300 with the index of the top node of the implant  head
F  = np.zeros(K_bc.shape[0])                        # is the force 
F[topnodeTr] = 1.0                                  # upward force at trebecular
F[topnodeHead] = -1.0                               # downward force at the head

u = sp.sparse.linalg.spsolve(K_bc, F) #calculate the force matrix then we need to plot u

print("done")
