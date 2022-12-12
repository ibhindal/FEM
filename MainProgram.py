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
#from matplotlib.tri import Triangulation

#Load in mesh
Meshfilename = 'data.mat'
mat = spio.loadmat(Meshfilename, squeeze_me=(True)) 
print("\n\n\nfile Loaded")

# Assigning variables to the data imported from MATLAB
globalNodeCoor = mat['nodecoor']
elemconnect    = mat['elemconnect'] - 1
nnodes         = globalNodeCoor.shape[0]  # Total number of nodes in mesh  
nelem          = elemconnect.shape[0]     # Total number of elements in mesh
elemNodeCoor   = np.zeros((nelem,4,2))    # Node coordinates of a 4 noded element. Assumes all elements in the mesh have 4 nodes.
print("data imported")                    # Array storing xy coordinates of the 4 nodes in an element

# Plot the Mesh and output to user
colour = {0 : 'b', 1 : 'r', 2 : 'g', 3 : 'c', 4 : 'y'}
nx, ny = np.zeros(4), np.zeros(4) 
plt.plot(globalNodeCoor[:, 0],globalNodeCoor[:, 1],'ro', markersize=0.5) # plots the points 
for i in range(nelem):
    color = colour[elemconnect[i,4]]
    for j in range(4):
        nx[j] = globalNodeCoor[elemconnect[i,j], 0]            
        ny[j] = globalNodeCoor[elemconnect[i,j], 1]
        if j > 0 :
            plt.plot([nx[j-1],nx[j]],[ny[j-1],ny[j]], color)   #plot a line of the quadrangular element
        if j == 3 :
            plt.plot([nx[0],nx[3]],[ny[0],ny[3]], color)       #plot the line of the quadrangular element from the first node to the last

#plt.plot(nx,ny) # trying to plot the mesh with element lines 
plt.plot(globalNodeCoor[:, 0],globalNodeCoor[:, 1],'ro', markersize=0.5) # plots the points 
plt.title("Mesh of implant")
plt.show()

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
print("Material properties determined")

# initializes variables     
shapefundx = shapefunde = jacob = []                # shape function (dx), shape function (de), Jacobian (will store the jacibian, inverse jacobian and determinate)        
Kg         = sp.sparse.csr_matrix((ndof,ndof))      # geometric (initial stress) stiffness matrix
StressDB = []

for i in range(nelem) :                             # for each element                  
    MatNo = elemconnect[i,4]    
    E = Mat_Prop_Dict[Material[MatNo]][0]           # Youngs modulus of current material
    v = Mat_Prop_Dict[Material[MatNo]][1]           # Poissions ratio of current material
    nodeXY = globalNodeCoor[elemconnect[i,0:4].T,:] # finds the node coordinates of current element
    ke,DdotB = kecalc(npts,E,v,nodeXY)                    # calculates the element siffness matrix and D dot B (for stress)
    
    StressDB.append(DdotB)                                  # Saves D dot B
                                                    # ke: elastic stiffness matrix 
    Kg = assembleSys(Kg,ke,elemconnect[i,0:4])      # geometric (initial stress) stiffness matrix

plt.spy(Kg, markersize=0.1)                         # plots the mattrix showing sparcity?
plt.title("Stiffness matrices")
plt.show()
print('Stiffness Matrices calculated')

##### Dirichlet BC (built-in edge y=(min_y)) ####### 
min_y   = np.min(globalNodeCoor[:,1])
nodesbc = np.where(globalNodeCoor[:,1] == min_y)           # find the nodes on the bottom edge y=(min_y)
dofbc   = np.c_[2*nodesbc[0], 2*nodesbc[0]+1].flatten()    #
K_bc    = DirichletBC(Kg,dofbc)                            # system matrix after boundary conditions

# Plot the sparsity of the system stiffness matrix
plt.spy(K_bc, markersize=0.1)
plt.show()
plt.title("Sparsity of system stiffness matrix")
flnmfig = "sparsity_K_beforeBC.png"
plt.savefig(flnmfig)

a = float('-inf')                                   # holds the index of the top node of the trebecular material
for i in range(nelem):                              # for each element
    if (3 == elemconnect[i,4]):                     # check element is the correct material
        try:
            if (globalNodeCoor[a,1] < globalNodeCoor[i,1]): # check element is higher than the last 
                a = i                                       # update a to have the new highest found element
        except:                                         # first trebecular element found 
            a = i
            continue
topnodeTr        = 2*a+1                            # index of the top node of the trebecular bone #Isaac: i think *2 for force fx and fy, +1 is for the y component
topnodeHead      = 2*np.max(elemconnect[:,1]) + 1   # top node of the implant head == top node
F                = np.zeros(K_bc.shape[0])          # Global Force vector  
F[topnodeTr]     =  1.0                             # upward force at trebecular
F[topnodeHead]   = -1.0                             # downward force at the head
F[topnodeHead-1] =  1.0                             # force in x direction at the head
print("Forces and boundary conditions determined")

u = sp.sparse.linalg.spsolve(K_bc, F)               # Calculate the force matrix then we need to plot u #isaac:What?
print("Deformation solved")                         # Calulates the displacement/ deformation



# plot the deformation, u on the mesh
EF = 1                                              # Exageration Factor
nx, ny, ux, uy  = np.zeros(4), np.zeros(4), np.zeros(4), np.zeros(4)
u_x = [num for i, num in enumerate(u) if i % 2 == 0]    # x component deformations for each node
u_y = [num for i, num in enumerate(u) if i % 2 == 1]    # y component deformations for each node
colourlinspace = [0,  1]
#maxfound, minfound = 0,1

ue = []

for i in range(nelem):
    for j in range(4):
        nx[j] = globalNodeCoor[elemconnect[i,j], 0]            
        ny[j] = globalNodeCoor[elemconnect[i,j], 1]
        ux[j] =            u_x[elemconnect[i,j]]
        uy[j] =            u_y[elemconnect[i,j]]
        
        #c     = ux[j] + uy[j]
        #if c > maxfound :
        #    maxfound = c
        #if c < minfound :
        #    minfound = c
        mini = -1.1668924815353183e-11      #minimum deformation found in abisheks mesh with normal force, hardcoded colourmap
        maxi = 3.926118295407239e-09        #maximum deformation found in abisheks mesh with normal force
        ratio = 2 * (ux[j]+uy[j]-mini) / (maxi - mini)
        b = min(1, max(0, (1 - ratio)))
        r = min(1, max(0, (ratio - 1)))
        g = min(1, max(0, 1 - b - r))

        if j > 0 :
            plt.plot([nx[j-1]+ux[j-1]*EF,nx[j]+ux[j]*EF],[ny[j-1]+uy[j-1]*EF,ny[j]+uy[j]*EF], color = (r, g, b))   #plot a line of the quadrangular element
        if j == 3 :
            plt.plot([nx[0]+ux[0]*EF,nx[3]+ux[3]*EF],[ny[0]+uy[0]*EF,ny[3]+uy[3]*EF], color = (r, g, b))           #plot the line of the quadrangular element from the first node to the last
    
    ue.append(np.array([ux[0], uy[0], ux[1], uy[1], ux[2], uy[2], ux[3], uy[3]]))       # Save u values in format needed for stress calculation
                                                                                                                        # set colour scheme to be a measure of deformation (green to red coloours bar), RGB tuple format
#print(maxfound)
#print(minfound)
plt.plot(u_x + globalNodeCoor[:,0], u_y + globalNodeCoor[:,1], 'ro', markersize = 0.5) # plots deformations + initial position
plt.title("Deformation of mesh")
plt.show()

# solve stresses
stress = []

for i in range(nelem):
    ueT = np.transpose(ue[i])       # Transpose ue to become column vector
    stressdb = StressDB[i]
    stress.append(np.dot(stressdb, ueT))    # Dot of DdotB and ue to get stresses
print("Stresses solved")

print("End of Program\n\n\n")
