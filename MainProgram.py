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
nx = ny = np.zeros(4) 
for i in range(nelem):
    for j in range(4):
        nx[j] = globalNodeCoor[elemconnect[i,j], 0]
        ny[j] = globalNodeCoor[elemconnect[i,j], 1]
        if j > 0 :
            plt.plot([nx[j-1],nx[j]],[ny[j-1],ny[j]])   #plot a line of the quadrangular element
        if j == 3 :
            plt.plot([nx[0],nx[3]],[ny[0],ny[3]])       #plot the line of the quadrangular element from the first node to the last

#plt.plot(nx,ny) # trying to plot the mesh with element lines 
plt.plot(globalNodeCoor[:, 0],globalNodeCoor[:, 1],'ro', markersize=0.5) # plots the points 
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
print("Matrial properties determined")

# initializes variables     
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
min_y   = np.min(globalNodeCoor[:,1])
nodesbc = np.where(globalNodeCoor[:,1] == min_y)           # find the nodes on the bottom edge y=(min_y)
dofbc   = np.c_[2*nodesbc[0], 2*nodesbc[0]+1].flatten()    #
K_bc    = DirichletBC(Kg,dofbc)                            # system matrix after boundary conditions

# Plot the sparsity of the system stiffness matrix
plt.spy(K_bc, markersize=0.1)
plt.show()
flnmfig = "sparsity_K_beforeBC.png"
plt.savefig(flnmfig)

a = float('-inf')                                       # holds the index of the top node of the trebecular material
for i in range(nelem):                                  # for each element
    if (3 == elemconnect[i,4]):                         # check element is the correct material
        try:
            if (globalNodeCoor[a,1] < globalNodeCoor[i,1]): # check element is higher than the last 
                a = i                                       # update a to have the new highest found element
        except:                                         # first trebecular element found 
            a = i
            continue
topnodeTr      = 2*a+1                                  # index of the top node of the trebecular bone #Isaac: i think *2 for force fx and fy, +1 is for the y component
topnodeHead    = 2*np.max(elemconnect[:,1]) + 1         # top node of the implant head == top node
F              = np.zeros(K_bc.shape[0])                # Global Force vector  
F[topnodeTr]   =  1.0                                   # upward force at trebecular
F[topnodeHead] = -1.0                                   # downward force at the head
print("Forces and boundary conditions determined")

u = sp.sparse.linalg.spsolve(K_bc, F)                   # Calculate the force matrix then we need to plot u #isaac:What?
print("Deformation solved")                             # calulates the displacement/ deformation

st = 1
print("Stress solved")  #Stress = D x strain 

# plot the deformation, u
EF = 1                                                  # Exageration Factor
u_x = [num for i, num in enumerate(u) if i % 2 == 0]    # x component deformations
u_y = [num for i, num in enumerate(u) if i % 2 == 1]    # y component deformations
plt.plot(u_x + globalNodeCoor[:,0], u_y + globalNodeCoor[:,1], 'ro', markersize = 0.5) # plots deformations + initial position
plt.show()

# plot the stress, st


# extract necessary data




"""
###### Abishek Lab 3, Code to output plot of deformation #########
from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import meshio
filenm = 'data.mat'
mesh = meshio.read(filenm)                                  # optionally specify file_format
mesh.cell_data                                              # reveals all physical tags
mesh.cell_data['line']['gmsh:physical']                     # following shows the physical tags (1: bc, 2: load)
meshvtk = meshio.Mesh(mesh.points, mesh.cells)              # write the mesh in vtk format for visualization in Paraview (for e.g.)
meshio.write("simpleModMeshio.vtk", meshvtk)
triang = Triangulation(mesh.points[:,0], mesh.points[:,1])
circx, circy = 0.4, 0.0
min_radius = 0.15
triang.set_mask(np.hypot((mesh.points[:,0]-circx)[triang.triangles].mean(axis=1),(mesh.points[:,1]-circy)[triang.triangles].mean(axis=1))< min_radius)
dax, day = 0.0, 0.0
zarbit = np.hypot(mesh.points[:,0] - dax, mesh.points[:,1] - day)
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.use_sticky_edges = False                                 # Enforce the margins, and enlarge them to give room for the vectors.
ax.margins(0.07)
ax.triplot(triang, lw=0.5, color='1.0')                     # The blank mesh
levels = np.arange(0., 1., 0.025)
cmap = cm.get_cmap(name='terrain', lut=None)
ax.tricontourf(triang, zarbit, levels = levels, cmap = cmap)
ax.tricontour (triang, zarbit, levels = levels, colors = ['0.25', '0.5', '0.5', '0.5', '0.5'], linewidths = [1.0, 0.5, 0.5, 0.5, 0.5])
plt.show()
"""
print("End of Program\n\n\n")
