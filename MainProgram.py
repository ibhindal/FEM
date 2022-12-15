# Main Program
import numpy as np
import pandas as pd
import scipy as sp
import scipy.io as spio
import matplotlib.pyplot as plt
import assembleSys
from ke import kecalc
from assembleSys import assembleSys
from DirichletBC import DirichletBC
#from matplotlib.tri import Triangulation

##### Pythagorous a**2 + b**2 = c**2, used for x & y nodal stress/deformation calculations #####
def PythagDist(i):
    if len(i) == 2:
        return (sum([i[0]**2, i[1]**2])**0.5)
    else:
        return (sum([i[0]**2, i[1]**2, i[2]**2])**0.5)

##### Importing Data from MATLAB #####
try:
    Meshfilename = 'data.mat'                   # .mat file containing mesh data
    mat = spio.loadmat(Meshfilename, squeeze_me=(True)) 
    globalNodeCoor = mat['nodecoor']            # Assigning Node Coordinates to a new variable
    elemconnect    = mat['elemconnect'] - 1     # Assigning Element Connectivity matrix to a new variable
except:
    data = pd.read_excel("GoodHipPos2.xlsx", sheet_name=0)      # run this code to test group made meshes
    globalNodeCoor = data.to_numpy()
    data = pd.read_excel("GoodHipQuads2.xlsx", sheet_name=0)
    elemconnect = data.to_numpy()
    for i in range(elemconnect.shape[0]):
        elemconnect[i,0] = elemconnect[i,0] -1 
        elemconnect[i,1] = elemconnect[i,1] -1 
        elemconnect[i,2] = elemconnect[i,2] -1 
        elemconnect[i,3] = elemconnect[i,3] -1 
print("\n\n\nfile Loaded")    

nnodes         = globalNodeCoor.shape[0]  # Total number of nodes in mesh  
nelem          = elemconnect.shape[0]     # Total number of elements in mesh
elemNodeCoor   = np.zeros((nelem,4,2))    # Node coordinates of a 4 noded element. Assumes all elements in the mesh have 4 nodes.
print("data imported")                    # Array storing xy coordinates of the 4 nodes in an element

##### Plot Mesh and Output To User #####
colour_dict = {0 : 'b', 1 : 'r', 2 : 'g', 3 : 'c', 4 : 'y', 61: 'b', 62: 'r', 63: 'g', 64: 'c', 65: 'y'} # Defines a diffrent colour for each material
nx, ny = np.zeros(4), np.zeros(4) 
plt.plot(globalNodeCoor[:, 0], globalNodeCoor[:, 1], 'bo', markersize=0.5) # Plots all the nodes from the mesh
for i in range(nelem):                                                     # For each element 
    color = colour_dict[elemconnect[i,4]]                                  # Assigning colour to each node depending on material type 
    for j in range(4):
        nx[j] = globalNodeCoor[elemconnect[i,j], 0] # Finding coresponding x coordinates of the node            
        ny[j] = globalNodeCoor[elemconnect[i,j], 1] # Finding corresponding y coordinate of the node
        if j > 0 :
            plt.plot([nx[j-1],nx[j]],[ny[j-1],ny[j]], color = color, linewidth=0.5)   # Plot a line of the quadrangular element
        if j == 3 :
            plt.plot([nx[0],nx[3]],[ny[0],ny[3]], color = color, linewidth=0.5)       # Plot the line of the quadrangular element from the first node to the last

plt.title("Mesh of implant and bone")
plt.show()

##### Initialising Variables and Constants ##### 
dofpn    = 2                            # Degrees of freedom per node, x co-or and y co-or
npts     = 4                            # Number of points per element
ndof     = nnodes*dofpn                 # Total number of dof
nodeCoor = []                           # Node coordinate, holding value for the current node                             

##### Material Properties #####
"""
E = Young's Modulus
nu = Poisson's Ratio

"""
#Default - Co-Cr_Mo alloy
E_head        = 2.1e11
nu_head       = 0.3
#Default - Ti-6 Al-4V alloy
E_stem        = 1.14e11
nu_stem       = 0.3
#Cortical Bone
E_cortical    = 1.6e10
nu_cortical   = 0.3
#Trebacular Bone
E_trebecular  = 1.0e9
nu_trebecular = 0.3
#Bone Marrow
E_marrow      = 3.0e8
nu_marrow     = 0.45
Mat_Prop_Dict = {"Head" : [E_head, nu_head], "Stem" : [E_stem, nu_stem], "Cortical":[E_cortical,nu_cortical], "Trebecular":[E_trebecular,nu_trebecular], "Marrow":[E_marrow,nu_marrow]} # Dictionary containing material properties for the 5 different material types
Material      = {0 : "Head", 1 : "Stem", 2 : "Cortical", 3 : "Trebecular", 4 : "Marrow",61: "Head", 62: "Stem", 63: "Cortical", 64: "Trebecular", 65: "Marrow"} # Assinging a number to each part
print("Material properties determined")

##### Initializing Variables #####   
shapefundx = shapefunde = jacob = []                # Shape function (dx), shape function (de), Jacobian (will store the jacibian, inverse jacobian and determinate)        
Kg         = sp.sparse.csr_matrix((ndof,ndof))      # Geometric (initial stress) stiffness matrix
StressDB   = np.zeros([nelem, 3, 8])                # Stress matrix initalization

##### Calculating global stiffness matrix #####
for i in range(nelem) :                             # For each element                  
    MatNo = elemconnect[i,4]    
    E = Mat_Prop_Dict[Material[MatNo]][0]           # Youngs modulus of current material
    v = Mat_Prop_Dict[Material[MatNo]][1]           # Poissions ratio of current material
    nodeXY = globalNodeCoor[elemconnect[i,0:4].T,:] # Finds the node coordinates of current element
    ke, DDotB = kecalc(npts,E,v,nodeXY)             # Calculates the element siffness matrix and D dot B (for stress)
    StressDB[i, :, :] = DDotB                       # Strain Displacement Matrix dot product with Elasticity Tensor for each element                                 
    Kg = assembleSys(Kg,ke,elemconnect[i,0:4])      # Geometric (initial stress) stiffness matrix

##### System stiffness matrix sparsity plot #####
plt.spy(Kg, markersize=0.1)                         # creates a plot of all non zero values using the equivalent, non sparce form, indecies                      
plt.title("Stiffness matrices")
plt.show()
print('Stiffness Matrices calculated')

##### Dirichlet Boundary Condition (built-in edge y=(min_y)) ####### 
min_y   = np.min(globalNodeCoor[:,1])
nodesbc = np.where(globalNodeCoor[:,1] == min_y)           # Find the nodes on the bottom edge y=(min_y)
dofbc   = np.c_[2*nodesbc[0], 2*nodesbc[0]+1].flatten()    
K_bc    = DirichletBC(Kg,dofbc)                            # System matrix after boundary conditions, sets the bottom edge as built in

##### System stiffness matrix sparsity plot ######
plt.spy(K_bc, markersize=0.1)
plt.show()
plt.title("Sparsity of system stiffness matrix")
flnmfig = "sparsity_K_beforeBC.png"
plt.savefig(flnmfig)

##### Determining Forces ##### 
a, a2, a3 = 0, 0, 0                                        # Holds the index of the top node(s) of the trebecular material
first, sec, third = True, True, True
for i in range(nelem):                                     # For each element
    if (64 == elemconnect[i,4] or 2 == elemconnect[i,4]):  # Check element is the correct material (bone) 
        for j in range(4):
            if (globalNodeCoor[a,1] < globalNodeCoor[elemconnect[i,j],1]) or first or sec or third: # Check element is higher than the last 
                a, a2, a3 = elemconnect[i,j], a2, a3                                                # Update a to have the new highest found node
                if first:
                    first = False
                elif sec:
                    sec = False
                elif third:
                    third = False

b, b2, b3 = 0, 0, 0                                 # Holds the index of the top node of the head
for i in range(nnodes):
    if globalNodeCoor[i,1] > globalNodeCoor[b,1]:   # Checks if node y-coordinate is higher than the last
        b = i                                       # Sets b as the index of the top node
    elif globalNodeCoor[i,1] > globalNodeCoor[b2,1]:
        b2 = i
    elif globalNodeCoor[i,1] > globalNodeCoor[b3,1]:
        b3 = i

topnodeTr        = 2*a + 1                          # Index of the top node of the trebecular bone # *2 for force fx and fy, +1 is for the y component
#topnodeTr2       = 2*a2 + 1
#topnodeTr3       = 2*a3 + 1

topnodeHead      = 2*b + 1                          # Top node of the implant head == top node
#topnodeHead2      = 2*b2 + 1
#topnodeHead3      = 2*b3 + 1

F                = np.zeros(K_bc.shape[0])          # Global Force vector  
F[topnodeTr]     = +1607#/3                         # Upward force at trebecular
#F[topnodeTr2]     = +1607/3
#F[topnodeTr3]     = +1607/3 

F[topnodeHead]   = -1607#/3                         # Downward force at the head
F[topnodeHead-1] = +373#/3                          # Force in x direction at the head

#F[topnodeHead2]   = -1607/3
#F[topnodeHead2-1] = +373/3

#F[topnodeHead3]   = -1607/3
#F[topnodeHead3-1] = +373/3

print("Forces and boundary conditions determined")

u = sp.sparse.linalg.spsolve(K_bc, F)                
u_x = [num for i, num in enumerate(u) if i % 2 == 0] # x component deformations for each node
u_y = [num for i, num in enumerate(u) if i % 2 == 1] # y component deformations for each node
print("Deformation solved")                          # Calulates the displacement/ deformation

##### Solving and Plotting Deformation #####
EF = 1                                               # Exageration Factor
nx, ny, ux, uy = np.zeros(4), np.zeros(4), np.zeros(4), np.zeros(4)
DeformationSum = [PythagDist(i) for i in zip(u_x, u_y)] # Sum of deformation of each node
mini = min(DeformationSum)                           # Minimum deformation found in mesh
maxi = max(DeformationSum)                           # Maximum deformation found in mesh

for i in range(nelem):
    for j in range(4):
        nx[j] = globalNodeCoor[elemconnect[i,j], 0]  # Global original x coordinate of node    
        ny[j] = globalNodeCoor[elemconnect[i,j], 1]  # Global original y coordinate of node
        ux[j] =            u_x[elemconnect[i,j]]     # x-direction deformation of node
        uy[j] =            u_y[elemconnect[i,j]]     # y-direction deformation of node
        ratio = 2 * (PythagDist([ux[j],uy[j]])-mini) / (maxi - mini) # For deformation colour scale
        b = min(1, max(0, (1 - ratio)))
        r = min(1, max(0, (ratio - 1)))
        g = min(1, max(0, 1 - b - r))

        if j > 0 :
            plt.plot([nx[j-1]+ux[j-1]*EF,nx[j]+ux[j]*EF],[ny[j-1]+uy[j-1]*EF,ny[j]+uy[j]*EF], color = (r, g, b))   # Plot a line of the quadrangular element
        if j == 3 :
            plt.plot([nx[0]+ux[0]*EF,nx[3]+ux[3]*EF],[ny[0]+uy[0]*EF,ny[3]+uy[3]*EF], color = (r, g, b))           # Plot the line of the quadrangular element from the first node to the last
                                                                                                                   # set colour scheme to be a measure of deformation (green to red coloours bar), RGB tuple format
print('Maximum deformation {}'.format(maxi))
print('Minimum deformation {}'.format(mini))

#plt.plot(globalNodeCoor[:, 0],globalNodeCoor[:, 1], 'bo', markersize = 0.5) # plots points of nodes where they would be witout deformation
plt.plot(u_x    + globalNodeCoor[:, 0], u_y    + globalNodeCoor[:, 1], 'o', markersize = 0.5) # plots each node with deformations + initial position
#plt.plot([u_x[a] + globalNodeCoor[a, 0]], [u_y[a] + globalNodeCoor[a, 1]], 'ro', markersize = 5) # plots a marker for where a force is applied 
#plt.plot([u_x[b] + globalNodeCoor[b, 0]], [u_y[b] + globalNodeCoor[b, 1]], 'ro', markersize = 5) # plots a marker for where b force is applied 

plt.title("Deformation of mesh")
plt.show()

##### Solving and Plotting Stress #####
stress = np.zeros(nelem*3)  # Initialize stress field 
ueT = np.zeros(8)           # Initialize array of x & y displacement for each node in an element
for i in range(nelem):
    ueT = np.transpose([u_x[elemconnect[i,0]], u_y[elemconnect[i,0]],
                        u_x[elemconnect[i,1]], u_y[elemconnect[i,1]],
                        u_x[elemconnect[i,2]], u_y[elemconnect[i,2]],
                        u_x[elemconnect[i,3]], u_y[elemconnect[i,3]]])  # Transpose ue to become column vector
    stress[i*3 : i*3 + 3] = (np.dot(StressDB[i, :, :], ueT))            # Dot of DdotB and ue to get stresses
print("Stresses solved")

EF = 1                                                      # Exageration Factor
st_x  = [num for i, num in enumerate(stress) if i % 3 == 0] # x component deformations for each node
st_y  = [num for i, num in enumerate(stress) if i % 3 == 1] # y component deformations for each node
st_xy = [num for i, num in enumerate(stress) if i % 3 == 2] 
stressSum = [PythagDist(i) for i in zip(st_x, st_y, st_xy)] # Sum of stresses
mini = min(stressSum)                                       # Minimum stress found in mesh 
maxi = max(stressSum)                                       # Maximum stress found in mesh 

for i in range(nelem):                                      # For number of elements
    for j in range(4):                                  
        nx[j] = globalNodeCoor[elemconnect[i,j], 0]         # Global original x coordinate of node    
        ny[j] = globalNodeCoor[elemconnect[i,j], 1]         # Global original y coordinate of node 
        ux[j] =            u_x[elemconnect[i,j]]            # x-direction deformation of node
        uy[j] =            u_y[elemconnect[i,j]]            # y-direction deformation of node
        ratio = 2 * (stressSum[i] - mini) / (maxi - mini)   # For stress colour scale
        b = min(1, max(0, (1 - ratio)))
        r = min(1, max(0, (ratio - 1)))
        g = min(1, max(0, 1 - b - r))

        if j > 0 :
            plt.plot([nx[j-1]+ux[j-1]*EF,nx[j]+ux[j]*EF],[ny[j-1]+uy[j-1]*EF,ny[j]+uy[j]*EF], color = (r, g, b))   # Plot a line of the quadrangular element
        if j == 3 :
            plt.plot([nx[0]+ux[0]*EF,nx[3]+ux[3]*EF],[ny[0]+uy[0]*EF,ny[3]+uy[3]*EF], color = (r, g, b))           # Plot the line of the quadrangular element from the first node to the last
                                                                                                                   # set colour scheme to be a measure of deformation (green to red coloours bar), RGB tuple format
print('Maximum Stress {}'.format(maxi))
print('Minimum Stress {}'.format(mini))
plt.plot(u_x + globalNodeCoor[:,0], u_y + globalNodeCoor[:,1], 'bo', markersize = 0.5) # Plots stress at each node + initial position
plt.title("Stress of mesh")
plt.show()

print("End of Program\n\n\n")