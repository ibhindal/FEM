import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
import numpy as np
import scipy as sp
from scipy import sparse, io
from scipy.sparse import linalg

import meshio

# user defined subroutines
from assembleSys import assembleSys
from DirichletBC import DirichletBC


################################################################
# read the mesh file
filenm = 'plateHole_linear.vtk'
mesh = meshio.read(filenm)  # optionally specify file_format


dofpn = 3                               # dof per node
nelem = mesh.cells_dict['quad'].shape[0]     # number of elements
npe = mesh.cells_dict['quad'].shape[1]       # nodes per element
con_mat = mesh.cells_dict['quad']            # connectivity matrix
nodeCoor = mesh.points                  # node coordinate matrix
nnodes = mesh.points.shape[0]           # number of nodes
ndof = nnodes*dofpn                     # total number of dof


# write the mesh in vtk format for visualization in Paraview (for e.g.)
meshvtk = meshio.Mesh(mesh.points, mesh.cells)

# distance from the origin
meshvtk.point_data = {'distR': (mesh.points**2).sum(axis=1),
    'distRSin': np.sin(mesh.points**2).sum(axis=1) }

meshvtk.cell_data = mesh.cell_data
meshio.write("plateHole_linear.vtk", meshvtk)

flnm = "ke.mat"
fl_list = sp.io.loadmat(flnm)
Kel = fl_list['ke']

K = sp.sparse.csr_matrix((ndof,ndof))
K = assembleSys(K,Kel,con_mat)

# Plot the sparsity of the system stiffness matrix
plt.spy(K, markersize=0.1)
plt.show()



##### Dirichlet BC (built-in edge y=0) #######
nodesbc = np.where(nodeCoor[:,1] == 0)[0]   # find the nodes on edge y=0
dofbc = np.c_[3*nodesbc, 3*nodesbc+1, 3*nodesbc+2].flatten()
K_bc = DirichletBC(K,dofbc)    # system matrix after boundary conditions


# Plot the sparsity of the system stiffness matrix
plt.spy(K_bc, markersize=0.1)
plt.show()
flnmfig = "sparsity_K_beforeBC.png"
plt.savefig(flnmfig)


# Advanced plot, requires Latex/Miktex
# Setting some Plot parameters
font = {'family' : 'times',
        'weight' : 'normal',
        'size'   : 16}
rc('font', **font)
rc('text', usetex=True)
rc('font', family='times')
#rc('xtick')
#rc('ytick')
rc('legend')


fig = plt.figure()
ax = fig.add_subplot(111)
ax.spy(K_bc, markersize=0.1)
ax.grid(False)
plt.xlabel(r'columns')
plt.ylabel(r'rows')
# plt.show()             # uncomment to show the plot
flnmfig = "sparsity_K_afterBC.png"
plt.savefig(flnmfig)





################################################################################
######## CALCULATING THE EIGENVALUES AND VECTORS OF THE SYSTEM MATRIX ##########
################################################################################
nmodes = 6      # number of modes to calculate
eigVal, eigVec = sp.sparse.linalg.eigsh(K_bc, k=nmodes, which='SM')



# write the mesh in vtk format for visualization in Paraview (for e.g.)
meshvtk = meshio.Mesh(mesh.points, mesh.cells)

for ii in range(nmodes) :
    nm = "eVec%d" %(ii+1)
    meshvtk.point_data[nm] = np.c_[np.zeros_like(eigVec[::3,ii]),
                np.zeros_like(eigVec[::3,ii]), eigVec[::3,ii]]


meshvtk.cell_data = mesh.cell_data
meshio.write("plateHole_linear_U.vtk", meshvtk)





