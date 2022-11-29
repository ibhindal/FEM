#trial1
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import scipy as sp
from scipy import sparse, io
from scipy.sparse import linalg

import meshio



# read the file
filenm = 'plateHole_linear.msh'
mesh = meshio.read(filenm)  # optionally specify file_format


dofpn = 3                               # dof per node
nelem = mesh.cells['quad'].shape[0]     # number of elements
npe = mesh.cells['quad'].shape[1]       # nodes per element
con_mat = mesh.cells['quad']            # connectivity matrix
nodeCoor = mesh.points                  # node coordinate matrix
nnodes = mesh.points.shape[0]           # number of nodes
ndof = nnodes*dofpn                     # total number of dof


######### SECTION FOR WRITING VTK FILE FOR PARAVIEW VISUALIZATION ##########
# write the mesh in vtk format for visualization in Paraview (for e.g.)
meshvtk = meshio.Mesh(mesh.points, mesh.cells)

# PLOT distance from the origin at each node
meshvtk.point_data = {'distR': (mesh.points**2).sum(axis=1),
    'distRSin': np.sin(mesh.points**2).sum(axis=1) }

meshvtk.cell_data = mesh.cell_data
meshio.write("plateHole_linear.vtk", meshvtk)



############# LOAD ELEMENT STIFFNESS MATRIX #############

flnm = "ke.mat"
fl_list = sp.io.loadmat(flnm)
Kel = fl_list['ke']

K = sp.sparse.csr_matrix((ndof,ndof))


# write a subroutine assembleSys.py to assemble the matrix

plt.spy(K, markersize=0.1)
plt.show()



##### IMPLEMENT Dirichlet BC #######






