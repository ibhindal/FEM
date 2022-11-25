#Main Program



import os
import numpy as np
import meshio
import areaFun
import B_function
import GaussQuad
import JacobianMat
import shapeFun
import shapefunDeri
import strainDisp2D
import mainmesh

#############################################


filenm = 'CW_P1_Prothesis.txt'

mainmesh.mainmesh(filenm)

mesh = meshio.read("linear_Mesh.vtk") 
