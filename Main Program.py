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

#############################################


filenm = 'plateHole.msh'
mesh = meshio.read(filenm) # optionally specify file_format
