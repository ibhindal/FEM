# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 12:12:35 2022

@author: Ben
"""
import scipy.io as spio
import numpy as np
import meshio

#filenm = 'plateHole_linear.msh'
#mesh = meshio.read(filenm)

#nodecoor = mesh.points

#Load in MATLAB .mat file containing node coordinates and element connectivity matrices
mat = spio.loadmat('data.mat', squeeze_me=(True))

#Assigning variables to the data imported from MATLAB
globalNodeCoor = mat['nodecoor']
elemconnect = mat['elemconnect']

#Node coordinates of a 4 noded element. Assumes all elements in the mesh are quadrangular.
nne = globalNodeCoor.shape[0] #Number of nodes 
ne = int(nne/4) #Number of elements
elemNodeCoor = np.zeros((ne,4,2))

x = -1

for j in range(ne):
    for i in range(4):
        x = x + 1
        for k in range(2):
            elemNodeCoor[j,i,k] = globalNodeCoor[x,k]
   
    
            
       



