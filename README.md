# FEM
This is the code respiritory for group 4 FEA module. 

This program displays, on a graph, the stress field and deformations in a 2D representation of a hip replacement. 
This program takes a 2D quadrangular mesh with up to 5 materials, applies loads to the top nodes of the implant head and cortical bone and applies a built in boundary condition to the bottom edge. The mesh is then solved for deformation and stress using predetermined material properties and forces. 

This code has been tested in Python 3.10.8

Required Modules and their tested versions
- numpy
- panda
- scipy
- matplotlib

Input files used to test program
- data.mat
- GoodHipPos.xls with GoodHipQuads.xls
- GoodHipPos2.xls with GoodHipQuads2.xls

The input mesh files are accepted in 2 formats: 
  a .mat file containing a strut that contains two matricies under the name of 'nodecoor' and 'elemconnect'. 
  2 .xls files containing a table of the equivalent 'nodecoor' and 'elemconnect matricies

nodecoor matrix should contain the global coordinates of each node in x, y and/or z components in their respective columns.  
elemconnect matrix should contain the element connectivity matrix i.e. the index ID of each node that each element is comprised of, seperated into columns, wwith an aditional column specifiying material. 

The materials should be identified under the following key
0, 61 : Implant Head
1, 62 : Implant Stem
2, 63 : Cortical Bone
3, 64 : Trebecular Bone
4, 65 : Bone Marrow

The forces are hardcoded point loads, placed at the top of the cortical bone and the top of the implant head
The boundary conditiion is a hardcoded and gives each node with the minimum y co-ordinate a built in boundary condition 

The program will supply output of;
  minimum and maximum deformation and stress to the terminal
  a graphical plot of the mesh with colours representing each material
  a graphical plot of the mesh with deformation with colour map scaled to deformation
  a graphical plot of the mesh with deformation with colour map scaled to stress
  a plot of the population of the global stiffness matrix
