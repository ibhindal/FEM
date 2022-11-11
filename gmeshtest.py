import meshio
import os
os.chdir("C:\\Users\\Ibrahim\\Desktop\\fem")

filenm = 'plateHole.msh'
mesh = meshio.read(filenm) # optionally specify file_format
