#this is to calculate the ke 
import numpy as np
from GaussQuad import GaussQuad
from shapefunDeri import shapefunDeri
from JacobianMat import ajacob
from B_function import B_function

#kecalc(npts,D,xyel)
xyel=np.array([8,2])
#xyel=np.array([[0,0] ,[0,1], [1,0], [1,1],[1,0], [2,0], [2,1],[1,1]])
xyel=np.array([[0,0] ,[0,1], [1,0], [1,1]])
npts=4

point, weit = GaussQuad(npts)

ke= np.array([])

count=0

for ii in range(npts) :
    for jj in range(npts):
        count+=1
        #print(count)
        xi = point[ii]
        eta = point[jj]
        wti = weit[ii]
        wtj = weit[jj]
        sn, dndx, dnde = shapefunDeri(xi, eta)
        ai, detj, I = ajacob(dndx,dnde,xyel)
        B= B_function(sn, dndx, dnde ,xyel)
        
        ke = ke + B.T*D*B*detJ*wti*wtj
print('wow')
    #return ke