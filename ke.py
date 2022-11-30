#this is to calculate the ke 
import numpy as np
from GaussQuad import GaussQuad

#kecalc(npts,D,xyel)
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
        #print(xi, eta, wti, wtj)
        arr1=np.array([[eta-1, 1-eta, 1+eta, -1-eta], [xi-1, -1-xi, 1+xi, 1-xi]])
        print(arr1)
        
       # [B,detJ]=JBMat(xyel,sdv)
        #ke = ke + B.'*D*B*detJ*wti*wtj

    #return ke