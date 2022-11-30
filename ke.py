#this is to calculate the ke 
import numpy as np
from GaussQuad import GaussQuad
from shapefunDeri import shapefunDeri
from JacobianMat import ajacob
from B_function import B_function

def kecalc(npts,MatNo,xyel):
    #xyel=np.array([8,2])
    #xyel=np.array([[0,0] ,[0,1], [1,0], [1,1],[1,0], [2,0], [2,1],[1,1]])
    #xyel=np.array([[0,0] ,[0,1], [1,0], [1,1]])
    #npts=4

    point, weit = GaussQuad(npts)



    E_head=2.1e11
    nu_head=0.3
    E_stem=1.14e11
    nu_stem=0.3
    E_cortical=1.6e10
    nu_cortical=0.3
    E_trebecular=1e9
    nu_trebecular=0.3
    E_marrow=3e8
    nu_marrow=0.45

    Mat_prop=[[E_head,nu_head],[E_stem,nu_stem],[E_cortical,nu_cortical],[E_trebecular,nu_trebecular],[E_marrow,nu_marrow]]

    def dmat(val):
        E=val[0]
        v= val[1]
        dmat=E/((1-v)**2) * np.array([[1, v, 0],
                                    [v, 1, 0], 
                                    [0, 0,((1-v)/2)]])
        return dmat



    De = {}
    for i, val in enumerate(Mat_prop):
        De[i] = dmat(val)



    ke= np.zeros([8,8])

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
            xee=De[0]
            print(xee)
            ke = ke + B.T.dot(De[MatNo]).dot(B) * detj * wti * wtj
        print('wow!')
    return ke, De