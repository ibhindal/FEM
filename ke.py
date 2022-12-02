#this is to calculate the ke 
import numpy as np
from GaussQuad import GaussQuad
from shapefunDeri import shapefunDeri
from JacobianMat import ajacob
from B_function import B_function

# Function for calculating the element siffness matrix 

def kecalc(npts,MatNo,xyel):
    
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

    Mat_prop = [[E_head,nu_head],[E_stem,nu_stem],[E_cortical,nu_cortical],[E_trebecular,nu_trebecular],[E_marrow,nu_marrow]]
    
    Mat_Prop_Dict = {"Head" : [E_head, nu_head], "Stem" : [E_stem, nu_stem], "Cortical":[E_cortical,nu_cortical], "Trebecular":[E_trebecular,nu_trebecular], "Marrow":[E_marrow,nu_marrow]}
    Material = {1 : "Head", 2 : "Stem", 3 : "Cortical", 4 : "Trebecular", 5 : "Marrow"}
    
    M = 1
    E = Mat_Prop_Dict[Material[M]][0]
    v = Mat_Prop_Dict[Material[M]][1]

    dmat=E/((1-v)**2) * np.array([[1, v, 0],
                                [v, 1, 0], 
                                [0, 0,((1-v)/2)]])
    
    D = {}
    for i, val in enumerate(Mat_prop):
        D[i] = dmat(val)

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
            ke = ke + B.T.dot(D[MatNo]).dot(B) * detj * wti * wtj
        print('KeCalc - Completed!')
    return ke, D