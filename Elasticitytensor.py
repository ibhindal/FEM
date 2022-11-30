import numpy as np
if __name__ == "__main__" :
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


