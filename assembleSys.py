from numpy import *
import scipy as sp

from scipy import sparse


def assembleSys(K,Kel,con_el) :

    dofpe = Kel.shape[0]
    npe = con_el.shape[1]     # nodes per element
    dofpn = dofpe/npe         # dof per node
    ndof = K.shape[0]
    
    Kelsparse = sp.sparse.csr_matrix(Kel)

    nelem = con_el.shape[0]                # number of elements
    for ii in range(nelem) :
        ix = con_el[ii]
        Rind = c_[dofpn*ix, dofpn*ix+1, dofpn*ix+2].flatten()
        Cind = arange(dofpe)
        indV = ones(dofpe)
        Emat = sp.sparse.csr_matrix((indV, (Rind,Cind)), shape=(ndof,dofpe))

        K = K + Emat.dot(Kelsparse).dot(Emat.T)


    return K


