from numpy import *
import scipy as sp
from scipy import sparse

def assembleSys(K,Kel,con_el) :
    
    """
    Inputs:
        k: Empty stiffness matrix in sparse form
        Kel: Element stiffness matrix
        con_el: Element connectivity matrix
    Outputs:
        k: Filled stiffness matrix in sparse form
    """
    
    dofpe     = Kel.shape[0]        # DoF per element
    npe       = con_el.shape[0]     # Nodes per element
    dofpn     = dofpe/npe           # DoF per node
    ndof      = K.shape[0]          # Total number of DoF
    Kelsparse = sp.sparse.csr_matrix(Kel) # Converting stiffness matrix Kel into sparse form

    ix   = con_el                                                        # ?
    Rind = c_[dofpn*ix, dofpn*ix+1].flatten()                            # ?
    Cind = arange(dofpe)                                                 # ?
    indV = ones(dofpe)                                                   # ?
    Emat = sp.sparse.csr_matrix((indV, (Rind,Cind)), shape=(ndof,dofpe)) # ?
    K    = K + Emat.dot(Kelsparse).dot(Emat.T)                           # ?

    return K


