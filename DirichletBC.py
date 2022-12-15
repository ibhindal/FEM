from numpy import *
import scipy as sp
from scipy import sparse



def DirichletBC(K,dof_bc) :
    
    """
    Calculates the system matrix with boundary conditions applied 
    Inputs:
        K: Stiffness matrix
        dof_bc: 
    Outputs:
        K_bc: Stiffness matrix with boundary conditions applied
            
    """

    K0_bc = sp.sparse.csr_matrix(K, copy=True)
    ndofBC = dof_bc.shape[0]
    csr_rows_set_nz_to_val(K0_bc, dof_bc)       # Removing the columns
    K0_bc_int = K0_bc.T.tocsr()                 # Transposing makes a CSR matrix to CSC
    csr_rows_set_nz_to_val(K0_bc_int, dof_bc)
    K_bc = K0_bc_int.T.tocsr()                  # Transposing makes a CSR matrix to CSC
    K_bc.setdiag(K.diagonal())

    return K_bc

def csr_row_set_nz_to_val(csr, row, value=0):
    """Set all nonzero elements (elements currently in the sparsity pattern)
    to the given value. Useful to set to 0 mostly.
    """
    if not isinstance(csr, sp.sparse.csr_matrix):
        raise ValueError('Matrix given must be of CSR format.')
    csr.data[csr.indptr[row]:csr.indptr[row+1]] = value

def csr_rows_set_nz_to_val(csr, rows, value=0):
    for row in rows:
        csr_row_set_nz_to_val(csr, row)
    if value == 0:
        csr.eliminate_zeros()

