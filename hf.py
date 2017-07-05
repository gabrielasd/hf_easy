import numpy as np
from horton import *


def ints_eval(xyz_file, basis):
    """
    Compute Gaussian integrals.

    Parameters:
    -----------
    xyz_file : molecule coordinate file 
    basis : basis set 
	
    """
     
    mol = IOData.from_file(context.get_fn(xyz_file))
    # Create a Gaussian basis set.
    obasis = get_gobasis(mol.coordinates, mol.numbers, basis)

    # Create a linalg factory.
    lf = DenseLinalgFactory(obasis.nbasis)

    # Compute One-electron Integrals.
    olp = obasis.compute_overlap(lf) # overlap matrix
    olp_m = olp._array
    kin = obasis.compute_kinetic(lf)
    kin_m = kin._array
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    na_m = na._array
    er = obasis.compute_electron_repulsion(lf)
    er_m = er._array

    return {'olp':olp_m, 'kin':kin_m, 'na':na_m, 'er':er_m}


def ortogonalize(olp_matrix):
    """
    Generate overlap ortogonalizing matrix

    Parameters:
    -----------
    olp_matrix : array
                 Overlap Matrix

    """

    diag, unitary_m = np.linalg.eigh(olp_matrix)
    diag_invsqr = np.power(diag, -0.5)
    diag_invsqr_m = diag_invsqr * np.eye(2)
    ortog_m = np.dot(unitary_m, np.dot(diag_invsqr_m, unitary_m.T))

    #print "Diagonalized matrix", X
    #print "Unitary matrix", U
    #XInv = np.linalg.inv(X)
    #Xinsqrt = np.power(Xinv, 0.5)
    #print "Xinsqrt2", Xinsqrt2
    #print Xfin
    #print U.T
    #print

    return {'diag_mtx':diag, 'unitary_mtx':unitary_m, 'ortog_mtx':ortog_m}


def fock_guess(kin_m, na_m, ortog_m):
    # Fock initial guess
    H = kin_m + na_m
    Fzero = np.dot(ortog_m.T, np.dot(H, ortog_m))

    #print "Initial Fock", Fzero
    return Fzero, H


def density_init_guess(Fzero, ortog_m, occ):
    """
    Construct initial guess density matrix
    
    """
	
    # Coefficient matrix initial guess
    # Diagonalize the initial Fock_m
    epsilon, Cp = np.linalg.eigh(Fzero)
    
    # Initial SCF eigenvecto matrix in the original basis
    Czero = np.dot(ortog_m, Cp)
    #print "Czero nos", Czero
    #print "epsilon", epsilon

    # Initial density matrix
    D = np.zeros((2,2))
    for m in range(len(Czero)):
        for v in range(len(Czero)):
            for i in range(len(Czero)):
                D[m,v] += occ[i]*Czero[m,i]*Czero[v,i]
   
    # Otra via de obtener la misma D_inicial
    #D2 = np.dot(occ*Czero, Czero.T)
    
    #print "Density matrix", D
    #print "Density matrix", D2
    return D


def scf_iter(H, D, er_m, E_nuc):
    """
    Perform the SCF iterations
    
    """
	
    # Construct new Fock matrix
	# this includes one and two electron terms
    F = np.zeros((2,2))
	
    for m in range(2):
        for v in range(2):
            F[m, v] = H[m, v]
            for p in range(2):
                for s in range(2):
                    F[m, v] += D[p, s] * (2 * er_m[m,p,v,s] - er_m[m,v,p,s])
    #print "Fock", F

    # Initialize electronic energy
    Energy_0 = 0

    for m in range(2):
        for v in range(2):
            Energy_0 += D[m,v]*(H[m,v] + F[m,v])

    
    E_tot = Energy_0 + E_nuc
    print "Total Energy", Energy_0 + E_nuc

    E_i = 0
    while abs(E_i - E_tot) > 10e-6:
        E_tot = E_i
		
        # Transform Fock matrix to the
		# orthonormal basis
        F = np.dot(Ssqt.T, np.dot(F, Ssqt))

        # Diagonalize Fock matrix
        epsilon, Cp = np.linalg.eigh(F)

        # construct new eigenvector matrix
        C = np.dot(Ssqt, Cp)

        # Form new density matrix
        for m in range(len(C)):
            for v in range(len(C)):
                for i in range(len(C)):
                    D[m,v] += occ[i]*C[m,i]*C[v, i]
        
        # New Fock matrix
        for m in range(2):
            for v in range(2):
                F[m, v] = H[m, v]
                for p in range(2):
                    for s in range(2):
                        F[m, v] += D[p, s] * (2 * er_m[m,p,v,s] - er_m[m,v,p,s])
        # New energy
        for m in range(2):
            for v in range(2):
                E_i += D[m,v]*(H[m,v] + F[m,v])

    print "SCF done. HF Energy is ", E_i


if __name__ == "__main__":

    # Construct a molecule from scratch
    ints_eval('test/h2.xyz', 'sto-3g')
    
    ort_overlap(olp_matrix)
	
    # Occupation matrix
    occ = np.array([1, 0])
	
    # Nuclear Repulsion Energy
    E_nuc = 0.713176830593
