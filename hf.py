import numpy as np
from horton import *


# Construct a molecule from scratch
mol = IOData.from_file(context.get_fn('test/h2.xyz'))

# # Create a Gaussian basis set
obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')

# # Create a linalg factory
lf = DenseLinalgFactory(obasis.nbasis)

# Basiset info
#print "contraction coefficients"

# # Compute Gaussian integrals
olp = obasis.compute_overlap(lf)
olp_m = olp._array
kin = obasis.compute_kinetic(lf)
kin_m = kin._array
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
na_m = na._array
er = obasis.compute_electron_repulsion(lf)
er_m = er._array

print "kin_m", kin_m
print
print "na_m", na_m
print
print "er_m", er_m

X, U = np.linalg.eigh(olp_m)
#print "Diagonalized matrix", X
#print "Unitary matrix", U
#XInv = np.linalg.inv(X)
#Xinsqrt = np.power(Xinv, 0.5)
Xinsqrt2 = np.power(X, -0.5)
Xfin = Xinsqrt2*np.eye(2)

#print "Xinsqrt2", Xinsqrt2
#print Xfin
#print U.T
#print

Ssqt = np.dot(U, np.dot(Xfin, U.T))

# Fock initial guess
H = kin_m + na_m
Fzero = np.dot(Ssqt.T, np.dot(H, Ssqt))

#print "core Hamiltonian", H
#print "Initial Fock", Fzero
#print

# Coeff initial guess
epsilon, Cp = np.linalg.eigh(Fzero)

Czero = np.dot(Ssqt, Cp)
#print "Cp", Cp
#print
#print "Czero nos", Czero
#print
#print "epsilon", epsilon

occ = np.array([1, 0])

# Initial density matrix
D = np.zeros((2,2))
for m in range(len(Czero)):
    for v in range(len(Czero)):
        for i in range(len(Czero)):
            D[m,v] += occ[i]*Czero[m,i]*Czero[v, i]

#D2 = np.dot(occ*Czero, Czero.T)

#print "Density matrix", D
#print
#print "Density matrix", D2

# New Fock matrix
F = np.zeros((2,2))

for m in range(2):
    for v in range(2):
        F[m, v] = H[m, v]
        for p in range(2):
            for s in range(2):
                F[m, v] += D[p, s] * (2 * er_m[m,p,v,s] - er_m[m,v,p,s])


#print "Fock", F

# Energy
Energy_0 = 0

for m in range(2):
    for v in range(2):
        Energy_0 += D[m,v]*(H[m,v] + F[m,v])

enuc = 0.713176830593
Etot = Energy_0 + enuc
#print "Energy", Energy_0 + enuc


Etemp = 0

while abs(Etot - Etemp) > 10e-6:
    Etot = Etemp
    # Transform Fock matrix
    F = np.dot(Ssqt.T, np.dot(F, Ssqt))

    # diagonalize Fock matrix
    epsilon, Cp = np.linalg.eigh(F)

    # construct new eigenvector matrix
    C = np.dot(Ssqt, Cp)

    # form new density matrix
    for m in range(len(C)):
        for v in range(len(C)):
            for i in range(len(C)):
                D[m,v] += occ[i]*C[m,i]*C[v, i]
    
    # new Fock matrix
    for m in range(2):
        for v in range(2):
            F[m, v] = H[m, v]
            for p in range(2):
                for s in range(2):
                    F[m, v] += D[p, s] * (2 * er_m[m,p,v,s] - er_m[m,v,p,s])
    # New energy
    for m in range(2):
        for v in range(2):
            Etemp += D[m,v]*(H[m,v] + F[m,v])

#print "HF Energy", Etemp
print "END"
