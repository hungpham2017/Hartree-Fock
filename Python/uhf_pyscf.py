'''
Unrestricted Hartree-Fock using PySCF intergrals
Author: Hung Q. Pham (UMN, phamx494@umn.edu)
'''

import numpy as np
import scipy as sp
from pyscf import gto, scf, ao2mo
from functools import reduce

mol = gto.M()
mol.atom =("""
 O                 -1.82412712    0.02352916    0.00000000
 O                 -0.72539390    0.07630510    0.00000000
 N                 -2.14114302    1.07685766    0.00000000
""")

mol.spin = 1
mol.basis = '6-31g'
mol.build()
Norb = mol.nao_nr()

#UHF convergent criteria
e_conv = 1.e-8
d_conv = 1.e-8
if (mol.nelectron %2 != 0):
	nel_a = int(0.5 * (mol.nelectron + 1))
	nel_b = int(0.5 * (mol.nelectron - 1))
else:
	nel_a = int(0.5 * mol.nelectron)
	nel_b = nel_a

damp_value = 0.20
damp_start = 1
ncycl = 100

# Core Hamiltonian
V = mol.intor_symmetric('cint1e_nuc_sph')
T = mol.intor_symmetric('cint1e_kin_sph')
H = T + V

S = mol.intor_symmetric('cint1e_ovlp_sph') 

#Get two-e integral (TEI) from PySCF format and convert it to chemist format 
g_format = mol.intor('cint2e_sph')
g = np.zeros([Norb, Norb, Norb, Norb]) 
for cn1 in range(Norb):
	for cn2 in range(Norb):
		for cn3 in range(Norb):
			for cn4 in range(Norb):
				g[cn1][cn2][cn3][cn4] = g_format[cn1*Norb+cn2][cn3*Norb+cn4]


# print(S.shape)
# print(I.shape)

A = S
A = sp.linalg.fractional_matrix_power(A, -0.5)

# print(A @ S @ A)


# Diagonalize Core H
def diag(F, A):
	Fp = reduce(np.dot, (A.T, F, A))
	eps, Cp = np.linalg.eigh(Fp)
	C = A.dot(Cp)
	return eps, C


eps_a, C_a = diag(H, A)
eps_b, C_b = diag(H, A)
Cocc_a, Cocc_b = C_a[:, :nel_a], C_b[:, :nel_b]
D_a, D_b = Cocc_a.dot(Cocc_a.T), Cocc_b.dot(Cocc_b.T)
D = D_a + D_b


#SCF procedure
E_old = 0.0
Fa_old = None
Fb_old = None
maxcyc = None
for iteration in range(ncycl):

    J = np.einsum("pqrs,rs->pq", g, D)
    Ka = np.einsum("prqs,rs->pq", g, D_a)
    Kb = np.einsum("prqs,rs->pq", g, D_b)
	
	
    Fa_new = H + J - Ka
    Fb_new = H + J - Kb
	
    # conditional iteration > start_damp
    if iteration >= damp_start:
        Fa = damp_value * Fa_old + (1.0 - damp_value) * Fa_new
        Fb = damp_value * Fb_old + (1.0 - damp_value) * Fb_new		
    else:
        Fa = Fa_new
        Fb = Fb_new
		
    Fa_old = Fa_new
    Fb_old = Fb_new	
	

    # Build the AO gradient
    grada = reduce(np.dot,(Fa, D_a, S)) - reduce(np.dot, (S, D_a, Fa))
    grada_rms = np.mean(grada ** 2) ** 0.5

    gradb = reduce(np.dot,(Fb, D_b, S)) - reduce(np.dot, (S, D_b, Fb))
    gradb_rms = np.mean(gradb ** 2) ** 0.5
    grad_rms = 	grada_rms + gradb_rms
	
    # Build the energy
    E_electric = 0.5 * np.sum( H * D + D_a * Fa + D_b * Fb)
    E_total = E_electric + mol.energy_nuc()

    E_diff = E_total - E_old
    E_old = E_total
    print("Iter=%3d  E = % 16.12f  E_diff = % 8.4e  D_diff = % 8.4e" %
            (iteration, E_total, E_diff, grad_rms))

    # Break if e_conv and d_conv are met
    if (abs(E_diff) < e_conv) and (grad_rms < d_conv):
        maxcyc = iteration
        break
		
    if (iteration == ncycl - 1):
        maxcyc = ncycl

    eps_a, C_a = diag(Fa, A)
    eps_b, C_b = diag(Fb, A)
    Cocc_a, Cocc_b = C_a[:, :nel_a], C_b[:, :nel_b]
    D_a, D_b = Cocc_a.dot(Cocc_a.T), Cocc_b.dot(Cocc_b.T)
    D = D_a + D_b

if (maxcyc == ncycl):
	print("SCF has not finished!\n")
else:
	print("SCF has finished!\n")
	
#PYSCF solution
'''
Note: PySCF uses a combination of many SCF techniques, this include DIIS, generating initial guess 1-rdm 
based on ANO basis'minao',...
For RHF solution, usually there is no difference when using/not using these techniques.
However, this is not the case for UHF solution. To compare with simple algorithm UHF above, I used damping technque
for both codes. This is not the best technique for UHF convergences

Will updated the DIIS in this SCF
'''


mf = scf.UHF(mol)
mf.init_guess = '1e'
mf.max_cycle = 200
mf.diis_start_cycle=200
mf.diis = False
mf.damp = damp_value #using damping instead of DIIS
mf.kernel()


pyscf_energy = mf.e_tot
print("Energy matches PySCF %s" % np.allclose(pyscf_energy, E_total))