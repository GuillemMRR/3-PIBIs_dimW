# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 19:02:03 2023

@author: user-pc
"""

import numpy as np
import scipy as sp
import scipy.optimize as sp_opt
from ncpol2sdpa import *

def gm(d):
    def gellmann(j, k, d):  
        if j > k:
            gjkd = np.zeros((d, d), dtype=np.complex128)
            gjkd[j - 1][k - 1] = 1
            gjkd[k - 1][j - 1] = 1
        elif k > j:
            gjkd = np.zeros((d, d), dtype=np.complex128)
            gjkd[j - 1][k - 1] = -1.j
            gjkd[k - 1][j - 1] = 1.j
        elif j == k and j < d:
            gjkd = np.sqrt(2/(j*(j + 1)))*np.diag([1 + 0.j if n <= j
                                                   else (-j + 0.j if n == (j + 1)
                                                         else 0 + 0.j)
                                                   for n in range(1, d + 1)])
        else:
            gjkd = np.diag([1 + 0.j for n in range(1, d + 1)])
        return gjkd
    gm = [gellmann(1+i, 1+j, d) for i in range(d) for j in range(d)]
    return np.array([gm[i] for i in range(d**2)])


#%% Bound qubits POVM -- variational

basis = gm(3*2)
def b1_qubits(gen):
    
    V0 = sp.linalg.expm(1j*np.einsum("i,ijk->jk",gen[0:(3*2)**2], basis))[0:2]
    V1 = sp.linalg.expm(1j*np.einsum("i,ijk->jk",gen[(3*2)**2:2*(3*2)**2], basis))[0:2]
    
    pi_00 = V0 @ np.diag([1,1,0,0,0,0]) @ V0.conjugate().transpose()  # From Naimark dillation
    pi_10 = V0 @ np.diag([0,0,1,1,0,0]) @ V0.conjugate().transpose()
    pi_01 = V1 @ np.diag([1,1,0,0,0,0]) @ V1.conjugate().transpose()
    pi_11 = V1 @ np.diag([0,0,1,1,0,0]) @ V1.conjugate().transpose()
    
    B1 = pi_00 + pi_11 + pi_01 + pi_10 - (pi_00 - pi_11) @ (pi_00 - pi_11) - (pi_01 - pi_10) @ (pi_01 - pi_10)
    return np.linalg.eigvalsh(B1)[0]


xopt, fopt, itera, funcalls, warnflag = sp_opt.optimize.fmin(b1_qubits, np.random.randn(2*(3*2)**2), full_output=True) # Nelder-Mead
print("qubits_variational = ", fopt)  

# We obtain -0.25, altough sometimes it gets stuck in a local minimum. 

#%% Quantum bound 

pi = generate_operators('pi', 4, hermitian=True) # 00,11,10,01

B1 = pi[0] + pi[1] - (pi[0] - pi[1])**2 + pi[2] + pi[3] - (pi[2] - pi[3])**2
# B1 = pi[0] * pi[1] + pi[1] * pi[0] + pi[2] * pi[3] + pi[3] * pi[2] # Proj

# 3-outcome POVM
inequalities_nc = [1 - pi[0] - pi[2]]
inequalities_nc += [1 - pi[1] - pi[3]]
inequalities_nc += [pi[0]] + [pi[1]] + [pi[2]] + [pi[3]] 
substitutions = {pi[0]**2: pi[0]} | {pi[1]**2: pi[1]} | {pi[2]**2: pi[2]} | {pi[3]**2: pi[3]}  # Proj


sdpRelaxation = SdpRelaxation(pi, verbose = 1)
sdpRelaxation.get_relaxation(level=2, objective=B1, inequalities=inequalities_nc)
sdpRelaxation.solve("mosek")
print("SDP status = ", sdpRelaxation.status)
print("bound = ", sdpRelaxation.primal)  

# One can add the projective equality constraints (and remove the other inequalities) solve level = 1 for more efficiency.

