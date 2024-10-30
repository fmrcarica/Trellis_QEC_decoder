# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 16:31:14 2024

@author: filip
"""

import numpy as np
from itertools import product
import itertools
import time
import os
import copy
from tqdm import tqdm

d = 7

assert(d % 2 == 1)

# Generate 6.6.6 planar lattice
lattice_N = d + int((d-1)/2)
lattice_sites = [np.empty(i+1, dtype=str) for i in range(lattice_N)]

for j in range(lattice_N):
    for i in range(j + 1):
        lattice_sites[j][i] = "d"
        if j == 1 and i == 0:
            lattice_sites[j][i] = "s"
        #if j == 0 and i == 0:
        #    lattice_sites[j][i] = "d"
        #if i > 0 and (lattice_sites[j][i-1] == "s"):
        #    lattice_sites[j][i] == "d"
        #print(j, i)
        #print(i > 1 and (lattice_sites[j][i-1] == "d" and lattice_sites[j][i-2] == "d"))
        if i > 1 and (lattice_sites[j][i-1] == "d" and lattice_sites[j][i-2] == "d"):
            lattice_sites[j][i] = "s"
        if (i == 0 and j > 1) and (lattice_sites[j-1][i] == "d" and lattice_sites[j-2][i] == "d"):
            lattice_sites[j][i] = "s"
        if i == 1 and (lattice_sites[j][i-1] == "d" and lattice_sites[j-1][i-1] == "d"):
            lattice_sites[j][i] = "s"
            
num_qubits = 0
num_stabs = 0


for j in range(lattice_N):
    for i in range(j + 1):
        if lattice_sites[j][i] == "d":
            num_qubits += 1
            
            
for j in range(lattice_N):
    for i in range(j + 1):
        if lattice_sites[j][i] == "s":
            num_stabs += 1
        
        
stabilizers = []
num_lattice_sites = num_qubits + num_stabs
            
for j in range(lattice_N):
    for i in range(j + 1):
        if lattice_sites[j][i] == "s":
            stabilizer = np.zeros(num_lattice_sites, dtype = int)
            # Test if on the boundary
            if i == 0 or i == j or j == (lattice_N - 1):
                if i == 0:
                    # Fill in stabilizer sites
                    # upper row
                    n = sum(np.arange(j)) + i - 1
                    stabilizer[n + 1] = 1
                    # middle row
                    n = sum(np.arange(j+1)) + i - 1
                    stabilizer[n + 2] = 1
                    # bottom row
                    n = sum(np.arange(j + 2)) + i
                    stabilizer[n] = 1
                    stabilizer[n + 1] = 1
                if i == j:
                    # Fill in stabilizer sites
                    # upper row
                    n = sum(np.arange(j)) + i - 1
                    stabilizer[n] = 1
                    # middle row
                    n = sum(np.arange(j + 1)) + i - 1
                    stabilizer[n] = 1
                    # bottom row
                    n = sum(np.arange(j + 2)) + i
                    stabilizer[n] = 1
                    stabilizer[n + 1] = 1
                if j == (lattice_N - 1):
                    # Fill in stabilizer sites
                    # upper row
                    n = sum(np.arange(j)) + i - 1
                    stabilizer[n] = 1
                    stabilizer[n + 1] = 1
                    # middle row
                    n = sum(np.arange(j + 1)) + i - 1
                    stabilizer[n] = 1
                    stabilizer[n + 2] = 1
                    
                
            else:
                # Fill in stabilizer sites
                # upper row
                n = sum(np.arange(j)) + i - 1
                stabilizer[n] = 1
                stabilizer[n + 1] = 1
                # middle row
                n = sum(np.arange(j + 1)) + i - 1
                stabilizer[n] = 1
                stabilizer[n + 2] = 1
                # bottom row
                n = sum(np.arange(j + 2)) + i
                stabilizer[n] = 1
                stabilizer[n + 1] = 1
                
            stabilizers.append(stabilizer)
            
stab_indices = []
for j in range(num_lattice_sites):
    if all(stab[j] == 0 for stab in stabilizers):
        stab_indices.append(j)
        
for j in reversed(stab_indices):
    stabilizers = [np.delete(stab, j) for stab in stabilizers]
    
# Generate logicals
logicalx = np.zeros(num_lattice_sites, dtype=int)

for j in range(lattice_N):
    n = sum(np.arange(j + 1))
    logicalx[n] = 1
    
for j in reversed(stab_indices):
    logicalx = np.delete(logicalx, j)
    
stabilizers_X = [np.append(stab, np.zeros(num_qubits, dtype= int)) for stab in stabilizers]
stabilizers_Z = [np.append(np.zeros(num_qubits, dtype= int), stab) for stab in stabilizers]

logicalz = np.append(np.zeros(num_qubits, dtype= int), logicalx)
logicalx = np.append(np.zeros(num_qubits, dtype= int), logicalx)

stabilizers = stabilizers_X + stabilizers_Z