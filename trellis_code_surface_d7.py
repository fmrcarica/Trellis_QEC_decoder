# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 14:29:19 2024

@author: filip
"""

"""
This code generates a syndrome trellis for a quantum code. It takes stabilizer 
generators and logical operators as vectors in a symplectic vector space
 of a code and constructs the given 0-syndrome trellis. 

Old code: The trellis is constructed by iterating over all n-Pauli operators in S⊥. For
each given n-qubit Pauli it finds the vertex-edge-vertex pairs of the trellis,
iterating over all qubits. The vertices are intermediate syndromes, here
encoded as integers ranging from 0 to 2^(n_stabilizers)-1. The edges are
single qubit Pauli operators, encoded with two bits [i, j] as X^i * Z^j

New code: The trellis is constructed of Wolf "Trellis-Based Decoding of
Binary Linear Block Codes" (1978). For each generator of S⊥, so 
[stabilizers, logicals], construct the binary trellis, i.e. trellis with only
two paths for [identity, Generator from S⊥]. Now, construct the full trellis
as products of all binary trellises.

The resulting trellis is saved as .txt file as a list with length number of
qubits. Each element is a list of each valid vertex-edge-vertex pair, given as
a list of length 4 
[int(previous_syndrome), int(X), int(Z), int(next_syndrome)]

The logical decoder for this project is written in julia.
"""

import numpy as np
from itertools import product
import itertools
import time
import os
import copy
from tqdm import tqdm


def symplectic_inner_product(u, v):
    n = len(u) // 2
    product = 0
    for i in range(n):
        product += u[i] * v[n + i] + u[n + i] * v[i]
    return product % 2

def is_orthogonal(v,vec_array):
    return all(symplectic_inner_product(v, s) % 2 == 0 for s in vec_array)

# Remove potential duplicates
def remove_duplicates(array_list):
    seen = set()
    unique_arrays = []
    for array in array_list:
        # Convert the array to a tuple so it can be hashed
        array_tuple = tuple(array)
        if array_tuple not in seen:
            unique_arrays.append(array)
            seen.add(array_tuple)
    return unique_arrays



def binary_vector_to_int(binary_vector):
    # Convert the binary vector (numpy array) to a binary string
    binary_str = ''.join(str(bit) for bit in binary_vector)
    # Convert the binary string to an integer
    integer_value = int(binary_str, 2)
    return integer_value

def generate_overlapping_squares(d, num_squares):
    squares = []
    size = 2  # Size of the filled square (2x2)
    
    for i in range(num_squares):
        array = np.zeros((d, d), dtype=int)
        row = (i // (d - 1))
        col = (i % (d - 1))
        array[row:row+size, col:col+size] = 1
        squares.append(array)
        
    return squares

# Define distance code, has to be uneven for this implementation

mode = "Full" #"Full or ""X"
d = 9
code_type = "color666" # "surface" or "color666"

assert(d % 2 == True)
assert(mode == "X" or mode == "Full")


#surface code


directory_name = code_type

if mode == "Full":
    directory_name = directory_name + "_d" + str(d)
elif mode == "X":
    directory_name = directory_name + "_d" + str(d) + '_X'
        

if not os.path.exists(directory_name):
    os.makedirs(directory_name)
    
if code_type == "surface":
    num_stab = d**2 - 1
    num_qubits = d**2

    num_squares = (d-1)**2
    squares = generate_overlapping_squares(d, num_squares)
    stabilizers_X = []
    stabilizers_Z = []
    
    for k, arr in enumerate(squares):
        if (k // (d-1)) % 2:
            if k % 2:
                stabilizers_X.append(arr.reshape(-1))
            else:
                stabilizers_Z.append(arr.reshape(-1))
                
        else:
            if k % 2:
                stabilizers_Z.append(arr.reshape(-1))
            else:
                stabilizers_X.append(arr.reshape(-1))
            
    # This part of stabilizer generation only works for distance d uneven
    def generate_sliding_window_arrays(d):
        arrays = []
        size = 2  # Size of the 1s block
        
        for i in range(d - size + 1):
            array = np.zeros(d, dtype=int)
            array[i:i+size] = 1
            arrays.append(array)
            
        return arrays
    
    # Example usage for d=5
    arrays = generate_sliding_window_arrays(d)
    
    squares_w2 = []
    
    to_fill_in = np.zeros((d-1,d), dtype = int)
    
    for k, arr in enumerate(arrays):
        if k % 2:
            x_stab = np.vstack((arr, to_fill_in))
            stabilizers_X.append(x_stab.reshape(-1))
            
            z_stab = np.vstack((to_fill_in, arr))
            z_stab = z_stab.T
            stabilizers_Z.append(z_stab.reshape(-1))
        else:
            x_stab = np.vstack((to_fill_in, arr))
            stabilizers_X.append(x_stab.reshape(-1))
            
            z_stab = np.vstack((arr, to_fill_in))
            z_stab = z_stab.T
            stabilizers_Z.append(z_stab.reshape(-1))
        
    
        
    zeros_for_symplectic_order = np.zeros((d**2),dtype=int)
    
    stabilizers_X = [np.vstack((arr,zeros_for_symplectic_order)).reshape(-1) for arr in stabilizers_X]
    stabilizers_Z = [np.vstack((zeros_for_symplectic_order,arr)).reshape(-1) for arr in stabilizers_Z]
    
    logicalx = np.vstack((np.ones(num_qubits, dtype=int),np.zeros(num_qubits, dtype=int))).reshape(-1)
    
    logicalz = np.vstack((np.zeros(num_qubits, dtype=int),np.ones(num_qubits, dtype=int))).reshape(-1)

if code_type == "color666":
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
    num_stab = 0


    for j in range(lattice_N):
        for i in range(j + 1):
            if lattice_sites[j][i] == "d":
                num_qubits += 1
                
                
    for j in range(lattice_N):
        for i in range(j + 1):
            if lattice_sites[j][i] == "s":
                num_stab += 1
            
            
    stabilizers = []
    num_lattice_sites = num_qubits + num_stab
                
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

'''
# Steane code
num_stab = 6
num_qubits = 7

# Given stabilizers
stabilizer1 = np.array([1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
stabilizer2 = np.array([0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
stabilizer3 = np.array([0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])

stabilizer4 = np.array([ 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0])
stabilizer5 = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0])
stabilizer6 = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,])

logicalx = np.array([1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0])

logicalz = np.array([ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1])

'''

# List of all possible binary vectors of dimension 7

stabilizers = stabilizers_X + stabilizers_Z
logicals = [logicalx, logicalz]

destabilizers = []

def int_to_binary_vector(n, length=2*num_qubits):
    # Convert integer to binary string, remove '0b' prefix, and pad with zeros up to the specified length
    binary_str = bin(n)[2:].zfill(length)
    # Convert binary string to a numpy array of integers (0 and 1)
    binary_vector = np.array([int(bit) for bit in binary_str])
    return binary_vector

#orthogonal_vectos = [np.array(v) for v in binary_vectors if is_orthogonal(v, stabilizers)]

def generate_unit_vector(n_qubits, i):
    # Create a zero vector of length 2 * n_qubits
    vector = np.zeros(2 * n_qubits, dtype=int)
    
    # Set the i-th position to 1
    vector[i] = 1
    
    return vector

# Calculate syndrome as binary vector and convert to integer
def calculate_syndrome(stabilizers, vec):
    syndrome = 0
    for k, stabilizer in enumerate(stabilizers):
        s = symplectic_inner_product(vec, stabilizer)
        syndrome = syndrome + 2**k * s
    return syndrome


print("Constructing destabilizer basis")
# Find destabilizers for given set of stabilizer generators and logicals

for stabilizer in tqdm(stabilizers):
    for k in range(2*num_qubits):
        if symplectic_inner_product(generate_unit_vector(num_qubits, k), stabilizer):
            destabilizers.append(generate_unit_vector(num_qubits, k))
            break
        
stabilizer_tableau = stabilizers + [logicalx] + destabilizers + [logicalz]
real_destabilizers = []

for l, destabilizer in tqdm(enumerate(destabilizers)):
    for k in range(l+ 1, num_qubits):
        
        
        if symplectic_inner_product(stabilizer_tableau[k], destabilizer):
            destabilizers[l] = (destabilizers[l] + stabilizer_tableau[num_qubits + k]) % 2
            
        if symplectic_inner_product(stabilizer_tableau[k + num_qubits], destabilizer):
            destabilizers[l] = (destabilizers[l] + stabilizer_tableau[k]) % 2
            
for l, destabilizer in tqdm(enumerate(destabilizers)):
    for k in range(0, l):
        
        if symplectic_inner_product(stabilizer_tableau[k], destabilizer):
            destabilizers[l] = (destabilizers[l] + destabilizers[k]) % 2
            
for l, destabilizer in tqdm(enumerate(destabilizers)):
    for k in range(num_stab):
        
        if symplectic_inner_product(destabilizers[k], destabilizer) and k != l:
            destabilizers[l] = (destabilizers[l] + stabilizers[k]) % 2
            
for l, destabilizer in tqdm(enumerate(destabilizers)):
    for k in range(num_stab):
        
        if symplectic_inner_product(destabilizers[k], destabilizer) and k != l:
            destabilizers[l] = (destabilizers[l] + stabilizers[k]) % 2
            
for l, destabilizer in tqdm(enumerate(destabilizers)):
    for k in range(len(logicals)):
        
        if symplectic_inner_product(logicals[k], destabilizer):
            destabilizers[l] = (destabilizers[l] + logicals[(k+1)%2]) % 2
            
for l, destabilizer in tqdm(enumerate(destabilizers)):
    for k in range(num_stab):
        
        if symplectic_inner_product(destabilizers[k], destabilizer) and k != l:
            destabilizers[l] = (destabilizers[l] + stabilizers[k]) % 2
            
        
list_of_lists = [arr.tolist() for arr in stabilizers]

# Convert the list of lists to a string
list_str = str(list_of_lists)

name_txt = directory_name + '/stabilizers' + '.txt'

# Save the string to a .txt file
with open(name_txt, 'w') as f:
    f.write(list_str)

list_of_lists = [arr.tolist() for arr in destabilizers]

# Convert the list of lists to a string
list_str = str(list_of_lists)

name_txt = directory_name + '/destabilizers' + '.txt'

# Save the string to a .txt file
with open(name_txt, 'w') as f:
    f.write(list_str)

list_of_lists = [arr.tolist() for arr in logicals]

# Convert the list of lists to a string
list_str = str(list_of_lists)

name_txt = directory_name + '/logicals' + '.txt'

# Save the string to a .txt file
with open(name_txt, 'w') as f:
    f.write(list_str)
 

# Trellis list with vertex-edge-vertex pairs for each qubit
trellis = [[] for _ in range(num_qubits)]

# Fill in trellis
# Iterate over n-qubit Paulis
if mode == "Full":
    combined_list = stabilizers + logicals
elif mode == "X":
    combined_list = stabilizers_Z + [logicalz]
    
#coefficients = [0, 1]

#iterations = [ ]

binary_trellises = []


print("Constructing binary trellises")
for Pauli in tqdm(combined_list):
    
    binary_trellis = [[] for _ in range(num_qubits)]
    
    for vec in [Pauli, np.zeros(num_qubits*2)]:
        vec_it = np.zeros(2*num_qubits)
        s = calculate_syndrome(stabilizers, vec_it)
        for k in range(num_qubits):
            
            # Define previous vertex value (past syndrome)
            s_prev = s
            
            # Define edge (Single qubit Pauli)
            vec_it[k] = vec[k]
            vec_it[k+num_qubits] = vec[k+num_qubits]
            
            # Define next vertex value (future syndrome)
            s = calculate_syndrome(stabilizers, vec_it)
            
            # Write into trellis
            trellis_state = np.array([int(s_prev),int(vec_it[k]),int(vec_it[k+num_qubits]),int(s)])
            
            binary_trellis[k].append(trellis_state)
            
    for k in range(len(binary_trellis)):
        binary_trellis[k] = remove_duplicates(binary_trellis[k])
        
    binary_trellises.append(binary_trellis)
        

# Binary trellis_products
trellis = copy.deepcopy(binary_trellises[0])


print("Constructing full trellis from binary trellises")
for binary_trellis in tqdm(binary_trellises[0:]):
    
    
    for k in range(num_qubits):
        for bin_trel_state in binary_trellis[k]:
            for l in range(len(trellis[k])):
                trellis[k].append(np.array([int(bin_trel_state[0]) ^ int(trellis[k][l][0]), 
                                 int(trellis[k][l][1] + bin_trel_state[1]) % 2, 
                                 int(trellis[k][l][2] + bin_trel_state[2]) % 2 , 
                                 int(bin_trel_state[3]) ^ int(trellis[k][l][3])]))
                
                
        trellis[k] = remove_duplicates(trellis[k])
    

'''
print("Finding S⊥")
for coeffs in tqdm(itertools.product(coefficients, repeat=len(combined_list[:14]))):
    start_time = time.time()
    # Linear combination using the current set of coefficients
    vec = (sum(c * arr for c, arr in zip(coeffs, combined_list)))
    
    vec = vec % 2
    
    iterations.append(binary_vector_to_int(vec))
    

print("Constructing trellis")
for it in tqdm(iterations):
    start_time = time.time()
    # Convert into binary vector
    vec = int_to_binary_vector(it)
    #print(vec)
    # Check if Pauli is in S⊥
    if is_orthogonal(vec, stabilizers):
        
        # Define intermediate Pauli for iteration over all qubits
        vec_it = np.zeros(2*num_qubits)
        s = calculate_syndrome(stabilizers_X, vec_it)
        for k in range(num_qubits):
            #print(s)
            
            # Define previous vertex value (past syndrome)
            s_prev = s
            
            # Define edge (Single qubit Pauli)
            vec_it[k] = vec[k]
            vec_it[k+num_qubits] = vec[k+num_qubits]
            
            # Define next vertex value (future syndrome)
            s = calculate_syndrome(stabilizers_X, vec_it)
            
            # Write into trellis
            trellis_state = np.array([s_prev,vec_it[k],vec_it[k+num_qubits],s]).astype(int)
            
            trellis[k].append(trellis_state)
            
    else:
        print("Error occured in trellis construction")
        break
    # Record the end time
    end_time = time.time()
    
    # Calculate the runtime
    runtime = end_time - start_time
    #print(runtime)
    
    '''
    
'''
    start_time = time.time()
    #print(vec)
    vec_it = np.zeros(2*num_qubits)
    s = calculate_syndrome(stabilizers, vec_it)
    for k in range(num_qubits):
        
        # Define previous vertex value (past syndrome)
        s_prev = s
        
        # Define edge (Single qubit Pauli)
        vec_it[k] = vec[k]
        vec_it[k+num_qubits] = vec[k+num_qubits]
        
        # Define next vertex value (future syndrome)
        s = calculate_syndrome(reversed(stabilizers), vec_it)
        
        #print(np.array([s_prev,vec_it[k],vec_it[k+num_qubits],s]))
        
        # Write into trellis
        trellis_state = np.array([s_prev,vec_it[k],vec_it[k+num_qubits],s]).astype(int)
        
        trellis[k].append(trellis_state)
    # Record the end time
    end_time = time.time()
    
    # Calculate the runtime
    runtime = end_time - start_time
    print(runtime)
    # Do something with the linear_combination
    #print(f"Coefficients: {coeffs}, Linear Combination: {linear_combination}")
'''
'''
print("Constructing trellis network")
for it in tqdm(range(2**(num_qubits))):
    start_time = time.time()
    # Convert into binary vector
    vec = int_to_binary_vector(it)
    #print(vec)
    # Check if Pauli is in S⊥
    if is_orthogonal(vec, stabilizers):
        
        # Define intermediate Pauli for iteration over all qubits
        vec_it = np.zeros(2*num_qubits)
        s = calculate_syndrome(stabilizers, vec_it)
        for k in range(num_qubits):
            #print(s)
            
            # Define previous vertex value (past syndrome)
            s_prev = s
            
            # Define edge (Single qubit Pauli)
            vec_it[k] = vec[k]
            vec_it[k+num_qubits] = vec[k+num_qubits]
            
            # Define next vertex value (future syndrome)
            s = calculate_syndrome(stabilizers, vec_it)
            
            # Write into trellis
            trellis_state = np.array([s_prev,vec_it[k],vec_it[k+num_qubits],s]).astype(int)
            
            trellis[k].append(trellis_state)
    # Record the end time
    end_time = time.time()
    
    # Calculate the runtime
    runtime = end_time - start_time
    print(runtime)
'''


for k in range(len(trellis)):
    trellis[k] = remove_duplicates(trellis[k])
    
zero_trellis =  copy.deepcopy(trellis)


def array_to_string(arr):
    return '[' + ', '.join(map(str, arr.tolist())) + ']'

# Convert each numpy array to string and organize them
formatted_data = [[array_to_string(arr) for arr in sublist] for sublist in trellis]

# Join the arrays with proper formatting
formatted_string = '[' + '], ['.join([', '.join(sublist) for sublist in formatted_data]) + ']'

formatted_string = '[' + formatted_string + ']'

name_txt = directory_name + '/trellis' + '.txt'

# Write to file
with open(name_txt, 'w') as f:
    f.write(formatted_string)
    
