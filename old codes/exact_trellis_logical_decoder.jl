#=

Logical decoder for rotated surface - quantum codes. Using representation for tensor networks coming 
representation of the quantum error correcting code as a zero-syndrome trellis, following the implementation in
(doi.org/10.48550/arXiv.2106.08251). This code is to be used accompanied with the trellis constructing code in
python.

The zero-syndrome trellis is imported as a .txt file in format 
[Qubit 1 trellis: [Vertex-Edge-Vertex pairs: [int(past_syndrome), bool(X), bool(Z), int(fut_syndrome)], ...], Qubit 2 trellis, ...]
Given the full stabilizer tableu (stabilizers, destabilizers, logicals), the code constructs the tensor network in
an ITensor implementation, the logicals are implemented in the detector picture, using the Hadamard trick to avoid 
high-weight parity.

=#
using ITensors
using Random
using LinearAlgebra
using AbstractAlgebra
using ProgressBars

mode = "Full" #"Full or ""X"
d = 3 # Only works for d unveven

@assert d % 2 == true 
@assert mode == "X" || mode == "Full" 


# Define number of qubits and stabilizers
n_stabs = Int(d^2 - 1)
if mode == "X"
    global n_stabs = Int(n_stabs / 2)
end
n_qubits = d^2

if mode == "Full"
    directory_name = "surface_d" * string(d)
end
if mode == "X"
    directory_name ="surface_d" * string(d) * "_X"
end

# Open and read in the trellis
filename = directory_name * "/trellis.txt"
file_content = open(filename) do file
    read(file, String)
end
function safe_eval(content)
    # Replace large integers (example for very large numbers)
    content = replace(content, r"(\d+)" => s -> try
        n = parse(Int64, s)
        return s  # Return the number as-is if it fits in Int64
    catch e
        return "BigInt($s)"  # Convert to BigInt if it overflows Int64
    end)
    
    return eval(Meta.parse(content))
end

#zero_syndrome_trellis = safe_eval(file_content)
zero_syndrome_trellis = eval(Meta.parse(file_content))

# Open and read in the stabilizers
filename = directory_name * "/stabilizers.txt"
file_content = open(filename) do file
    read(file, String)
end
stabilizers = eval(Meta.parse(file_content))

stabilizers_z = [ ]

if mode == "X"
    global stabilizers_z = stabilizers[(n_stabs+1):(2*n_stabs)]
    #global stabilizers = stabilizers[1:n_stabs]
end



# Open and read in the logicals
filename = directory_name * "/logicals.txt"
file_content = open(filename) do file
    read(file, String)
end
logical_operators = eval(Meta.parse(file_content))

logical_operator_x = logical_operators[1]
logical_operator_z = logical_operators[2]

# Open and read in the destabilizers
filename = directory_name * "/destabilizers.txt"
file_content = open(filename) do file
    read(file, String)
end
destabilizers = eval(Meta.parse(file_content))

if mode == "X"
    global destabilizers = destabilizers[1:n_stabs]
end

function symplectic_product(v1, v2, n)
    # The symplectic product of two binary vectors v1 and v2
    return sum(v1[i] * v2[n+i] + v1[n+i] * v2[i] for i in 1:n) % 2
end


logical_weigth = count(x -> x == 1, logical_operator_x)


# Define the future and past intermediate syndromes
syndrome_values_past = [[vec[1] for vec in arr] for arr in zero_syndrome_trellis]
syndrome_values_future = [[vec[4] for vec in arr] for arr in zero_syndrome_trellis]


# cut out duplicates
syndrome_values_future = [unique(arr) for arr in syndrome_values_future]
syndrome_values_past = [unique(arr) for arr in syndrome_values_past]


# First index
s0 = Index(1) # 000000-syndrome
# intermediate syndromes
s = [Index(length(syndrome_values_future[i])) for i in 1:n_qubits]
s = pushfirst!(s, s0)

# Helper function to get syndrome
function check_orthogonality(vectors::Vector{Vector{Int64}}, vector::Vector{Int64})
    n = size(vectors, 1)
    orthogonality = Vector{Bool}(undef, n)
    
    for i in 1:n
        orthogonality[i] = symplectic_product(vectors[i], vector,n_qubits)
    end
    
    return orthogonality
end


function Pauli_error(p::Float64, n::Int64)
    Pauli_error = Vector{Int}(undef, 2*n)
    p_Pauli = p
    if mode == "X"
        for i in 1:n
            Pauli_error[i], Pauli_error[i+n] = 0, 0
            r = rand()
            if r < p_Pauli
                Pauli_error[i + n_qubits] = 1
            end
        end
    end
    
    if mode == "Full"
        for i in 1:n
            Pauli_error[i], Pauli_error[i+n] = 0, 0
            r = rand()
            if r < p_Pauli
                Pauli_error[i] = 1
            else
                if r <= 2*p_Pauli
                    Pauli_error[i+n] = 1
                else
                    if r <= 3*p_Pauli
                        Pauli_error[i] = 1
                        Pauli_error[i+n] = 1
                    end
                end
            end
        end
    end
    
    return Pauli_error
end

p_qubit = 0.1

p = p_qubit/3

if mode == "X"
    global p = p_qubit
end


occured_error = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
#Pauli_error(p, n_qubits)

s_meas = check_orthogonality(stabilizers, occured_error)


function find_Pauli_give_syndrome(stabilizers, s_measured,n)
    poss = 2^(2*n)-1
    v = 0
    padded_vector = 0
    for i in 0:poss
        v = digits(i, base = 2)
        padded_vector = fill(0, 2*n)
        padded_vector[1:length(v)] = v
        if check_orthogonality(stabilizers, padded_vector) == s_measured
            break
        end
    end
    return padded_vector
end

#PS = find_Pauli_give_syndrome(stabilizers, s_meas, n_qubits)
PS = zeros(Int, 2*n_qubits)

for i in 1:n_stabs
    if s_meas[i]
        global PS = (PS .+ destabilizers[i]) .% 2
    end
end

x = [Index(2) for i in 1:n_qubits]
z = [Index(2) for i in 1:n_qubits]

f = [Array{Int}(undef, length(syndrome_values_past[k]), 2, 2, length(syndrome_values_future[k])) for k in 1:length(syndrome_values_future)]
#push!(f,Array{Int}(undef, length(syndrome_values_past[7]), 2, 2, 64))





function binary_vector_to_reversed_int(binary_vector)::Int64
    # Reverse the binary vector
    
    # Convert the reversed binary vector to an integer
    result = 0
    for (i, bit) in enumerate(binary_vector)
        result += bit * (2^(i-1))
    end
    
    return result
end

# intialize variables
intermediate_syndrome_past = 0
intermediate_syndrome_future = 0

# Fill in indicator network from trellis
println("Constructing trellis network")
for trellis_pos in ProgressBar(1:length(zero_syndrome_trellis))
    
    global intermediate_syndrome_past = intermediate_syndrome_future

    intermediate_vector = fill(0, 2*n_qubits)
    intermediate_vector[1:trellis_pos] = PS[1:trellis_pos]
    intermediate_vector[n_qubits + 1:(n_qubits + trellis_pos)] = PS[n_qubits + 1:(n_qubits + trellis_pos)]

    global intermediate_syndrome_future = binary_vector_to_reversed_int(check_orthogonality(stabilizers, intermediate_vector))
    
    for (i, past_syndrome) in enumerate(syndrome_values_past[trellis_pos])
        for j in 1:2
            for k in 1:2
                for (l, future_syndrome) in enumerate(syndrome_values_future[trellis_pos])
                    # Condition: Sum of indices is even
                    if [past_syndrome, (j-1) ⊻ PS[trellis_pos] , (k-1) ⊻ PS[n_qubits + trellis_pos], future_syndrome] in zero_syndrome_trellis[trellis_pos]
                        f[trellis_pos][i, j, k, l] = 1;
                    else
                        f[trellis_pos][i, j, k, l] = 0;
                    end
                end
            end
        end
    end
end

zero_syndrome_trellis = nothing

f_tensor = [ITensor(f[i], s[i], x[i], z[i], s[i+1]) for i in 1:n_qubits] 

beginning_node = Array{Int}(undef, 1)
beginning_node[1] = 1 
f0 = ITensor(beginning_node, s0)

end_node = Array{Int}(undef, 1)
end_node[1] = 1
flast = ITensor(end_node, last(s))

x_P = [Index(2,"P_x_$i") for i in 1:n_qubits]
z_P = [Index(2,"P_z_$i") for i in 1:n_qubits]

#three edge identity
ident_3 = Array{Float64}(undef, 2, 2, 2)
for ind_i in 1:2
    for ind_j in 1:2
        for ind_k in 1:2
            if ind_i == ind_j && ind_i == ind_k
                ident_3[ind_i, ind_j, ind_k] = 1
            else
                ident_3[ind_i, ind_j, ind_k] = 0
            end
        end
    end
end

#two edge identity
ident_2 = Array{Float64}(undef, 2, 2)
for ind_i in 1:2
    for ind_j in 1:2
        if ind_i == ind_j
            ident_2[ind_i, ind_j] = 1
        else
            ident_2[ind_i, ind_j] = 0
        end
    end
end


lx = [ ]
lz = [ ]

for i in 1:n_qubits
    if logical_operator_x[i] == 1
        push!(lx, Index(2, "lx_$i") )
    end
    if logical_operator_z[i + n_qubits] == 1
        push!(lz, Index(2, "lz_$i"))
    end
end

qubit_x = [ ]
qubit_z = [ ]

log_counter_x = 0
log_counter_z = 0
for i in 1:n_qubits
    if logical_operator_x[i] == 1
        global log_counter_x += 1
        push!(qubit_x, ITensor(ident_3,x_P[i], x[i], lx[log_counter_x]))
    else 
        push!(qubit_x, ITensor(ident_2,x_P[i], x[i]))
    end

    if logical_operator_z[i + n_qubits] == 1
        global log_counter_z += 1
        push!(qubit_z, ITensor(ident_3,z_P[i], z[i], lz[log_counter_z]))
    else 
        push!(qubit_z, ITensor(ident_2,z_P[i], z[i]))
    end
end

# Define probability matrix for depolarizing channel
Probability_matrix = Array{Float64}(undef, 2, 2)
Probability_matrix[1,1], Probability_matrix[2, 2] = 1-3*p, p
Probability_matrix[1,2], Probability_matrix[2, 1] = p, p

if mode == "X"
    global Probability_matrix[1,1], Probability_matrix[2, 2] = 1-p, 0
    global Probability_matrix[1,2], Probability_matrix[2, 1] = p, 0
end
P = [ITensor(Probability_matrix,x_P[i], z_P[i]) for i in 1:n_qubits]

Hadamard_matrix = Array{Float64}(undef, 2, 2)
Hadamard_matrix[1,1], Hadamard_matrix[1,2] = 1, 1
Hadamard_matrix[2,1], Hadamard_matrix[2,2] = 1, -1

indices_hadamard_x = [Index(2, "Hadamard_x_$i") for i in 1:logical_weigth]
Hadamard_tensors_x = [ITensor(Hadamard_matrix, lx[i],indices_hadamard_x[i]) for i in 1:logical_weigth]

indices_hadamard_z = [Index(2, "Hadamard_z_$i") for i in 1:logical_weigth]
Hadamard_tensors_z = [ITensor(Hadamard_matrix, lz[i],indices_hadamard_z[i]) for i in 1:logical_weigth]

#logical_class_x = [Index(2, "logical_x_$i") for i in 1:n_qubits]

# Create an IndexSet with these indices
#indices_x = IndexSet(push!(indices_hadamard_x, logical_class_x))
logical_x = [ITensor(indices_hadamard_x[i]) for i in 1:logical_weigth]
logical_z = [ITensor(indices_hadamard_z[i]) for i in 1:logical_weigth]

# Contraction steps
# initialize contracted Tensor
f_contracted = f0


Logical_probabilities = Array{Float64}(undef, 2, 2)

log_counter_x = 0
log_counter_z = 0

println("Contracting TN")
for i in ProgressBar(1:length(Logical_probabilities))
    global f_contracted = f0
    global log_counter_x = 0
    global log_counter_z = 0
    #println(f_contracted)
    for j in ProgressBar(1:n_qubits)
        
        # Update logical counters
        if logical_operator_x[j] == 1
            global log_counter_x += 1
        end
        if logical_operator_z[j + n_qubits] == 1
            global log_counter_z += 1
        end

        if i == 1 
            logical_x[log_counter_x][1], logical_x[log_counter_x][2] = 1, 0
            logical_z[log_counter_z][1], logical_z[log_counter_z][2] = 1, 0
        end
        if i == 2
            logical_x[log_counter_x][2], logical_x[log_counter_x][1] = 1, 0
            logical_z[log_counter_z][1], logical_z[log_counter_z][2] = 1, 0
        end
        if i == 3 
            logical_x[log_counter_x][1], logical_x[log_counter_x][2] = 1, 0
            logical_z[log_counter_z][2], logical_z[log_counter_z][1] = 1, 0
        end
        if i == 4
            logical_x[log_counter_x][2], logical_x[log_counter_x][1] = 1, 0
            logical_z[log_counter_z][2], logical_z[log_counter_z][1] = 1, 0
        end
        
        #=
        logical_x[j][1] = (logical_x[j][1] + PS[j]) % 2
        logical_x[j][2] = (logical_x[j][2] + PS[j]) % 2
        logical_z[j][1] = (logical_z[j][1] + PS[j + n_qubits]) % 2
        logical_z[j][2] = (logical_z[j][2] + PS[j + n_qubits]) % 2
        =#
        if logical_operator_x[j] == 1
            global f_contracted = f_contracted*f_tensor[j]*qubit_x[j]* qubit_z[j]* P[j] * Hadamard_tensors_x[log_counter_x] * logical_x[log_counter_x]
        else
            global f_contracted = f_contracted*f_tensor[j]*qubit_x[j]* qubit_z[j]* P[j]
        end

        if logical_operator_z[j + n_qubits] == 1
            global f_contracted = f_contracted *  Hadamard_tensors_z[log_counter_z] * logical_z[log_counter_z]
        end
            #println(f_contracted)
        #Logical_probabilities[i] = f_contracted[1]

    end
    global f_contracted = f_contracted* flast
    Logical_probabilities[i] = f_contracted[1]
end




println("Physical Error")
println(occured_error)

println("Detected syndrome")
println(s_meas)

# Hadamard transform of data
Logical_probabilities = Hadamard_matrix * Logical_probabilities * Hadamard_matrix

# Logical error of maximum probability
index_max = argmax(Logical_probabilities)

correcting_operation = zeros(Int, n_qubits*2)

println("Most likely error to have occured")
if index_max == CartesianIndex(1, 1)
    println("No logical error")
end
if index_max == CartesianIndex(2, 1)
    println("Logical x error")
end
if index_max == CartesianIndex(1, 2)
    println("Logical z error")
end
if index_max == CartesianIndex(2, 2)
    println("Logical y error")
end

println("Correction operation:")
if index_max == CartesianIndex(1, 1)
    global correcting_operation = PS
    println(PS)
end
if index_max == CartesianIndex(2, 1)
    global correcting_operation = (PS .+ logical_operator_x) .% 2 
    println((PS .+ logical_operator_x) .% 2 )
end
if index_max == CartesianIndex(1, 2)
    global correcting_operation = (PS .+ logical_operator_z) .% 2 
    println((PS .+ logical_operator_z) .% 2 )
end
if index_max == CartesianIndex(2, 2)
    global correcting_operation = (PS .+ logical_operator_z .+ logical_operator_x) .% 2 
    println((PS .+ logical_operator_x .+ logical_operator_z) .% 2 )
end


function int_to_binary_vector(n, length=9)
    # Check if the input number is within the valid range
    if n < 0 || n >= 2^length
        error("The number is out of the valid range 0 to 2^$length - 1")
    end
    
    # Convert integer to binary string, remove "0b" prefix, and pad with zeros up to the specified length
    binary_str = bitstring(n)[end-length+1:end] |> reverse
    
    # Convert binary string to a BitVector
    binary_vector = BitVector([c == '1' for c in binary_str])
    
    return binary_vector
end
#=
global comb = BitVector(fill(0, 2*n_qubits))
# check if correcting operation was correcting
for i in 0:(2^n_stabs - 1)
    vec = int_to_binary_vector(i, n_stabs)
    global comb = BitVector(fill(0, 2*n_qubits))
    for j in 1:n_stabs
        if vec[j] == 1
            global comb = (comb .+ stabilizers[j]) .% 2  # XOR operation for mod 2 addition
        end
    end
    if comb == ((correcting_operation .+ occured_error) .% 2)
        println("Correction worked")
        println(((correcting_operation .+ comb) .% 2))
        break
    end
end
=#

# Function to reduce a matrix over GF(2) to its RREF
function rref_gf2(A::Matrix{Int})::Matrix{Int}
    m, n = size(A)
    A = A .% 2  # Ensure the matrix is in GF(2)
    lead = 1
    for r in 1:m
        if n < lead
            return A
        end
        i = r
        while A[i, lead] == 0
            i += 1
            if i > m
                i = r
                lead += 1
                if n < lead
                    return A
                end
            end
        end

        A[i, :], A[r, :] = A[r, :], A[i, :]  # Swap rows
        A[r, :] .= A[r, :] ./ A[r, lead]  # Scale the leading entry to 1

        for i in 1:m
            if i != r
                A[i, :] .= (A[i, :] - A[i, lead] .* A[r, :]) .% 2  # Subtract to clear the column
            end
        end

        lead += 1
    end
    return A
end

# Function to check if the vector is in the span
function is_in_span(test_op::Vector{Int}, vectors)::Bool
    # Ensure vectors and test_op are in GF(2)
    test_op = test_op .% 2
    vectors = vectors .% 2
    
    # Append the test_op to the vectors as the last row
    augmented_matrix = vcat(vectors, test_op')
    # Row-reduce the augmented matrix
    rref_matrix = rref_gf2(augmented_matrix) .% 2

    # Check the last row in the RREF matrix
    return all(x -> x == 0, rref_matrix[end, :])
end

# Example usage

# Random binary vectors of length 2*n_qubits
vectors = transpose(reduce(hcat, stabilizers))
if mode == "X"
    global vectors = transpose(reduce(hcat, stabilizers_z))
end
test_op = (correcting_operation .+ occured_error) .% 2

# Check if test_op is in the span
in_span = is_in_span(test_op, vectors)
if in_span
    println("Correction successfull")
else
    println("Correction not successfull")
end