#=

Logical decoder for rotated surface - quantum codes. Using representation for tensor networks coming 
representation of the quantum error correcting code as a zero-syndrome trellis, following the implementation in
(doi.org/10.48550/arXiv.2106.08251). This code is to be used accompanied with the trellis constructing code in
python.

The zero-syndrome trellis is imported as a .txt file in format 
[Qubit 1 trellis: [Vertex-Edge-Vertex pairs: [int(past_syndrome), bool(X), bool(Z), int(fut_syndrome)], ...], Qubit 2 trellis, ...]
Given the full stabilizer tableu (stabilizers, destabilizers, logicals), the code constructs the tensor network in
an ITensor implementation, the logicals are implemented in the detector picture, using the Hadamard trick to avoid 
high-weight parity checks.

=#
using ITensors
using Random
using LinearAlgebra
using AbstractAlgebra
using ProgressBars

include("helper_functions.jl")

path = "C:/Users/filip/Documents/ETH/Quantum Engineering MSc/Semester Project - QIT/"


mode = "Full" #"Full or ""X"
d = 7 # Only works for d unveven

cutoff_svd = 1e-5
truncate_dim = 2^20

local_cutoff_svd = 0
local_truncate_dim = 2^15

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

directory_name = path * directory_name

# Open and read in the trellis
filename = directory_name * "/trellis.txt"
file_content = open(filename) do file
    read(file, String)
end
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

<
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
s = [Index(length(syndrome_values_future[i]), "s_$i") for i in 1:n_qubits]
s = pushfirst!(s, s0)


p_qubit = 0.08

p = p_qubit/3

if mode == "X"
    global p = p_qubit
end


occured_error = Pauli_error(p, n_qubits)

s_meas = check_orthogonality(stabilizers, occured_error)



#PS = find_Pauli_give_syndrome(stabilizers, s_meas, n_qubits)
PS = zeros(Int, 2*n_qubits)

#=
for i in 1:n_stabs
    if s_meas[i]
        global PS = (PS .+ destabilizers[i]) .% 2
    end
end
=#

f_x = [Index(2, "trellis_p_x$i") for i in 1:n_qubits]
f_z = [Index(2, "trellis_p_z$i") for i in 1:n_qubits]

x = [Index(2, "x_edge_$i") for i in 1:n_qubits]
z = [Index(2, "z_edge_$i") for i in 1:n_qubits]

f = [Array{Int}(undef, length(syndrome_values_past[k]), 2, 2, length(syndrome_values_future[k])) for k in 1:length(syndrome_values_future)]




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
                    #=
                    if [past_syndrome, (j-1) ⊻ PS[trellis_pos] , (k-1) ⊻ PS[n_qubits + trellis_pos], future_syndrome] in zero_syndrome_trellis[trellis_pos]
                        f[trellis_pos][i, j, k, l] = 1;
                    else
                        f[trellis_pos][i, j, k, l] = 0;
                    end
                    =#
                    if [past_syndrome, (j-1) , (k-1) , future_syndrome] in zero_syndrome_trellis[trellis_pos]
                        f[trellis_pos][i, j, k, l] = 1;
                    else
                        f[trellis_pos][i, j, k, l] = 0;
                    end
                end
            end
        end
    end
end

f_tensor = [ITensor(f[i], s[i], f_x[i], f_z[i], s[i+1]) for i in 1:n_qubits] 

parity_matrix_3 = Array{Float64}(undef, 2, 2, 2)
for ind_i in 1:2
    for ind_j in 1:2
        for ind_k in 1:2
            if (((ind_i - 1) + (ind_j - 1) + (ind_k - 1)) % 2) == 1 
                parity_matrix_3[ind_i, ind_j, ind_k] = 1
            else
                parity_matrix_3[ind_i, ind_j, ind_k] = 0
            end
        end
    end
end

error_inputs_x = [Index(2,"error_x_$i") for i in 1:n_qubits]
error_inputs_z = [Index(2,"error_z_$i") for i in 1:n_qubits]

#=
error_inputs = [ ]

for i in 1:(n_qubits)
    push!(error_inputs, error_inputs_x[i])
    push!(error_inputs, error_inputs_z[i])
end
=#

error_inputs = [Index(4, "error_q$i") for i in 1:n_qubits]

error_tranlation_matrix = Array{Float64}(undef, 4, 2, 2)
for ind_1 in 1:4
    for ind_2 in 1:2
        for ind_3 in 1:2
            if (ind_1-1) == ((ind_2-1) + 2*(ind_3-1))
                error_tranlation_matrix[ind_1, ind_2, ind_3] = 1
            else
                error_tranlation_matrix[ind_1, ind_2, ind_3] = 0
            end
        end
    end
end
error_translation_tensor = [ITensor(error_tranlation_matrix, error_inputs[i], error_inputs_x[i], error_inputs_z[i]) for i in 1:n_qubits]





Error_check_tensors_x = [ITensor(parity_matrix_3, f_x[i], error_inputs_x[i], x[i]) for i in 1:n_qubits]
Error_check_tensors_z = [ITensor(parity_matrix_3, f_z[i], error_inputs_z[i], z[i]) for i in 1:n_qubits]

# Fill in tensors for errors for testing and adjust contraction process
#=
error_matrix = [Array{Float64}(undef, 2) for i in 1:(2*n_qubits)]
for i in 1:(2*n_qubits)
    if Bool(PS[i])
        error_matrix[i][1], error_matrix[i][2] = 1, 0
    else
        error_matrix[i][1], error_matrix[i][2] = 0, 1
    end
end

Test_error_tensors_x = [ITensor(error_matrix[i], error_inputs_x[i]) for i in 1:n_qubits]
Test_error_tensors_z = [ITensor(error_matrix[i+n_qubits], error_inputs_z[i]) for i in 1:n_qubits]
=#

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
logical_x = [ITensor(indices_hadamard_x[i]) for i in 1:logical_weigth]
logical_z = [ITensor(indices_hadamard_z[i]) for i in 1:logical_weigth]

# Contraction steps
# initialize contracted Tensor
#f_contracted = f0


Logical_probabilities = Array{Float64}(undef, 2, 2)

log_counter_x = 0
log_counter_z = 0

trellis_MPS_logical_classes = [ ]

trellis_mps = nothing


println("Constructing Logical Decoder as MPS network")
for i in ProgressBar(1:length(Logical_probabilities))

    global trellis_mps = nothing

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
        #=
        if logical_operator_x[j] == 1
            global f_contracted = f_contracted*f_tensor[j]*Error_check_tensors_x[j]*Error_check_tensors_z[j]*qubit_x[j]* qubit_z[j]* P[j] * Hadamard_tensors_x[log_counter_x] * logical_x[log_counter_x]
        else
            global f_contracted = f_contracted*f_tensor[j]*qubit_x[j]* qubit_z[j]* P[j]
        end

        if logical_operator_z[j + n_qubits] == 1
            global f_contracted = f_contracted *  Hadamard_tensors_z[log_counter_z] * logical_z[log_counter_z]
        end
        =#
        
        #local_MPS_copy = trellis_mps
        f_contracted = nothing
        #=
        if j == 1
            f_contracted = f0
        else
            #=
            f_contracted = trellis_mps[1]
            for k in 2:length(trellis_mps)
                f_contracted = f_contracted * trellis_mps[k]
            end
            =#
            f_contracted = last(trellis_mps)
        end
        =#

        if logical_operator_x[j] == 1
            f_contracted = f_tensor[j]*Error_check_tensors_x[j]*Error_check_tensors_z[j]*error_translation_tensor[j]*qubit_x[j]* qubit_z[j]* P[j] * Hadamard_tensors_x[log_counter_x] * logical_x[log_counter_x]
        else
            f_contracted = f_tensor[j]*qubit_x[j]* qubit_z[j]* P[j]
        end

        if logical_operator_z[j + n_qubits] == 1
            f_contracted = f_contracted *  Hadamard_tensors_z[log_counter_z] * logical_z[log_counter_z]
        end
        #=
        if j != 1
            len = length(trellis_mps)
            for k in 1:(len - 1)
                f_contracted = trellis_mps[len - k] * f_contracted
            end
        end
        =#

        if j == n_qubits
            f_contracted = f_contracted * flast
        end

        if j == 1
            f_contracted = f0 * f_contracted
        end

        
        # TODO: Fix ordering of mps tensors
        if j != 1
            
            local_tensor = last(trellis_mps) * f_contracted
            trellis_mps = trellis_mps[1:(length(trellis_mps) - 1)]

            # Reorder indices
            indis = [ind for ind in inds(local_tensor)]
            for l in 1:length(indis)
                if indis[l] == s[j + 1]
                    ind = indis[l]
                    deleteat!(indis, l)
                    push!(indis, ind)
                end
            end

            local_tensor = permute(local_tensor, indis)
            indis = nothing

            local_mps = MPS(local_tensor, inds(local_tensor); cutoff= local_cutoff_svd, maxdim = local_truncate_dim)

            local_mps[2] = local_mps[1] * local_mps[2]

            #trellis_mps[length(trellis_mps)] = local_mps[1]
            local_mps = local_mps[2:length(local_mps)]

            trellis_mps_tensors = [trellis_mps[i] for i in 1:(length(trellis_mps))]

            local_mps_tensors = [local_mps[i] for i in 1:(length(local_mps))]

            global trellis_mps = MPS(vcat(trellis_mps_tensors, local_mps_tensors))
            trellis_mps_tensors = nothing
            local_mps_tensors = nothing
            local_mps = nothing
            local_tensor = nothing
        end
        if j == 1
            global trellis_mps = MPS(f_contracted, vcat(error_inputs[j], s[j+1]); cutoff=local_cutoff_svd, maxdim = local_truncate_dim)
        end

        #println(trellis_mps)

        #println(trellis_mps)
        #=
        if j == n_qubits
            global trellis_mps = MPS(f_contracted, vcat(error_inputs[1:2*j]); cutoff=5e-4, maxdim = 50)
        end
        =#
        f_contracted = nothing
            #println(f_contracted)
        #Logical_probabilities[i] = f_contracted[1]
    end
    #global f_contracted = f_contracted * flast
    #Logical_probabilities[i] = f_contracted[1]

    truncate!(trellis_mps; cutoff = cutoff_svd, maxdim = truncate_dim)
    push!(trellis_MPS_logical_classes, trellis_mps)

end



using Serialization

# Delete the file
if isfile(directory_name * "/logical_decoder_mps.jls")
    rm(directory_name * "/logical_decoder_mps.jls")
    rm(directory_name * "/error_input_edges.jls")
end


open(directory_name * "/logical_decoder_mps.jls", "w") do io
    serialize(io, trellis_MPS_logical_classes)
end

open(directory_name * "/error_input_edges.jls", "w") do io
    serialize(io, error_inputs)
end

open(directory_name * "/error_probability.jls", "w") do io
    serialize(io, p_qubit)
end