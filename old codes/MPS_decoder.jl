
using ITensors
using LinearAlgebra
using AbstractAlgebra
using Serialization
using ProgressBars

include("helper_functions.jl")

path = "C:/Users/filip/Documents/ETH/Quantum Engineering MSc/Semester Project - QIT/"

mode = "Full" #"Full" or ""X"
d = 7 # Only works for d unveven

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

# Open and read in the destabilizers
filename = directory_name * "/destabilizers.txt"
file_content = open(filename) do file
    read(file, String)
end
destabilizers = eval(Meta.parse(file_content))

# Open and read in the logicals
filename = directory_name * "/logicals.txt"
file_content = open(filename) do file
    read(file, String)
end
logical_operators = eval(Meta.parse(file_content))

logical_operator_x = logical_operators[1]
logical_operator_z = logical_operators[2]

error_inputs = open(directory_name * "/error_input_edges.jls", "r") do io
    deserialize(io)
end




# Load trellis MPS TN
trellis_MPS_logical_classes = open(directory_name * "/logical_decoder_mps.jls", "r") do io
    deserialize(io)
end

p = open(directory_name * "/error_probability.jls", "r") do io
    deserialize(io)
end

if mode == "Full"
    p /= 3.
end


#trellis_MPS_logical_classes = deepcopy(original_trellis_MPS_logical_classes)

occured_error = Pauli_error(p, n_qubits)

# Determine weight of the error

s_meas = check_orthogonality(stabilizers, occured_error)

Hadamard_matrix = Array{Float64}(undef, 2, 2)
Hadamard_matrix[1,1], Hadamard_matrix[1,2] = 1, 1
Hadamard_matrix[2,1], Hadamard_matrix[2,2] = 1, -1

Logical_probabilities = Array{Float64}(undef, 2, 2)

#PS = find_Pauli_give_syndrome(stabilizers, s_meas, n_qubits)
PS = zeros(Int, 2*n_qubits)

for i in 1:n_stabs
    if s_meas[i]
        global PS = (PS .+ destabilizers[i]) .% 2
    end
end

#=
P_S = Vector{Int}(undef, 2*n_qubits)
for i in 1:n_qubits
    P_S[2*i-1] = PS[i]
    P_S[2*i] = PS[i+n_qubits]
end
=#

error_matrix = [Array{Float64}(undef, 4) for i in 1:n_qubits]

for i in 1:n_qubits
    for j in 1:4
        if (PS[i] + 2 * PS[i + n_qubits]) == j - 1
            error_matrix[i][j] = 1
        else
            error_matrix[i][j] = 0
        end
    end
end

Test_error_tensors = [ITensor(error_matrix[i], error_inputs[i]) for i in 1:(n_qubits)]





#println("Physical Error")
#println(occured_error)

#println("Detected syndrome")
#println(s_meas)


for i in 1:length(Logical_probabilities)
    for j in 1:(n_qubits)
        global trellis_MPS_logical_classes[i][j] = trellis_MPS_logical_classes[i][j] * Test_error_tensors[j]
    end
end

#println("Decoding logical probabilities")
for j in 1:length(Logical_probabilities)
    V = ITensor(1.)
    for i in 1:(n_qubits)
        V *= (trellis_MPS_logical_classes[j][i])
    end
    Logical_probabilities[j] = scalar(V)
end 


# Hadamard transform of data
Logical_probabilities = Hadamard_matrix * Logical_probabilities * Hadamard_matrix

# Logical error of maximum probability
index_max = argmax(Logical_probabilities)

correcting_operation = zeros(Int, n_qubits*2)

#println("Most likely error to have occured")
if index_max == CartesianIndex(1, 1)
    #println("No logical error")
end
if index_max == CartesianIndex(2, 1)
    #println("Logical x error")
end
if index_max == CartesianIndex(1, 2)
    #println("Logical z error")
end
if index_max == CartesianIndex(2, 2)
    #println("Logical y error")
end

#println("Correction operation:")
if index_max == CartesianIndex(1, 1)
    global correcting_operation = PS
    #println(PS)
end
if index_max == CartesianIndex(2, 1)
    global correcting_operation = (PS .+ logical_operator_x) .% 2 
    #println((PS .+ logical_operator_x) .% 2 )
end
if index_max == CartesianIndex(1, 2)
    global correcting_operation = (PS .+ logical_operator_z) .% 2 
    #println((PS .+ logical_operator_z) .% 2 )
end
if index_max == CartesianIndex(2, 2)
    global correcting_operation = (PS .+ logical_operator_z .+ logical_operator_x) .% 2 
    #println((PS .+ logical_operator_x .+ logical_operator_z) .% 2 )
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
#=
if in_span
    global n_true += 1
end
=#