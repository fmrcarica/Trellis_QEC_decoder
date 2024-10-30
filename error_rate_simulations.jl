using ITensors
using LinearAlgebra
using AbstractAlgebra
using Serialization
using ProgressBars

include("helper_functions.jl")
include("trellis_mps_constructor_function.jl")
include("MPS_decoder_function.jl")

mode = "Full" #"Full" or ""X"
#d = 3 # Only works for distance d unveven
code_type = "surface" # "surface" or "color666"
exact_decoding = true # bool
permuted = false # bool

#@assert d % 2 == true 
@assert mode == "X" || mode == "Full" 


# Error simulation for exact Decoding

# Define physical error rates
error_probabilities = [0.06]#range(0.06, stop=0.15, length=11)

# Number of shots
n_shots = 10000

distances = [3,5,7,9,11]


for d in distances
    logical_error_rates = []
    file_name = "sim_err_rate" * "_d" * string(d) * ".txt"

    if isfile(file_name)
        # File exists, clear its contents
        open(file_name, "w") do file
            # This automatically truncates the file, so no need to write anything
        end
    else
        # File does not exist, create it
        open(file_name, "w") do file
        end
    end

    for p in error_probabilities
        # Create trellis MPS
        println("Create TN")
        n_qubits, stabilizers, stabilizers_z, p_qubit  = create_trellis_MPS(mode, d, code_type, permuted, exact_decoding, p)
        n_true = 0
        logical_error_rate = 0
        p_q = p
        if mode == "Full"
            p_q = p/3
        end
        
        for i in 1:n_shots
            # Simulate error under depolarizing channel
            occured_error = Pauli_error(p_q, n_qubits)
            # Determine syndrome
            s_meas = check_orthogonality(stabilizers, occured_error)
            
            # Decode
            correcting_operation, p_qubit, log_class, prob_mat = MPS_decoder(mode, d, code_type, permuted, exact_decoding, s_meas)
            #correcting_operation_inv, p_qubit_inv, log_class_inv, prob_mat_inv = MPS_decoder_opposite_contraction(mode, d, code_type, permuted, exact_decoding, s_meas)

            # Determine if correction was succesfull by seeing if correction * error is in the stabilizer group
            vectors = transpose(reduce(hcat, stabilizers))
            if mode == "X"
                vectors = transpose(reduce(hcat, stabilizers_z))
            end
            test_op = (correcting_operation .+ occured_error) .% 2
            
            # Check if test_op is in the span
            in_span = is_in_span(test_op, vectors)
            
            if in_span
                n_true += 1
            end
            
            
        end
        # Determine logical error rate
        logical_error_rate = 1.0 - n_true / n_shots
        push!(logical_error_rates, logical_error_rate)

        
        println(logical_error_rate)
        

        # Deocoder errors

    end
    open(file_name, "w") do file
        # Write p values as the first row, separated by commas
        write(file, join(error_probabilities, ",") * "\n")
        
        # Write logical_error_rate values as the second row, separated by commas
        write(file, join(logical_error_rates, ",") * "\n")
    end
end