#=
This code runs the decoding procedure for the TN decoder. The decoder takes as arguments the
code specifics (mode, distance, codetype) and the syndrome
=#


function MPS_decoder(mode::String, d::Int, code_type::String, permuted::Bool, exact_decoding::Bool, s_meas)

    @assert d % 2 == true 
    @assert mode == "X" || mode == "Full" 
    
    directory_name = code_type
    
    if mode == "Full"
        directory_name = directory_name * "_d" * string(d)
    end
    if mode == "X"
        directory_name = directory_name * "_d" * string(d) * "_X"
    end
    
    # Open and read in the stabilizers
    filename = directory_name * "/stabilizers.txt"
    file_content = open(filename) do file
        read(file, String)
    end
    stabilizers = eval(Meta.parse(file_content))
    
    stabilizers_z = [ ]
    
    # Define number of qubits and stabilizers
    n_stabs = length(stabilizers)
    
    if mode == "X"
        n_stabs = Int(n_stabs / 2)
    end
    
    
    if mode == "X"
        stabilizers_z = stabilizers[(n_stabs+1):(2*n_stabs)]
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
    
    if !(exact_decoding)
        directory_name = directory_name * "/MPS_dim_" * string(truncate_dim)
    end
    
    error_inputs = open(directory_name * "/error_input_edges.jls", "r") do io
        deserialize(io)
    end
    
    
    
    n_qubits = Int(length(logical_operator_x)/2)
    
    # Load trellis MPS TN
    
    trellis_MPS_logical_classes = open(directory_name * "/logical_decoder_mps.jls", "r") do io
        deserialize(io)
    end
    
    
    #trellis_MPS_logical_classes = nothing
    Hadamard_matrix = nothing
    PS = nothing
    
    p = open(directory_name * "/error_probability.jls", "r") do io
        deserialize(io)
    end
    
    if mode == "Full"
        p /= 3.
    end
    
    #println(occured_error)

    # Determine weight of the error


    #s_meas = check_orthogonality(stabilizers, occured_error)

    Hadamard_matrix = Array{Float64}(undef, 2, 2)
    Hadamard_matrix[1,1], Hadamard_matrix[1,2] = 1, 1
    Hadamard_matrix[2,1], Hadamard_matrix[2,2] = 1, -1

    Logical_probabilities = Array{Float64}(undef, 2, 2)

    #PS = find_Pauli_give_syndrome(stabilizers, s_meas, n_qubits)
    PS = zeros(Int, 2*n_qubits)

    for i in 1:n_stabs
        if Bool(s_meas[i])
            PS = (PS .+ destabilizers[i]) .% 2
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

    if mode == "X"
        for i in 1:n_qubits
            error_matrix[i][1] = 0
            error_matrix[i][3] = 0
            if PS[i+n_qubits]==1
                error_matrix[i][4] = 1
                error_matrix[i][2] = 0
            else
                error_matrix[i][4] = 0
                error_matrix[i][2] = 1
            end
        end
    end
    #println(error_matrix)

    Test_error_tensors = [ITensor(error_matrix[i], error_inputs[i]) for i in 1:(n_qubits)]





    #println("Physical Error")
    #println(occured_error)

    #println("Detected syndrome")
    #println(s_meas)

    for i in 1:length(Logical_probabilities)
        for j in 1:(n_qubits)
            trellis_MPS_logical_classes[i][j] = trellis_MPS_logical_classes[i][j] * Test_error_tensors[j]
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

    #println(Logical_probabilities)
    # Hadamard transform of data
    Logical_probabilities = Hadamard_matrix * Logical_probabilities * Hadamard_matrix

    # Logical error of maximum probability
    index_max = argmax(Logical_probabilities)

    correcting_operation = zeros(Int, n_qubits*2)

    logical_class = 0

    #println("Correction operation:")
    if index_max == CartesianIndex(1, 1)
        correcting_operation = PS
        logical_class = 1
        #println(PS)
    end
    if index_max == CartesianIndex(2, 1)
        correcting_operation = (PS .+ logical_operator_x) .% 2 
        logical_class = 2
        #println((PS .+ logical_operator_x) .% 2 )
    end
    if index_max == CartesianIndex(1, 2)
        correcting_operation = (PS .+ logical_operator_z) .% 2 
        logical_class = 3
        #println((PS .+ logical_operator_z) .% 2 )
    end
    if index_max == CartesianIndex(2, 2)
        correcting_operation = (PS .+ logical_operator_z .+ logical_operator_x) .% 2 
        logical_class = 4
        #println((PS .+ logical_operator_x .+ logical_operator_z) .% 2 )
    end
    return correcting_operation, p, logical_class, Logical_probabilities
end


function MPS_decoder_opposite_contraction(mode::String, d::Int, code_type::String, permuted::Bool, exact_decoding::Bool, s_meas)

    @assert d % 2 == true 
    @assert mode == "X" || mode == "Full" 
    
    directory_name = code_type
    
    if mode == "Full"
        directory_name = directory_name * "_d" * string(d)
    end
    if mode == "X"
        directory_name = directory_name * "_d" * string(d) * "_X"
    end
    
    # Open and read in the stabilizers
    filename = directory_name * "/stabilizers.txt"
    file_content = open(filename) do file
        read(file, String)
    end
    stabilizers = eval(Meta.parse(file_content))
    
    stabilizers_z = [ ]
    
    # Define number of qubits and stabilizers
    n_stabs = length(stabilizers)
    
    if mode == "X"
        n_stabs = Int(n_stabs / 2)
    end
    
    
    if mode == "X"
        stabilizers_z = stabilizers[(n_stabs+1):(2*n_stabs)]
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
    
    if !(exact_decoding)
        directory_name = directory_name * "/MPS_dim_" * string(truncate_dim)
    end
    
    error_inputs = open(directory_name * "/error_input_edges.jls", "r") do io
        deserialize(io)
    end
    
    
    
    n_qubits = Int(length(logical_operator_x)/2)
    
    # Load trellis MPS TN
    
    trellis_MPS_logical_classes = open(directory_name * "/logical_decoder_mps.jls", "r") do io
        deserialize(io)
    end
    
    
    #trellis_MPS_logical_classes = nothing
    Hadamard_matrix = nothing
    PS = nothing
    
    p = open(directory_name * "/error_probability.jls", "r") do io
        deserialize(io)
    end
    
    if mode == "Full"
        p /= 3.
    end
    
    #println(occured_error)

    # Determine weight of the error


    #s_meas = check_orthogonality(stabilizers, occured_error)

    Hadamard_matrix = Array{Float64}(undef, 2, 2)
    Hadamard_matrix[1,1], Hadamard_matrix[1,2] = 1, 1
    Hadamard_matrix[2,1], Hadamard_matrix[2,2] = 1, -1

    Logical_probabilities = Array{Float64}(undef, 2, 2)

    #PS = find_Pauli_give_syndrome(stabilizers, s_meas, n_qubits)
    PS = zeros(Int, 2*n_qubits)

    for i in 1:n_stabs
        if Bool(s_meas[i])
            PS = (PS .+ destabilizers[i]) .% 2
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

    if mode == "X"
        for i in 1:n_qubits
            error_matrix[i][1] = 0
            error_matrix[i][3] = 0
            if PS[i+n_qubits]==1
                error_matrix[i][4] = 1
                error_matrix[i][2] = 0
            else
                error_matrix[i][4] = 0
                error_matrix[i][2] = 1
            end
        end
    end
    #println(error_matrix)

    Test_error_tensors = [ITensor(error_matrix[i], error_inputs[i]) for i in 1:(n_qubits)]





    #println("Physical Error")
    #println(occured_error)

    #println("Detected syndrome")
    #println(s_meas)

    for i in 1:length(Logical_probabilities)
        for j in n_qubits:-1:1
            trellis_MPS_logical_classes[i][j] = trellis_MPS_logical_classes[i][j] * Test_error_tensors[j]
        end
    end

    #println("Decoding logical probabilities")
    for j in 1:length(Logical_probabilities)
        V = ITensor(1.)
        for i in n_qubits:-1:1
            V *= (trellis_MPS_logical_classes[j][i])
        end
        Logical_probabilities[j] = scalar(V)
    end 

    #println(Logical_probabilities)
    # Hadamard transform of data
    Logical_probabilities = Hadamard_matrix * Logical_probabilities * Hadamard_matrix

    # Logical error of maximum probability
    index_max = argmax(Logical_probabilities)

    correcting_operation = zeros(Int, n_qubits*2)

    logical_class = 0

    #println("Correction operation:")
    if index_max == CartesianIndex(1, 1)
        correcting_operation = PS
        logical_class = 1
        #println(PS)
    end
    if index_max == CartesianIndex(2, 1)
        correcting_operation = (PS .+ logical_operator_x) .% 2 
        logical_class = 2
        #println((PS .+ logical_operator_x) .% 2 )
    end
    if index_max == CartesianIndex(1, 2)
        correcting_operation = (PS .+ logical_operator_z) .% 2 
        logical_class = 3
        #println((PS .+ logical_operator_z) .% 2 )
    end
    if index_max == CartesianIndex(2, 2)
        correcting_operation = (PS .+ logical_operator_z .+ logical_operator_x) .% 2 
        logical_class = 4
        #println((PS .+ logical_operator_x .+ logical_operator_z) .% 2 )
    end
    return correcting_operation, p, logical_class, Logical_probabilities
end