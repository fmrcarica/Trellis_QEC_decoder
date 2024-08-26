function check_orthogonality(vectors::Vector{Vector{Int64}}, vector::Vector{Int64})
    n = size(vectors, 1)
    orthogonality = Vector{Bool}(undef, n)
    
    for i in 1:n
        orthogonality[i] = symplectic_product(vectors[i], vector,n_qubits)
    end
    
    return orthogonality
end

function symplectic_product(v1, v2, n)
    # The symplectic product of two binary vectors v1 and v2
    return sum(v1[i] * v2[n+i] + v1[n+i] * v2[i] for i in 1:n) % 2
end

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

function binary_vector_to_reversed_int(binary_vector)::Int64
    # Reverse the binary vector
    
    # Convert the reversed binary vector to an integer
    result = 0
    for (i, bit) in enumerate(binary_vector)
        result += bit * (2^(i-1))
    end
    
    return result
end