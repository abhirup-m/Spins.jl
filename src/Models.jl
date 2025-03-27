function Heisenberg(
        numSpins::Int64;
        joinEnds::Bool=true
    )
    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]
    for i in 1:numSpins-1
        append!(hamiltonian, [("zz", [i, i+1], 0.25), ("zz", [i+1, i], 0.25)])
        append!(hamiltonian, [("+-", [i, i+1], 0.5), ("+-", [i+1, i], 0.5)])
    end
    if joinEnds
        append!(hamiltonian, [("zz", [1, numSpins], 0.25), ("zz", [numSpins, 1], 0.25)])
        append!(hamiltonian, [("+-", [1, numSpins], 0.5), ("+-", [numSpins, 1], 0.5)])
    end
    return hamiltonian
end


function LiebMattis(
        numSpinsHalf::Int64;
    )
    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]
    for i in 1:2:2*numSpinsHalf
        for j in 2:2:2*numSpinsHalf
            append!(hamiltonian, [("zz", [i, j], 0.25), ("zz", [j, i], 0.25)])
            append!(hamiltonian, [("+-", [i, j], 0.5), ("+-", [j, i], 0.5)])
        end
    end
    return hamiltonian
end
