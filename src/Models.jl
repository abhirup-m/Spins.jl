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
        globalField::Float64=0.,
    )
    hamiltonian = tuple{string, vector{int64}, float64}[]
    for i in 1:numSpinsHalf
        for j in numSpinsHalf+1:2*numSpinsHalf
            append!(hamiltonian, [("zz", [i, j], 0.25)])
            append!(hamiltonian, [("+-", [i, j], 0.5), ("+-", [j, i], 0.5)])
        end
    end
    if globalField ≠ 0
        for i in 1:2*numSpinsHalf
            append!(hamiltonian, [("z", [i], globalField/2)])
        end
    end
    return hamiltonian
end


function LiebMattis(
        numSpinsHalf::Int64,
        intraCoupling::Float64;
        globalField::Float64=0.,
    )
    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]
    for i in 1:numSpinsHalf
        for j in numSpinsHalf+1:2*numSpinsHalf
            append!(hamiltonian, [("zz", [i, j], 0.25)])
            append!(hamiltonian, [("+-", [i, j], 0.5), ("+-", [j, i], 0.5)])
        end
    end
    for i in 1:numSpinsHalf
        for j in 1:numSpinsHalf
            append!(hamiltonian, [("zz", [i, j], 0.25 * intraCoupling)])
            append!(hamiltonian, [("+-", [i, j], 0.5 * intraCoupling), ("+-", [j, i], 0.5 * intraCoupling)])
            append!(hamiltonian, [("zz", [i + numSpinsHalf, j + numSpinsHalf], 0.25 * intraCoupling)])
            append!(hamiltonian, [("+-", [i + numSpinsHalf, j + numSpinsHalf], 0.5 * intraCoupling), ("+-", [j + numSpinsHalf, i + numSpinsHalf], 0.5 * intraCoupling)])
        end
    end
    if globalField ≠ 0
        for i in 1:2*numSpinsHalf
            append!(hamiltonian, [("z", [i], globalField/2)])
        end
    end
    return hamiltonian
end
