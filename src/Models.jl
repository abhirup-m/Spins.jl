"""
    Heisenberg(numSpins; joinEnds=true, globalField=0.)

Hamiltonian definition for a spin-half nearest-neighbour 
1D Heisenberg model.
PBC is enabled by default, can be 
disabled by passing `joinEnds=false`.
"""
function Heisenberg(
        numSpins::Int64;
        joinEnds::Bool=true,
        globalField::Float64=0.,
    )
    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]
    for i in 1:numSpins
        if i < numSpins
            append!(hamiltonian, [("zz", [i, i+1], 0.25)])
            append!(hamiltonian, [("+-", [i, i+1], 0.5), ("+-", [i+1, i], 0.5)])
        elseif joinEnds && numSpins > 2
            append!(hamiltonian, [("zz", [1, numSpins], 0.25)])
            append!(hamiltonian, [("+-", [1, numSpins], 0.5), ("+-", [numSpins, 1], 0.5)])
        end
        if globalField ≠ 0
            append!(hamiltonian, [("z", [i], -globalField)])
        end
    end
    return hamiltonian
end


"""
    LiebMattis(numSpinsHalf; globalField=0.)

Hamiltonian definition for a spin-half Lieb-Mattis model.
where A_i and B_j are the ith and jth spins of the A and B
sublattices. The total set of spin indices runs from
1 to 2 * numSpinsHalf, among which the first 1 ... numSpinsHalf
represent the A sublattice, and the other indices 2... numSpinsHalf
represents the other sublattice.
"""
function LiebMattis(
        numSpinsHalf::Int64;
        totSpinSqCoupling::Number,
        globalField::Number=0.,
    )
    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]
    for i in 1:2*numSpinsHalf
        for j in 1:2*numSpinsHalf
            if i ≤ numSpinsHalf && j ≥ numSpinsHalf + 1
                append!(hamiltonian, [("zz", [i, j], (1 - totSpinSqCoupling ) * 0.25)])
                append!(hamiltonian, [("+-", [i, j], (1 - totSpinSqCoupling ) * 0.5), 
                                      ("+-", [j, i], (1 - totSpinSqCoupling ) * 0.5)
                                     ])
            end
            if totSpinSqCoupling ≠ 0
                append!(hamiltonian, [("zz", [i, j], 0.5 * totSpinSqCoupling * 0.25)])
                append!(hamiltonian, [("+-", [i, j], 0.5 * totSpinSqCoupling * 0.5), 
                                      ("+-", [j, i], 0.5 * totSpinSqCoupling * 0.5)
                                     ])
            end
        end
        if globalField ≠ 0
            append!(hamiltonian, [("z", [i], globalField/2)])
        end
    end
    return hamiltonian
end
