using LinearAlgebra

function Spectrum(
        hamiltonian::Dict{Tuple{String,Vector{Int64}},Float64},
        basisStates::Dict{Float64,Vector{BitArray}}, 
    )
    hamiltonianMatrix = generalOperatorMatrix(basisStates, hamiltonian)
    eigenStates = Dict{Float64, Vector{Vector{Float64}}}()
    eigenValues = Dict{Float64, Vector{Float64}}()
    for (sector, matrix) in hamiltonianMatrix
        V, F = eigen(Hermitian(matrix))
        eigenValues[sector] = V
        eigenStates[sector] = [column for column in eachcol(F)]
    end
    return eigenValues, eigenStates
end
